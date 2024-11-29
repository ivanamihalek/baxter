#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb01_ucsc_bac_genomes
# in that case django will take care of the paths and also check for migrations and such
import re
from pprint import pprint

import pandas as pd
import requests

from pathlib import Path
from bs4 import BeautifulSoup

from bad_bac_exercise.models import UCSCAssembly
from django.db import connection


def douwnload_html_source(url):
    # Fetch the HTML content from the URL
    response = requests.get(url)

    if response.status_code != 200:
        print(f"Failed to retrieve content from {url}, status code: {response.status_code}")
        return

    return response


def download_ucsc_bac_genomes_info(ucsc_genomes_file):

    response = douwnload_html_source("https://hgdownload.soe.ucsc.edu/hubs/bacteria/index.html")

    # Parse the HTML content using BeautifulSoup
    cleaned_text = response.text.replace('&nbsp;', ' ')

    soup = BeautifulSoup(cleaned_text, 'html.parser')

    # Find the first table in the HTML content
    table = soup.find('table')

    if table is None:
        print("No table found in the provided HTML.")
        return

    # Use pandas to read the HTML table into a DataFrame
    df = pd.read_html(str(table))[0]

    # the original column names indicate the meaning of the hyperlinks they contain
    # - get rid of that
    df.rename(columns={'common name and view in UCSC browser[IGV browser]': 'common_name',
                       'scientific name and data download': 'scientific_name',
                       'assembly date, source link': 'assembly_date',
                       'common name': 'common_name',
                       'NCBI assembly': 'NCBI_assembly',
                       }, inplace=True)
    # get rid of the "count" column (it not a count of anything, it is a row number)
    del df['count']

    # get rid of the string that links to IGV in the original webpage
    df['common_name'] = df['common_name'].str.replace('[IGV]', '', regex=False)

    # Split 'NCBI assembly' into  assembly identifier and NCBI accession number
    split_cols = df['NCBI_assembly'].str.split('_', expand=True)
    df['RefSeq_assembly_id'] = split_cols[0] + '_' + split_cols[1]
    df['NCBI_accession_number'] = split_cols[2]
    del df['NCBI_assembly']

    # Save the DataFrame to a TSV file
    df.to_csv(ucsc_genomes_file, sep='\t', index=False)
    print(f"Wrote {ucsc_genomes_file} in TSV format.")


def run():

    # download the list of bacterial genomes from the UCSC website
    # skip if the table is already present
    ucsc_genomes_home = Path("/storage/databases/ucsc")
    ucsc_genomes_file = ucsc_genomes_home / "bac_genomes.tsv"
    if not ucsc_genomes_home.exists():
        print(f"{ucsc_genomes_home} not found")
        exit()

    if ucsc_genomes_file.exists() and ucsc_genomes_file.stat().st_size >= 0:
        print(f"found {ucsc_genomes_file}")
    else:
        download_ucsc_bac_genomes_info(ucsc_genomes_file)

    # read tsv into dataframe
    df = pd.read_csv(ucsc_genomes_file, sep="\t")
    # Change all column names to lowercase
    df.columns = df.columns.str.lower()
    # remove all rows that do not have a proper accession identifier
    df = df[df['ncbi_accession_number'].str.startswith('ASM')]
    # Convert date column to proper format
    df['assembly_date'] = pd.to_datetime(df['assembly_date']).dt.date
    # Convert DataFrame to dictionary format
    data_dict = df.to_dict(orient='records')

    # in the local database, clear the old data if present
    # - delete all records from the table
    UCSCAssembly.objects.all().delete()
    #
    # # Step 2: Reset the auto-increment index - if frequently done this should be handled by migration
    # but the hope here is that this is a one time only thing
    # note also that in PostgreSQL the symtax would be
    # "ALTER SEQUENCE myapp_mymodel_id_seq RESTART WITH 1"
    with connection.cursor() as cursor:
        cursor.execute("ALTER TABLE ucsc_assemblies AUTO_INCREMENT = 1")  # For MySQL
    # Bulk create instances in the database
    UCSCAssembly.objects.bulk_create([UCSCAssembly(**record) for record in data_dict])


#######################
if __name__ == "__main__":
    run()
