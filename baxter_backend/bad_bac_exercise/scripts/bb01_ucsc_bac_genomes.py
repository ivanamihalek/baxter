#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb01_ucsc_bac_genomes
# in that case django will take care of the paths and also check for migrations and such
import re

import pandas as pd
import requests

from pathlib import Path
from bs4 import BeautifulSoup


def douwnload_html_source(url):
    # Fetch the HTML content from the URL
    response = requests.get(url)

    if response.status_code != 200:
        print(f"Failed to retrieve content from {url}, status code: {response.status_code}")
        return

    return response

#
def download_ucsc_bac_genomes_info(ucsc_genomes_file):
    blah = 'common\xa0name\xa0and view\xa0in\xa0UCSC\xa0browser[IGV\xa0browser]'
    print(repr(blah))
    blah2 = blah.replace('\u00A0', ' ')
    print(repr(blah2))
    exit()


    response = douwnload_html_source("https://hgdownload.soe.ucsc.edu/hubs/bacteria/index.html")

    # Parse the HTML content using BeautifulSoup
    cleaned_text = response.text.replace('\u00A0', ' ')
    soup = BeautifulSoup(cleaned_text, 'html.parser')

    # Find the first table in the HTML content
    table = soup.find('table')

    if table is None:
        print("No table found in the provided HTML.")
        return

    # Use pandas to read the HTML table into a DataFrame
    df = pd.read_html(str(table))[0]
    print(df.columns)
    df.rename(columns={'common name and view in UCSC browser[IGV browser]': 'common name',
                       'scientific name and data download': 'scientific name',
                       'assembly date, source link': 'assembly date',
                       }, inplace=True)
    for colname in df.columns:
        print(repr(colname))
    exit()

    # Save the DataFrame to a TSV file
    df.to_csv(ucsc_genomes_file, sep='\t', index=False)
    print(f"Table has been written to {ucsc_genomes_file} in TSV format.")


def run():
    ucsc_genomes_home = Path("/storage/databases/ucsc")
    ucsc_genomes_file = ucsc_genomes_home / "bac_genomes.tsv"
    if not ucsc_genomes_home.exists():
        print(f"{ucsc_genomes_home} not found")
        exit()

    if not ucsc_genomes_file.exists() or  ucsc_genomes_file.stat().st_size == 0:
        download_ucsc_bac_genomes_info(ucsc_genomes_file)


#######################
if __name__ == "__main__":
    run()
