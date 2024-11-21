#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such
import requests
import pandas as pd

def run():
    # Define the base URL for the RCSB PDB Search API
    base_url = "https://search.rcsb.org/structure-search"

    # Construct your query
    query = {
        "query": {
            "type": "terminal",
            "service": "pdb",
            "parameters": {
                "gene": "rpsL",
                "ligand": "streptomycin"
            }
        }
    }

    # Make a request to the Search API
    response = requests.post(base_url, json=query)

    # Check if the request was successful
    if response.status_code == 200:
        results = response.json()
        # Extract PDB IDs from results
        pdb_ids = [entry['identifier'] for entry in results['result_set']]
    else:
        print("Error:", response.status_code)
        exit()

    # Convert to DataFrame for better visualization
    df = pd.DataFrame(pdb_ids, columns=["PDB ID"])
    print(df)

#######################
if __name__ == "__main__":
    run()
