#! /usr/bin/env python
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pprint import pprint
from typing import Optional, List

import json
import time
import re
from dataclasses import dataclass, asdict
from typing import List, Optional, Dict, Any
import requests
import sys
import argparse
from urllib.parse import quote_plus


import django

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings')
django.setup()

from sys import argv

# ---------- Config ----------
EUROPE_PMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
# Simple user-agent to avoid being rude
USER_AGENT = "env-resistance-finder/1.0 (+https://example.invalid/; contact: none)"

# Keywords for environmental contexts and data deposition markers
ENV_KEYWORDS = [
    "environment", "wastewater", "sewage", "river", "soil", "marine", "coastal",
    "surface water", "groundwater", "sludge", "sediment", "hospital effluent"
]
DEPOSITION_MARKERS = [
    "sequence read archive", "sra", "ena", "ebi", "bioproject", "prjna", "prjeb",
    "accession", "biosample", "run accession", "raw reads", "short read archive"
]
# simple regex to capture PRJ* style accessions e.g. PRJNA12345 PRJEB12345
PRJ_REGEX = re.compile(r"\b(PRJNA|PRJEB|PRJDB|PRJDA|PRJCA)\w*\d+\b", re.IGNORECASE)



# ---------- Data classes ----------

@dataclass
class PaperMatch:
    pmid: Optional[str]
    pmcid: Optional[str]
    title: str
    author_string: Optional[str]
    journal: Optional[str]
    pub_year: Optional[int]
    abstract: Optional[str]
    env_keywords_found: List[str]
    deposition_markers_found: List[str]
    prj_accessions: List[str]
    score: float  # heuristic score for relevance


# ---------- Step 2: search literature (Europe PMC) ----------
def build_europepmc_query(search_terms: List[str], env_terms: List[str], page_size: int = 25) -> str:
    """
    Build a Europe PMC search query combining protein/gene terms and environmental context.
    We'll look for papers that mention at least one protein term AND at least one env term.
    """
    prot_q = " AND ".join([f'"{t}"' for t in search_terms])
    env_q = " OR ".join([f'"{t}"' for t in env_terms])
    # require mention of sequencing deposition markers optionally by adding them as SHOULD terms
    dep_markers_q = " OR ".join([f'"{m}"' for m in DEPOSITION_MARKERS])
    # Combined: (prot) AND (env) AND (deposition markers)
    query = f"({prot_q}) AND ({env_q}) AND ({dep_markers_q})"

    return query

def search_europe_pmc(query, result_type='core', format='json', page_size=25):
    """
    Searches Europe PMC for articles based on a given query.

    MODIFIED:
    - Defaults to result_type='core' to ensure abstracts are fetched.
    - Filters the JSON response to return only results with a non-empty abstract.

    Args:
        query (str): The search query string.
        result_type (str): The type of results to retrieve (e.g., 'lite', 'core', 'full').
                           **Must be 'core' or 'full' for abstract filtering to work.**
        format (str): The desired output format ('json' or 'xml').
        page_size (int): The maximum number of results to return.

    Returns:
        dict or str: The parsed JSON response (with filtered results) or raw XML string.
    """
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {
        "query": query,
        "resultType": result_type,
        "format": format,
        "pageSize": page_size
    }

    # Add a warning if the user explicitly requests 'lite' with 'json'
    if result_type == 'lite' and format == 'json':
        print("Warning: result_type='lite' does not include abstracts. "
              "Filtering will likely result in an empty list. "
              "Consider using result_type='core' or 'full'.")

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raise an exception for HTTP errors

        if format == 'json':
            data = response.json()

            # --- START: Modification to filter results ---

            # Safely get the list of results, default to empty list if not found
            original_results = data.get('resultList', {}).get('result', [])

            # Use a list comprehension to filter for non-empty abstracts
            # .get('abstractText', '') handles missing keys safely
            # .strip() handles abstracts that contain only whitespace
            filtered_results = [
                article for article in original_results
                if article.get('abstractText', '').strip()
            ]

            # Replace the original result list with the filtered one
            if 'resultList' in data:
                data['resultList']['result'] = filtered_results
            # If 'resultList' was not in the original data, we just return the data as-is

            # --- END: Modification ---

            return data
        else:
            # Return raw XML for XML format (no filtering applied)
            return response.text

    except requests.exceptions.RequestException as e:
        print(f"Error making API request: {e}")
        return None


def extract_papermatch_from_record(rec: dict) -> PaperMatch:
    """
    Convert an Europe PMC record dict into a PaperMatch with heuristics to find env and deposition markers.
    """
    abstract = rec.get("abstractText") or rec.get("abstract") or ""
    title = rec.get("title") or ""
    author_string = rec.get("authorString")
    journal = rec.get("journalTitle")
    pub_year = None
    try:
        pub_year = int(rec.get("pubYear")) if rec.get("pubYear") else None
    except Exception:
        pub_year = None
    pmid = rec.get("pmid")
    pmcid = rec.get("pmcid")
    # find environment keywords
    env_found = []
    for kw in ENV_KEYWORDS:
        if kw.lower() in (abstract or "").lower() or kw.lower() in title.lower():
            env_found.append(kw)
    # deposition markers
    dep_found = []
    for marker in DEPOSITION_MARKERS:
        if marker.lower() in (abstract or "").lower() or marker.lower() in title.lower():
            dep_found.append(marker)
    # PRJ accessions
    prjs = PRJ_REGEX.findall((abstract or ""))
    # PRJ_REGEX.findall returns tuples for groups; flatten
    prj_accessions = []
    for m in PRJ_REGEX.finditer((abstract or "")):
        prj_accessions.append(m.group(0))
    # heuristic score: prefer records with both env and deposition markers
    score = 0.0
    score += 1.5 * len(env_found)
    score += 2.0 * len(dep_found)
    score += 3.0 * (1 if prj_accessions else 0)
    return PaperMatch(
        pmid=pmid,
        pmcid=pmcid,
        title=title,
        author_string=author_string,
        journal=journal,
        pub_year=pub_year,
        abstract=abstract,
        env_keywords_found=env_found,
        deposition_markers_found=dep_found,
        prj_accessions=prj_accessions,
        score=score
    )



def find_candidates(terms, page_size: int = 25) -> List[PaperMatch]:


    qry = build_europepmc_query(terms, ENV_KEYWORDS)
    print(qry)

    results = search_europe_pmc(qry, page_size=page_size)
    if not results:
        print(f"Warning: Europe PMC query failed for terms {terms}", file=sys.stderr)
        return {}

    hits = results.get("resultList", {}).get("result", []) if results else []

    paper_matches: List[PaperMatch] = []
    for rec in hits:
        pm = extract_papermatch_from_record(rec)
        if pm.score <= 0: continue

        paper_matches.append(pm)
    # sort by heuristics score descending
    paper_matches.sort(key=lambda x: x.score, reverse=True)

    # be polite
    # time.sleep(sleep_between_api)

    return paper_matches



def run():
    pass

def test_run():
    terms = ["gyrA", "S81L"]
    page_size = 100
    candidate_papers = find_candidates(terms, page_size)
    print(f"Number of candidate papers:", len(candidate_papers))

    for cp in candidate_papers[:5]:
        print(cp)
        print()



#######################
def main():
    if len(argv) < 2:
        exit("Tell me what to do - test or run?")
    if argv[1] == 'run':
        run()
    elif argv[1] == 'test':
        test_run()
    else:
        print(f"I don't know what is '{argv[1]}'. What should I do - test or run?")

#######################
if __name__ == "__main__":
    main()
