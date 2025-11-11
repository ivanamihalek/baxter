#! /usr/bin/env python
from __future__ import annotations
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

import os
import re
import sys
from time import sleep

import requests
from sys import argv
from dataclasses import dataclass
from typing import List, Optional, Iterable, Mapping, Dict, Tuple
from bs4 import BeautifulSoup
from lxml import etree
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from zeep import Client
from zeep.exceptions import Fault, TransportError
from zeep.transports import Transport

import django

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings')
django.setup()
from models.bad_bac_models import AntibioticResMutation


# ---------- Config ----------
EUROPE_PMC_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EUROPE_PMC_SEARCH_URL = f"{EUROPE_PMC_BASE}/search"

# Simple user-agent to avoid being rude
USER_AGENT = "env-resistance-finder/1.0 (+https://example.invalid/; contact: none)"

# Keywords for environmental contexts and data deposition markers
ENV_KEYWORDS = [
    "environment", "wastewater", "sewage", "river", "soil", "marine", "coastal",
    "surface water", "groundwater", "sludge", "sediment", "hospital effluent",
    "aquatic"
]
DEPOSITION_MARKERS = [
    "sequence read archive", "sra", "ena", "ebi", "bioproject", "prjna", "prjeb",
    "accession", "biosample", "run accession", "raw reads", "short read archive"
]
# simple regex to capture PRJ* style accessions e.g. PRJNA12345 PRJEB12345
PRJ_REGEX = re.compile(r"\b(PRJNA|PRJEB|PRJDB|PRJDA|PRJCA)\w*\d+\b", re.IGNORECASE)


BLOCK_TAGS = {
"address","article","aside","blockquote","canvas","dd","div","dl","dt","fieldset",
"figcaption","figure","footer","form","h1","h2","h3","h4","h5","h6","header","hr",
"li","main","nav","noscript","ol","p","pre","section","table","tfoot","ul","tr","td","th"
}

INLINE_STYLING_TAGS = ["em", "strong", "b", "i", "u", "mark", "small", "del",
                   "ins", "sub", "sup", "span", "code", "kbd", "samp",
                   "var", "abbr", "cite", "q", "italic"]

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

# ---------- Step 1: search Europe PMC API ----------
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

def search_europe_pmc(query, result_type='core', page_size=25) -> {}:
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
        "format": 'json',
        "pageSize": page_size
    }

    # Add a warning if the user explicitly requests 'lite' with 'json'
    if result_type == 'lite':
        print("Warning: result_type='lite' does not include abstracts. "
              "Filtering will likely result in an empty list. "
              "Consider using result_type='core' or 'full'.")

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raise an exception for HTTP errors
        data = response.json()

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
        return data

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

class EuropePMCFetchError(Exception):
    pass


def _build_session(
    total_retries: int = 5,
    backoff_factor: float = 0.5,
    timeout_seconds: int = 20,
) -> requests.Session:
    """
    Create a requests Session with retries and sensible defaults.
    """
    retry = Retry(
        total=total_retries,
        read=total_retries,
        connect=total_retries,
        backoff_factor=backoff_factor,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET", "HEAD"]),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session = requests.Session()
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update(
        {
            "User-Agent": "Outlier-ModelPlayground/1.0 (+https://www.outlier.ai) PythonRequests",
            "Accept": "application/xml,text/xml;q=0.9,application/xhtml+xml;q=0.8,*/*;q=0.7",
        }
    )
    # Attach a default timeout to session via a wrapper
    session.request = _with_timeout(session.request, timeout_seconds)  # type: ignore
    return session


def _with_timeout(request_func, timeout: int):
    """
    Wrap requests.Session.request to always apply a timeout unless one is explicitly provided.
    """
    def wrapper(method, url, **kwargs):
        if "timeout" not in kwargs:
            kwargs["timeout"] = timeout
        return request_func(method, url, **kwargs)
    return wrapper


def _normalize_pmcid(pmcid: str) -> str:
    """
    Normalize a PMCID string. Accepts values like 'PMC12345' or '12345' and returns 'PMC12345'.
    """
    pmcid = pmcid.strip().upper()
    if pmcid.startswith("PMC"):
        return pmcid
    if pmcid.isdigit():
        return f"PMC{pmcid}"
    raise ValueError(f"Invalid PMCID format: {pmcid!r}")


def strip_tags(xml_string: str) -> str:
    return re.sub(r"<[^>]*>", "", xml_string)


def fetch_full_text_xml(pmcid: str, session: requests.Session | None = None) -> str:
    """
    Fetch the full text (JATS XML) from Europe PMC and return extracted plain text.
    Prefers abstract + body content. Falls back to entire XML text if structure is missing.
    """
    pmcid = _normalize_pmcid(pmcid)
    close_session = False
    if session is None:
        session = _build_session()
        close_session = True

    try:
        # Primary: JATS XML
        xml_url = f"{EUROPE_PMC_BASE}/{pmcid}/fullTextXML"
        resp = session.get(xml_url)
        if resp.status_code == 404:
            raise EuropePMCFetchError(f"PMCID not found at Europe PMC: {pmcid}")
        if not resp.ok or not resp.text.strip():
            # Fallback: HTML, in case XML endpoint is unavailable for this article
            html_url = f"{EUROPE_PMC_BASE}/{pmcid}/fullTextHTML"
            html_resp = session.get(html_url)
            if not html_resp.ok or not html_resp.text.strip():
                raise EuropePMCFetchError(
                    f"Unable to fetch full text for {pmcid} (XML status {resp.status_code}, "
                    f"HTML status {html_resp.status_code})"
                )
            return strip_tags(html_resp.text)

        return strip_tags(resp.text)
    finally:
        if close_session:
            session.close()



def compile_term_patterns(terms: Iterable[str]) -> List[Tuple[str, re.Pattern]]:
    """
    Compile regex patterns for each term:
    - Single-word terms: use word boundaries
    - Multi-word or non-alphanumeric terms: escaped literal substring search (case-insensitive)
    Returns list of (original_term, compiled_pattern).
    """
    patterns: List[Tuple[str, re.Pattern]] = []
    for term in terms:
        t = term.strip()
        if not t:
            continue
        # Use casefold for better Unicode-insensitive matching through re.IGNORECASE
        if re.fullmatch(r"[A-Za-z0-9_]+", t):
            # Single token -> word boundaries
            pat = re.compile(rf"\b{re.escape(t)}\b", flags=re.IGNORECASE)
        else:
            # Phrase or contains punctuation -> plain escaped search
            pat = re.compile(re.escape(t), flags=re.IGNORECASE)
        patterns.append((term, pat))
    return patterns


def find_terms_in_text(text: str, terms: Iterable[str]) -> Dict[str, Dict[str, object]]:
    """
    Check presence of terms in text.
    Returns a dict: term -> {found: bool, count: int, spans: List[(start, end)]}
    """
    results: Dict[str, Dict[str, object]] = {}
    patterns = compile_term_patterns(terms)
    for original_term, pattern in patterns:
        spans: List[Tuple[int, int]] = [m.span() for m in pattern.finditer(text)]
        results[original_term] = {
            "found": bool(spans),
            "count": len(spans),
            "spans": spans,
        }
    return results


def check_terms_in_pmcid(pmcid: str, terms: Iterable[str]) -> Dict[str, Dict[str, object]]:
    """
    High-level function: fetch text for a PMCID and check for terms.
    """
    text = fetch_full_text_xml(pmcid)
    # we dont want a match in references
    return find_terms_in_text(text.replace("REFERENCES", "References").split('References')[0], terms)


def find_candidates(terms, page_size: int = 25, sleep_between_api=1.3 ) -> List[PaperMatch]:

    qry = build_europepmc_query(terms, ENV_KEYWORDS)
    results = search_europe_pmc(qry, page_size=page_size)
    if not results:
        print(f"Warning: Europe PMC query failed for terms {terms}", file=sys.stderr)
        return []

    hits = results.get("resultList", {}).get("result", []) if results else []

    paper_matches: List[PaperMatch] = []
    for rec in hits:
        pm = extract_papermatch_from_record(rec)
        if pm.score <= 0: continue
        paper_matches.append(pm)
    # sort by heuristics score descending
    paper_matches.sort(key=lambda x: x.score, reverse=True)

    # be polite
    sleep(sleep_between_api)

    return paper_matches


def run():
    page_size = 25
    ct = 0
    # TODO create paper entry in the database
    for abr_mutation in AntibioticResMutation.objects.all():
        mutation_terms = [abr_mutation.gene.name, abr_mutation.mutation]
        print(mutation_terms)
        candidate_papers = find_candidates(mutation_terms, page_size)
        print(f"Number of candidate papers:", len(candidate_papers))

        for cp in candidate_papers[:5]:

            sleep(2.7)
            # TODO if the pmid already in the database - continue
            try:
                results = check_terms_in_pmcid(cp.pmcid, mutation_terms + ENV_KEYWORDS + DEPOSITION_MARKERS)
            except (EuropePMCFetchError, ValueError) as e:
                print(f"Error: {e}")
                continue

            for term, info in results.items():
                if info['count'] == 0: continue
                print(f"- {term}: found={info['found']} count={info['count']}")

            terms_found = set(results.keys())
            mutations_found = terms_found.intersection(mutation_terms)
            if len(mutations_found) != len(mutation_terms): continue

            env_terms_found = terms_found.intersection(ENV_KEYWORDS)
            if len(env_terms_found) ==0: continue

            deposition_terms_found = terms_found.intersection(DEPOSITION_MARKERS)
            if len(deposition_terms_found) == 0: continue
            print()
            print(cp.pmid)
            print(cp.title)

            print(mutations_found)
            print(env_terms_found)
            print(deposition_terms_found)
            # TODO store the paper entry and the through table
        print()
        if (ct := ct +1)==5:  exit()


def test_run_0(pmcid = "PMC12284700"):
    # Define your list elsewhere; included here for demonstration.
    TERMS_LIST = ["gyrA", "S81L"]

    try:
        results = check_terms_in_pmcid(pmcid, TERMS_LIST + ENV_KEYWORDS + DEPOSITION_MARKERS)
    except (EuropePMCFetchError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(2)

    for term, info in results.items():
        if info['count'] == 0: continue
        print(f"- {term}: found={info['found']} count={info['count']}")

def test_run():
    mutation_terms = ["gyrA", "S81L"]
    page_size = 100
    candidate_papers = find_candidates(mutation_terms, page_size)
    print(f"Number of candidate papers:", len(candidate_papers))

    for cp in candidate_papers[:5]:
        if not cp.pmcid: continue
        print("--------------------------------------")
        print(cp.pmid, cp.pmcid)
        print(cp.title)
        try:
            results = check_terms_in_pmcid(cp.pmcid, mutation_terms + ENV_KEYWORDS + DEPOSITION_MARKERS)
        except (EuropePMCFetchError, ValueError) as e:
            print(f"Error: {e}")
            continue

        for term, info in results.items():
            if info['count'] == 0: continue
            print(f"- {term}: found={info['found']} count={info['count']}")
        sleep(2.7)


#######################
def main():
    if len(argv) < 2:
        exit("Tell me what to do - test or run?")
    if argv[1] == 'run':
        run()
    elif argv[1] == 'test':
        # for pmcid in ["PMC11984141", "PMC12284700"]:
        #     print(pmcid)
        #     test_run_0(pmcid)
        # print("******************************")
        test_run()
    else:
        print(f"I don't know what is '{argv[1]}'. What should I do - test or run?")

#######################
if __name__ == "__main__":
    main()
