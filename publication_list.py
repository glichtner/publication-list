#!/usr/bin/env python
# coding: utf-8
import os
import re
from typing import Any, Literal

import dotenv
import pandas as pd
import requests
from tqdm import tqdm

dotenv.load_dotenv()

JCR_YEAR = 2021
""" Year for which to fetch the JCR impact factor (should be the latest JCR edition """

JCR_API_URL = "https://api.clarivate.com/apis/wos-journals/v1"
""" API URL for the Journal Citation Reports (JCR) """

PUBMED_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
""" API URL for the PubMed database """

CLARIVATE_API_KEY = os.environ["CLARIVATE_API_KEY"]
"""
API key for the Clarivate API

Either add it here explicitly or add it to a .env file in the same directory as this script.
You can get one here: https://clarivate.com/products/web-of-science-group/
"""

MEDRXIV_API_URL = "https://api.medrxiv.org/details/medrxiv"
""" API URL for the medRxiv preprint server """


def fetch_pubmed_publications(
    author_name: list[str], max_results: int = 100
) -> pd.DataFrame:
    """
    Fetch publication metadata from PubMed for a given list of author names.

    :param author_name: Author name(s) to search for
    :param max_results: Maximum number of results to return from PubMed
    :return: DataFrame with publication metadata
    """
    if not isinstance(author_name, (list, tuple)) or not all(
        isinstance(item, str) for item in author_name
    ):
        raise ValueError("author_name must be a list of strings")

    query = "(" + ") OR (".join([f"{name}[Author]" for name in author_name]) + ")"

    r = requests.get(
        f"{PUBMED_API_URL}?db=pubmed&retmode=json&retmax={max_results}&term={query}"
    )

    res_ids = r.json()["esearchresult"]

    assert int(res_ids["count"]) <= int(
        res_ids["retmax"]
    ), "Found more publications than requested, need to paginate"

    print(f'Found {len(res_ids["idlist"])} publications')

    dfs = []
    keys = [
        "pubdate",
        "epubdate",
        "source",
        "title",
        "fulljournalname",
        "elocationid",
        "issn",
        "essn",
        "lastauthor",
        "sortfirstauthor",
        "volume",
        "issue",
        "pages",
    ]

    for id_ in tqdm(
        res_ids["idlist"], desc="Fetching publication metadata from PubMed"
    ):

        r = requests.get(
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={id_}&retmode=json"
        )

        res_paper = r.json()["result"][id_]

        articleid_doi = [
            aid["value"] for aid in res_paper["articleids"] if aid["idtype"] == "doi"
        ]
        if len(articleid_doi) == 1:
            doi = articleid_doi[0]
        elif len(articleid_doi) == 0:
            doi = ""
        else:
            raise ValueError("Multiple DOIs found")

        authors = [a["name"] for a in res_paper["authors"]]

        df = pd.DataFrame({k: res_paper[k] for k in keys}, index=[0])
        df["authors"] = ", ".join(authors)
        df["doi"] = doi
        df["id"] = id_
        df["pubtype"] = ", ".join(res_paper["pubtype"])

        dfs.append(df)

    df = pd.concat(dfs).reset_index(drop=True)

    author_search = [n.lower() for n in author_name]
    (
        (df["sortfirstauthor"].str.lower().isin(author_search))
        | (df["lastauthor"].str.lower().isin(author_search))
    ).sum()

    df["year"] = df["pubdate"].str.extract(r"(20\d\d)").astype(int)
    df["is_first_or_last_author"] = (
        df["sortfirstauthor"].str.lower().isin(author_search)
    ) | (df["lastauthor"].str.lower().isin(author_search))

    return df


def clarivate_request(
    endpoint: str, params: dict[str, str] | None = None
) -> dict[str, Any]:
    """
    Make a request to the Clarivate API.

    :param endpoint: endpoint to query
    :param params: parameters to pass to the request
    :return: JSON response
    """
    headers = {"X-ApiKey": CLARIVATE_API_KEY}

    r = requests.get(f"{JCR_API_URL}{endpoint}", params=params, headers=headers)

    return r.json()


def get_journal_id(query: str) -> int | None:
    """
    Get the journal ID for a given journal name from JCR.

    :param query: query string to search for
    :return: journal ID or None if no journal was found
    """
    res_jcr = clarivate_request(
        "/journals/",
        params={"edition": "SCIE", "q": query, "count": "100", "firstRecord": "1"},
    )

    if "error" in res_jcr:
        print(res_jcr)
        raise ValueError()

    n_hits = len(res_jcr["hits"])

    if n_hits == 0:
        return None

    assert n_hits == 1, f"Expected 1 hit, found {n_hits}"
    journal_id = res_jcr["hits"][0]["id"]

    return journal_id


def fetch_journal_impact_factor(
    essn: str, issn: str, years: list[str]
) -> list[dict[str, Any]]:
    """
    Fetch the impact factor for a journal for a given year.

    :param essn: Electronic International Standard Serial Number
    :param issn: International Standard Serial Number
    :param years: Years to fetch the impact factor for
    :return: List of dictionaries with the impact factor and the year
    """
    if essn:
        journal_id = get_journal_id(essn)
        if journal_id is None and essn != issn:
            print(f"Warning: couldn't find essn {essn}, trying issn {issn}")
            journal_id = get_journal_id(issn)
    elif issn:
        journal_id = get_journal_id(essn)

    if journal_id is None:
        print(f"Warning: couldn't find essn {essn}/issn {issn}")
        return [{"essn": essn, "jif_year": year} for year in years]

    res = []

    for year in years:
        res_jcr = clarivate_request(f"/journals/{journal_id}/reports/year/{year}")

        jif = float(res_jcr["metrics"]["impactMetrics"]["jif"])
        jname = res_jcr["journal"]["name"]

        res.append(
            {
                "essn": essn,
                "jif_year": res_jcr["year"],
                "journal_name": jname,
                "jif": jif,
            }
        )

    return res


def fetch_publications_with_impact_factor(
    author_name: list[str] | str, jcr_year: int
) -> pd.DataFrame:
    """
    Fetch publications from PubMed and add impact factor from JCR.

    Adds both the impact factor for the given year (jcr_year) and from that
    year of the publication (if available).

    :param author_name: list of author names to search for
    :param jcr_year: year to fetch impact factor for (publication year impact factors are fetched additionally)
    :return: dataframe with publications and impact factors
    """
    if isinstance(author_name, str):
        author_name = [author_name]

    df = fetch_pubmed_publications(author_name)

    # build list of essn, issn and years to fetch from JCR
    essn_years = df.groupby(["essn", "issn"])["year"].unique().explode()
    assert (
        essn_years.reset_index().groupby("essn")["issn"].nunique() == 1
    ).all(), "Found multiple different ISSN per one ESSN"

    dfs = []
    for (essn, issn), g in tqdm(
        essn_years.groupby(level=[0, 1]), desc="Fetching impact factors from JCR"
    ):
        # use only years that are less than the most recent JCR year and add that year explicitly
        years = g[g < jcr_year].values.tolist() + [jcr_year]

        dfs += fetch_journal_impact_factor(essn, issn, years=years)

    df_jif = pd.DataFrame(dfs)

    df_jif_jcr_year = df_jif.loc[
        df_jif["jif_year"] == jcr_year, ["essn", "journal_name", "jif"]
    ]

    # add JIF of {jcr_year} to the main dataframe
    df = pd.merge(
        df,
        df_jif_jcr_year.rename(columns={"jif": f"jif_{jcr_year}"}),
        on="essn",
        how="left",
    )

    # add JIF of the publication year to the main dataframe
    df = pd.merge(
        df,
        df_jif[["essn", "jif_year", "jif"]].rename(
            columns={"jif_year": "year", "jif": "jif_pubyear"}
        ),
        on=["essn", "year"],
        how="left",
    )

    return df.sort_values(by="year", ascending=False)


def fetch_preprints(author_name: list[str], medrxiv_doi: list[str]) -> pd.DataFrame:
    """
    Fetch preprints from medRxiv.

    :param author_name: list of author names
    :param medrxiv_doi: list of medRxiv DOIs
    :return: dataframe with preprints
    """
    res_medrxiv = []

    for doi in tqdm(medrxiv_doi, desc="Fetching preprints from medRxiv"):
        r = requests.get(f"{MEDRXIV_API_URL}/{doi}")

        res_arxiv = r.json()["collection"][0]
        m = re.findall(r"(\w+), ([\w\. ]+)", res_arxiv["authors"])
        authors = ", ".join(
            [ms[0] + " " + ms[1].replace(".", "").replace(" ", "") for ms in m]
        )
        res_medrxiv.append(
            {
                "title": res_arxiv["title"],
                "authors": authors,
                "pubdate": res_arxiv["date"],
                "doi": doi,
                "year": int(res_arxiv["date"][:4]),
                "source": res_arxiv["server"],
                "pubtype": "Preprint",
            }
        )

    df = pd.DataFrame(res_medrxiv)

    authors = df["authors"].str.split(", ").str
    df["is_first_or_last_author"] = authors[0].isin(author_name) | authors[-1].isin(  # type: ignore
        author_name
    )

    return df


def render_template(
    df: pd.DataFrame,
    author_name: list[str],
    template_filename: str,
    report_jif: Literal["latest", "publication_date"] = "latest",
) -> str:
    """
    Render a template with publications.

    :param df: dataframe with publications
    :param author_name: list of author names
    :param template_filename: name of the template file
    :param report_jif: whether to report JIF, and if so, which one
    :return: None
    """
    from jinja2 import Environment, PackageLoader, select_autoescape

    environment = Environment(
        loader=PackageLoader("publication_list"),
        autoescape=select_autoescape(),
    )

    if report_jif == "latest":
        jif_cols = [col for col in df.columns if re.match(r"jif_20\d\d", col)]
        if not jif_cols:
            raise ValueError("No JIF columns found in the dataframe")
        elif len(jif_cols) > 1:
            raise ValueError("Multiple JIF columns found in the dataframe")
        jif_column = jif_cols[0]
        jif_year = jif_column.split("_")[1]
    elif report_jif == "publication_date":
        jif_column = "jif_pubyear"
        jif_year = df["year"]
    else:
        raise ValueError(f"Unknown value for report_jif: {report_jif}")

    environment.tests["contains"] = lambda string, pattern: re.search(pattern, string)
    environment.tests["notempty"] = lambda x: (not pd.isna(x)) & (not x == "")

    template = environment.get_template(template_filename)

    columns = [
        "journal_short",
        "title",
        "authors",
        "doi",
        "pubtype",
        "year",
        "volume",
        "issue",
        "pages",
        "jif",
        "jif_year",
        "is_first_or_last_author",
    ]

    data = df.rename(columns={jif_column: "jif", "source": "journal_short"}).assign(
        jif_year=jif_year
    )

    data = data[columns].sort_values(by="year", ascending=False)

    return template.render(
        data=data.to_dict(orient="records"), author_names=author_name
    )


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Fetch publication data for authors from PubMed with Journal Impact Factor."
    )
    parser.add_argument(
        "-a",
        "--authors",
        nargs="+",
        type=str,
        help="List of author names",
        required=True,
    )
    parser.add_argument(
        "-y",
        "--jcr-year",
        dest="jcr_year",
        type=int,
        help="JCR year to use for JIF (apart from the respective publication years)",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file name (default: publication-list.csv)",
        default="publication-list.csv",
    )
    parser.add_argument(
        "--render-template",
        dest="render_template",
        action="store_true",
        help="Render a publication list",
        default=False,
    )
    parser.add_argument(
        "--template-filename",
        dest="template_filename",
        type=str,
        default="template.tex.jinja2",
        help="Name of template file to render",
    )
    parser.add_argument(
        "--template-jif-report",
        dest="template_jif_report",
        choices=["latest", "publication_date"],
        help="Which JIF to report in template file",
        default="latest",
    )
    parser.add_argument("--summary", action="store_true", help="Write JIF summary file")
    parser.add_argument(
        "--preprints",
        nargs="*",
        type=str,
        help="List of medRxiv DOIs to include in the publication list",
    )

    args = parser.parse_args()

    filename = Path(args.output)

    df = fetch_publications_with_impact_factor(args.authors, args.jcr_year)

    if args.preprints:
        df_preprints = fetch_preprints(
            medrxiv_doi=args.preprints, author_name=args.authors
        )
        df = pd.concat([df, df_preprints], ignore_index=True).reset_index(drop=True)

    df.to_csv(filename, index=False)

    if args.render_template or args.summary:
        idx_exclude = df["pubtype"].str.contains("Published Erratum|Letter")
        print(
            f"Excluding {idx_exclude.sum()} articles for template rendering and summary:"
        )
        print(df.loc[idx_exclude, ["source", "pubtype", "title"]])
        df = df[~idx_exclude]

    if args.render_template:
        rendered = render_template(
            df,
            author_name=args.authors,
            report_jif=args.template_jif_report,
            template_filename=args.template_filename,
        )

        with open(filename.with_suffix(".tex"), "w", encoding="utf-8") as f:
            f.write(rendered)

    if args.summary:
        df_summary = df.groupby("is_first_or_last_author")[
            [f"jif_{args.jcr_year}", "jif_pubyear"]
        ].agg(["sum", "count", "mean", "max"])
        print(df_summary)
        df_summary.to_csv(filename.with_suffix(".summary.csv"), index=True)
