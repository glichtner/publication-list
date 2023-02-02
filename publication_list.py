#!/usr/bin/env python
# coding: utf-8
import os
import re
import textwrap
from typing import Any, Literal

import dotenv
import pandas as pd
import requests
from tqdm import tqdm

dotenv.load_dotenv()

author_name = ["Lichtner G"]
# author_name = ["von Dincklage F", "Von Dincklage Falk", "Dincklage FV", "Dincklage Falk", "Dincklage Falk V"]


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
""" API URL for the medRxiv preprint server"""


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

    return df


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
        authors = (
            ", ".join(
                [ms[0] + " " + ms[1].replace(".", "").replace(" ", "") for ms in m]
            )
            + "."
        )
        res_medrxiv.append(
            {
                "title": res_arxiv["title"],
                "authors": authors,
                "pubdate": res_arxiv["date"],
                "doi": doi,
                "year": res_arxiv["date"][:4],
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


def format_latex(
    df: pd.DataFrame,
    author_name: list[str],
    template_filename: str,
    report_jif: Literal[False, "latest", "publication_date"] = "latest",
) -> str:
    """
    Format the dataframe into a LaTeX file using the given template

    :param df: dataframe with publications
    :param author_name: list of author names
    :param template_filename: name of the template file
    :param report_jif: whether to report JIF, and if so, which one
    :return: None
    """
    from jinja2 import Environment, FileSystemLoader, select_autoescape

    environment = Environment(
        loader=FileSystemLoader("templates/"),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = environment.get_template(template_filename)
    template.render(
        {
            "author_name": author_name,
            "first_last_author_publications": "--",
            "preprints": None,
            "coauthor_publications": "--",
        }
    )

    def format_year_volume_issue_pages(
        year: str, volume: str, issue: str, pages: str
    ) -> str:
        assert (
            (not volume and not issue and not pages)
            or (not issue and volume and pages)
            or (volume and issue and pages)
            or (volume and not issue and not pages)
        ), f"volume: '{volume}', issue: '{issue}', pages: '{pages}'"

        s = f"{year}"

        if volume:
            s += f";{volume}"
        if issue:
            s += f"({issue})"
        if pages:
            s += f":{pages}"

        return f"{s}."

    def row_to_string(row: pd.Series) -> str:
        s = r"\item" + "\n"
        s += row["authors"].replace(author_name, rf"\textbf{{{author_name}}}") + ". "
        s += rf"\textit{{{row['title']}}}" + ", "
        s += row["source"] + ". "
        s += format_year_volume_issue_pages(
            row["year"], row["volume"], row["issue"], row["pages"]
        )

        if not row["doi"] == "":
            s += rf" \studylink{{{row['doi']}}}"

        if not pd.isna(row[jif_column]) and not row[jif_column] == "":
            s += rf" \impact{{{JCR_YEAR if jif_column == 'jif' else row['year']}}}{{{row[jif_column]}}}"

        s += "\n\n"

        return s

    if report_jif == "latest":
        jif_cols = [col for col in df.columns if re.match(r"jif_20\d\d", col)]
        if not jif_cols:
            raise ValueError("No JIF columns found in the dataframe")
        elif len(jif_cols) > 1:
            raise ValueError("Multiple JIF columns found in the dataframe")
        jif_column = jif_cols[0]
    elif report_jif == "publication_date":
        jif_column = "jif_pubyear"
    elif not report_jif:
        jif_column = None
    else:
        raise ValueError(f"Unknown value for report_jif: {report_jif}")

    s = textwrap.dedent(
        r"""
    \subsection{Erst- und Letzt-Autorenschaften}

    \begin{enumerate}[nolistsep]
    \setlength\itemsep{1em}

    """
    )

    for _, row in df.query("is_first_or_last_author").iterrows():
        s += row_to_string(row)

    s += textwrap.dedent(
        r"""
    \end{enumerate}
    \vspace{0.1cm}
    {\footnotesize $\dagger$: gleichwertiger Beitrag}
    \vspace{0.1cm}
    """
    )

    df_preprints = df.query("pubtype == 'Preprint'")
    if not df_preprints.empty:
        s += textwrap.dedent(
            r"""
        \subsection{Preprints}
        \begin{enumerate}[nolistsep]
        \setlength\itemsep{0.9em}
        """
        )

        for _, row in df_preprints.iterrows():
            s += row_to_string(row)

        s += textwrap.dedent(
            r"""
        \newpage

        \end{enumerate}
        """
        )

    s += textwrap.dedent(
        r"""
    \subsection{Co-Autorenschaften}

    \begin{enumerate}[nolistsep]
    \setlength\itemsep{0.9em}
    """
    )

    for _, row in df.query("~is_first_or_last_author").iterrows():
        s += row_to_string(row)

    s += textwrap.dedent(
        r"""
    \end{enumerate}
    """
    )

    return s


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description="Fetch publication data for authors.")
    parser.add_argument(
        "-a", "--authors", nargs="+", type=str, help="List of author names"
    )
    parser.add_argument(
        "-y",
        "--jcr-year",
        dest="jcr_year",
        type=int,
        help="JCR year to use for JIF (apart from the respective publication years)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file name (default: publication-list.csv)",
        default="publication-list.csv",
    )
    parser.add_argument(
        "--latex", action="store_true", help="Write publication list LaTeX file"
    )
    parser.add_argument(
        "--latex-template",
        dest="latex_filename",
        type=str,
        default="latex-template.tex.j2",
        help="Filename of latex template file",
    )
    parser.add_argument(
        "--latex-jif-report",
        dest="latex_jif_report",
        choices=[False, "latest", "publication_date"],
        help="Which JIF to report in latex file",
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

    df = fetch_publications_with_impact_factor(args.author_name, args.jcr_year)

    if args.preprints:
        df_preprints = fetch_preprints(
            medrxiv_doi=args.preprints, author_name=args.author_name
        )
        df = pd.concat([df, df_preprints], ignore_index=True).reset_index(drop=True)

    df.to_csv(filename, index=False)

    if args.latex or args.summary:
        idx_exclude = df["pubtype"].str.contains("Published Erratum|Letter")
        print(f"Exluding {idx_exclude.sum()} articles for latex/summary:")
        print(df.loc[idx_exclude, ["source", "pubtype", "title"]])
        df = df[~idx_exclude]

    if args.latex:
        latex = format_latex(
            df,
            author_name=args.author_name,
            report_jif=args.latex_jif_report,
            template_filename=args.latex_filename,
        )

        with open(filename.with_suffix(".tex"), "w", encoding="utf-8") as f:
            f.write(latex)

    if args.summary:
        df_summary = df.groupby("is_first_or_last_author")[
            [f"jif_{args.jcr_year}", "jif_pubyear"]
        ].agg(["sum", "count", "mean", "max"])
        print(df_summary)
        df_summary.to_csv(filename.with_suffix(".summary.csv"), index=True)
