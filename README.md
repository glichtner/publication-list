# Publication List with Impact Factor
A python script that fetches all publications of an author from PubMed and adds journal impact factor of the latest Journal Citation Reports (JCR) edition and of the publication year.

## Requirements
- Python 3.x with requests, tqdm, python-dotenv and pandas libraries

## Installation
Clone or download the repository and install the required libraries using pip:

    pip install -r requirements.txt

## Usage

To run the script, use the following command:

    python publication_list.py --authors "<author_name>" [<author_name_2>, ...] --jcr_year 2021

Replace `<author_name>` with the name of the author you want to fetch publications for.

The script will output a CSV file named `publications.csv` in the current directory.
The file will contain amongst others the following columns:
- DOI
- Authors
- Journal Name
- Date of Publication
- Issue, Volume, Pages
- JCR Impact Factor (Latest Edition)
- JCR Impact Factor (Publication Year)

### Arguments

- `-a` / `--authors`: List of author names.
- `-y` / `--jcr_year`: JCR year to use for JIF (apart from the respective publication years).
- `-o` / `--output`: Output file name (default: `publication-list.csv`).
- `--latex-template`: Write publication list to LaTeX template.
- `--latex-jif-report`: Which JIF to report in latex file (default: `latest`).
    - `False`: Do not report JIF in latex file.
    - `latest`: Use JIF of the latest JCR edition.
    - `publication_date`: Use JIF of the year of publication.
- `--summary`: Write JIF summary file (descriptive statistics of author's JIFs).
- `--preprints`: List of medRxiv DOIs to include in the publication list.


### Examples

The following command will fetch publications of `Lichtner G` and `Lichtner Gregor` using JIF of 2021 and write the output to `custom_output.csv`. It will also write LaTeX file and JIF summary file:

    python publication_list.py -a "Lichtner G" "Lichtner Gregor" -y 2021 -o "custom_output.csv" --latex latex-template.tex.jinja --summary


The following command will fetch publications of `Lichtner G` and `Lichtner Gregor` using JIF of 2021 and write the output to `custom_output.csv`. It will also include preprints from medRxiv:

    python publication_list.py -a "Lichtner G" "Lichtner Gregor" -y 2021 -o "custom_output.csv" --preprints "10.1101/2021.01.01.21249673" "10.1101/2021.01.01.21249673"


## Configuration
In order to fetch JCR Impact Factors, you will need an API key from Clarivate Analytics.

Create a file named `.env` in the same directory as the script and add your Clarivate API key as follows:

    CLARIVATE_API_KEY=<your_api_key>

Replace `<your_api_key>` with your actual Clarivate API key.

You can copy the `sample.env` file in the repository and rename it to `.env`. Replace the placeholder text in the `.env` file with your actual API key.

The script will read the API key from the `.env` file when it runs.


## Limitations
- The JCR Impact Factors are based on the latest edition of the Journal Citation Reports, which is released annually by Clarivate Analytics.
- If a publication year or JCR Impact Factor for a specific journal is not available, the corresponding cells in the CSV file will be left blank.
