# Publication List with Impact Factor
A python script that fetches all publications of an author from PubMed and adds journal impact factor of the latest Journal Citation Reports (JCR) edition and of the publication year.

## Requirements
- Python >=3.10 with requests, tqdm, python-dotenv and pandas libraries
- Optionally jinja for rendering templates

## Installation
Clone or download the repository and install the required libraries using pip:

```bash
pip install -r requirements.txt
```

### Clarivate API Key
In order to fetch JCR Impact Factors, you will need an API key from Clarivate Analytics.

Create a file named `.env` in the same directory as the script and add your Clarivate API key as follows:

    CLARIVATE_API_KEY=<your_api_key>

Replace `<your_api_key>` with your actual Clarivate API key. You can copy the `sample.env` file in the repository and rename it to `.env`. Replace the placeholder text in the `.env` file with your actual API key.

## Usage

To run the script, use the following command:

```bash
python publication_list.py [-h] \
      -a AUTHORS [AUTHORS ...] \
      -y JCR_YEAR \
      [-o OUTPUT] \
      [--render-template] \
      [--template-filename TEMPLATE_FILENAME] \
      [--template-jif-report {latest,publication_date}] \
      [--summary] \
      [--preprints [PREPRINTS ...]]
```

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
```
  -h, --help            show this help message and exit
  -a AUTHORS [AUTHORS ...], --authors AUTHORS [AUTHORS ...]
                        List of author names
  -y JCR_YEAR, --jcr-year JCR_YEAR
                        JCR year to use for JIF (apart from the respective publication years)
  -o OUTPUT, --output OUTPUT
                        Output file name (default: publication-list.csv)
  --render-template     Render a publication list
  --template-filename TEMPLATE_FILENAME
                        Name of template file to render (default: templates/template.tex.jinja2)
  --template-jif-report {latest,publication_date}
                        Which JIF to report in template file (default: latest)
  --summary             Write JIF summary file
  --preprints [DOI ...]
                        List of medRxiv DOIs to include in the publication list
```

### Examples

The following command will fetch publications of `Lichtner G` and `Lichtner Gregor` using JIF of 2021 and write the output to `custom_output.csv`. It will also write LaTeX file and JIF summary file:

```bash
python publication_list.py \
  -a "Lichtner G" "Lichtner Gregor" \
  -y "2021" \
  -o "custom_output.csv" \
  --render-template \
  --summary
```

The following command will fetch publications of `Lichtner G` using JIF of 2021 and write the output
to `output/publication-list.csv`. It will also render a latex file, output a summary and
include preprints from medRxiv:

```bash
python publication_list.py \
  -a "Lichtner G" \
  -y 2021 \
  -o "output/publication-list.csv" \
  --render-template \
  --summary \
  --preprints "10.1101/2022.07.18.22277750" "10.1101/2022.05.12.22274089"
```




## Limitations
- The JCR Impact Factors are based on the latest edition of the Journal Citation Reports, which is released annually by Clarivate Analytics.
- If a publication year or JCR Impact Factor for a specific journal is not available, the corresponding cells in the CSV file will be left blank.
