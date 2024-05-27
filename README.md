**PubMed Article Fetcher with MeSH Classification**

This Python script allows users to fetch PubMed articles, extract MeSH (Medical Subject Headings) terms associated with the articles, classify them based on MeSH hierarchy, and store the results in a consolidated CSV file. This README provides an overview of the script and instructions for usage.

### Features
- Fetch PubMed articles based on specified queries and publication years.
- Extract MeSH terms from the fetched articles.
- Classify MeSH terms into three hierarchical levels.
- Store the article data along with MeSH classifications in a CSV file.
- Handle retries and errors gracefully.

### Requirements
- Python 3.x
- Required Python libraries: `metapub`, `tqdm`, `concurrent.futures`, `time`, `Bio`, `requests`, `xml.etree.ElementTree`, `csv`.

### Installation
1. Ensure you have Python installed on your system.
2. Install the required Python libraries using pip:
    ```
    pip install metapub tqdm requests biopython
    ```

### Usage
1. Clone or download the repository to your local machine.
2. Ensure you have the necessary permissions to access PubMed resources.
3. Run the script `PubMed_ter.py` using Python:
    ```
    python PubMed_ter.py
    ```
4. Follow the prompts to provide your email address (required for accessing PubMed) and to specify the search queries and publication years.
5. Once the script completes execution, you will find the consolidated article data in a CSV file named `sorted_articles.csv`.

### Notes
- The script may take some time to execute, depending on the number of articles fetched and the classification process.
- Ensure that you handle the fetched article data responsibly, respecting copyright and usage policies.
- For large-scale data extraction or modifications, consider parallelizing or optimizing the script further.
