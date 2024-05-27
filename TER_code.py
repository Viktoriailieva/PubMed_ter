import metapub
import csv
from tqdm import tqdm
import concurrent.futures
import time
from Bio import Entrez
import requests
import xml.etree.ElementTree as ET

start_time = time.time()
class MeshTerm:
    def __init__(self, name, definition, synonyms, tree_numbers):
        self.name = name
        self.definition = definition
        self.synonyms = synonyms
        self.tree_numbers = tree_numbers

def parse_mesh_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    mesh_terms = {}
    tree_number_to_name = {}

    for descriptor in root.findall('DescriptorRecord'):
        name = descriptor.findtext('DescriptorName/String', default='')

        definition = ''
        concept_list = descriptor.find('ConceptList')
        if concept_list is not None:
            definition = concept_list.findtext('Concept/ScopeNote', default='')

        synonyms = []
        if concept_list is not None:
            for term in concept_list.findall('Concept/Term'):
                if term.attrib.get('Type') == 'Synonym':
                    synonyms.append(term.text)

        tree_numbers = []
        tree_number_list = descriptor.find('TreeNumberList')
        if tree_number_list is not None:
            for tn in tree_number_list.findall('TreeNumber'):
                tree_number = tn.text
                tree_numbers.append(tree_number)
                tree_number_to_name[tree_number] = name

        mesh_terms[name] = MeshTerm(name, definition, synonyms, tree_numbers)

    return mesh_terms, tree_number_to_name

# Parse the downloaded XML file
mesh_terms, tree_number_to_name = parse_mesh_xml('desc2024.xml')

def update_classification(level_classes, level, term_name):
    if level not in level_classes:
        level_classes[level] = []
    level_classes[level].append(term_name)

def classify_mesh_terms(mesh_terms, tree_number_to_name):
    major_classes = {
        'A': 'Anatomy',
        'B': 'Organisms',
        'C': 'Diseases',
        'D': 'Chemicals and Drugs',
        'E': 'Analytical, Diagnostic and Therapeutic Techniques, and Equipment',
        'F': 'Psychiatry and Psychology',
        'G': 'Phenomena and Processes',
        'H': 'Disciplines and Occupations',
        'I': 'Anthropology, Education, Sociology, and Social Phenomena',
        'J': 'Technology, Industry, and Agriculture',
        'K': 'Humanities',
        'L': 'Information Science',
        'M': 'Named Groups',
        'N': 'Health Care',
        'V': 'Publication Characteristics',
        'Z': 'Geographicals'
    }

    level_1_classes = {}
    level_2_classes = {}
    level_3_classes = {}

    for term_name, term in mesh_terms.items():
        for tree_number in term.tree_numbers:
            segments = tree_number.split('.')
            for i in range(len(segments)):
                sub_tree_number = '.'.join(segments[:i+1])
                level_name = major_classes.get(sub_tree_number[0], tree_number_to_name.get(sub_tree_number, sub_tree_number))

                if i == 0:
                    update_classification(level_1_classes, level_name, term_name)
                elif i == 1:
                    second_level = tree_number_to_name.get('.'.join(segments[:i]), level_name)
                    update_classification(level_2_classes, second_level, term_name)
                elif i == 2:
                    third_level = tree_number_to_name.get('.'.join(segments[:i]), level_name)
                    update_classification(level_3_classes, third_level, term_name)

    return level_1_classes, level_2_classes, level_3_classes

def fetch_mesh_terms_and_classify(pmid):
    Entrez.email = "your@email.com"  # Provide your email here
    handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    
    try:
        mesh_terms_list = records['PubmedArticle'][0]['MedlineCitation']['MeshHeadingList']
        mesh_annotations = [mesh_term['DescriptorName'] for mesh_term in mesh_terms_list]
        classified_terms = classify_mesh_terms(
            {name: mesh_terms[name] for name in mesh_annotations if name in mesh_terms}, tree_number_to_name)
        return mesh_annotations, classified_terms
    except KeyError:
        return [], ({}, {}, {})

# Combine fetching article data and MeSH terms
def fetch_article_data_with_mesh(pmid, retries=5):
    for attempt in range(retries):
        try:
            fetch = metapub.PubMedFetcher()
            article = fetch.article_by_pmid(pmid)
            mesh_annotations, classified_terms = fetch_mesh_terms_and_classify(pmid)
            return [pmid, article.title, article.abstract, ', '.join(article.authors), ', '.join(mesh_annotations)], classified_terms
        except Exception as e:
            print(f"Error fetching article {pmid}: {e}. Retrying ({attempt + 1}/{retries})...")
            time.sleep(1)
    print(f"Failed to fetch article {pmid} after {retries} retries.")
    return None, ({}, {}, {})

# List of queries
queries = ['Lung cancer', 'Lung Neoplasms', 'Carcinoma, Non-Small-Cell Lung', 'Carcinoma, Small Cell']
years = range(2015, 2016)

# Function to fetch PMIDs for a given query and year
def fetch_pmids(query, year, retries=5):
    for attempt in range(retries):
        try:
            pmids = metapub.PubMedFetcher().pmids_for_query(query, since=year, until=year, retmax=int(1e6))
            return pmids
        except Exception as e:
            print(f"Error fetching PMIDs for query '{query}' in year {year}: {e}. Retrying ({attempt + 1}/{retries})...")
            time.sleep(1)
    print(f"Failed to fetch PMIDs for query '{query}' in year {year} after {retries} retries.")
    return []

# Collect all unique PMIDs and track repetitions
pmid_to_keywords = {}

for query in queries:
    for year in tqdm(years, desc=f'Processing query: {query}'):
        pmids = fetch_pmids(query, year)
        for pmid in pmids:
            if pmid in pmid_to_keywords:
                pmid_to_keywords[pmid].add(query)
            else:
                pmid_to_keywords[pmid] = {query}

unique_pmids = set(pmid_to_keywords.keys())
repeated_pmids = {pmid: keywords for pmid, keywords in pmid_to_keywords.items() if len(keywords) > 1}

print(f"Total unique PMIDs: {len(unique_pmids)}")
print(f"Total repeated PMIDs: {len(repeated_pmids)}")

# Fetch and write unique articles to a consolidated CSV
consolidated_csv_file = 'sorted_articles3.csv'

# Example usage:
with open(consolidated_csv_file, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Keywords', 'PMID', 'Title', 'Abstract', 'Authors', 'Mesh Terms', 'Level 1', 'Level 2', 'Level 3'])
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(fetch_article_data_with_mesh, pmid): pmid for pmid in unique_pmids}
        
        failed_pmids = []
        for future in concurrent.futures.as_completed(futures):
            pmid = futures[future]
            try:
                article_data, classified_terms = future.result()
                if article_data:
                    keywords = ', '.join(pmid_to_keywords[pmid])
                    article_data.insert(0, keywords)  # Insert keywords as the first element
                    level_1, level_2, level_3 = classified_terms
                    # Convert level lists to strings
                    level_1_str = ' ; \t '.join([f"  {class_name}: {', '.join(terms)}" for class_name, terms in level_1.items()])
                    level_2_str = ' ; \t'.join([f"  {class_name}: {', '.join(terms)}" for class_name, terms in level_2.items()])
                    level_3_str = ' ; \t'.join([f"  {class_name}: {', '.join(terms)}" for class_name, terms in level_3.items()])
                    # Write the row to the CSV file
                    writer.writerow(article_data + [level_1_str, level_2_str, level_3_str])
                else:
                    failed_pmids.append(pmid)
            except Exception as e:
                print(f"Failed to process PMID {pmid}: {e}")
                failed_pmids.append(pmid)

# After the concurrent processing section
failed_pmids_retry = []

for pmid in failed_pmids:
    try:
        article_data, classified_terms = fetch_article_data_with_mesh(pmid)
        if article_data:
            keywords = ', '.join(pmid_to_keywords[pmid])
            article_data.insert(0, keywords)  # Insert keywords as the first element
            level_1, level_2, level_3 = classified_terms
            # Convert level lists to strings
            level_1_str = ' ; \t '.join([f"  {class_name}: {', '.join(terms)}" for class_name, terms in level_1.items()])
            level_2_str = ' ; \t'.join([f"  {class_name}: {', '.join(terms)}" for class_name, terms in level_2.items()])
            level_3_str = ' ; \t'.join([f"  {class_name}: {', '.join(terms)}" for class_name, terms in level_3.items()])
            # Reopen the CSV file
            with open(consolidated_csv_file, 'a', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                # Write the row to the CSV file
                writer.writerow(article_data + [level_1_str, level_2_str, level_3_str])
        else:
            failed_pmids_retry.append(pmid)
    except Exception as e:
        print(f"Failed to process PMID {pmid}: {e}")
        failed_pmids_retry.append(pmid)

if failed_pmids_retry:
    print(f"Failed to fetch the following PMIDs on retry: {failed_pmids_retry}")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total time taken: {elapsed_time} seconds")
