from chromadb.config import Settings
from chromadb.utils.embedding_functions import OpenAIEmbeddingFunction
import chromadb

import streamlit as st
import pandas as pd
from metapub import PubMedFetcher
import backoff
import requests
from Bio import Entrez
from unidecode import unidecode
from bs4 import BeautifulSoup
import pickle

from utils import generate_random_emails
from updateDbHomologues import *
from geneVsIds import gene_vs_ids
from geneVsProject import project_vs_gene

fetch = PubMedFetcher()

# Create a ChromaDB client with the given settings
@st.cache_resource(show_spinner=False)
def create_client():
    settings = Settings(chroma_db_impl="duckdb+parquet", persist_directory=".chromadb/")
    chroma_client = chromadb.Client(settings)
    return chroma_client

#Find list of existing projects in collection
def find_project():
    chroma_client = create_client()
    collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
    #st.write(collection_db.count())
    data = collection_db.get()
    df = pd.DataFrame.from_dict(data)
    df = pd.concat([df.drop(['metadatas'], axis=1), df['metadatas'].apply(pd.Series)], axis=1)
    #st.dataframe(df)
    unique_project = list(set(df['analysis'].str.split(';').explode()))
    unique_project = [project for project in unique_project if project != ""]
    return unique_project 

def get_create_persist_collection(collection_name, openai_api_key, chroma_client, openai_model_name="text-embedding-ada-002"):
    """
    This function creates a persistent collection in ChromaDB with a given name and an OpenAI embedding function.
    
    Arguments:
    collection_name: A string representing the name of the collection.
    openai_api_key: A string representing the API key for accessing OpenAI.
    openai_model_name: An optional string representing the name of the OpenAI model to use for the embedding function (default is "text-embedding-ada-002").
    """
    # Create an OpenAI EmbeddingFunction object
    openai_ef = OpenAIEmbeddingFunction(api_key=openai_api_key, model_name=openai_model_name)

    # Create a ChromaDB client with the given settings
    #settings = Settings(chroma_db_impl="duckdb+parquet", persist_directory=".chromadb/")
    #chroma_client = chromadb.Client(settings)

    # Get or create a collection with the given name and embedding function
    collection = chroma_client.get_or_create_collection(name=collection_name, embedding_function=openai_ef)

    # Persist the changes to disk
    # chroma_client.persist()

    return collection

def query_builder(keywords):
    # split the keywords into individual words
    words = keywords.split()

    # build the query for each word
    queries = []
    for word in words:
        queries.append(f'"{word}"[All Fields]')

    # combine the queries with OR operators
    query = " AND ".join(queries)

    # add the MeSH term query at the beginning
    query = f'"{keywords}"[MeSH Terms] OR ({query}) OR "{keywords}"[All Fields]'

    return query

def search_with_homologues(file,format,pcov,pident):
    """
    This function searches for homologous genes based on a given gene file (txt or fasta), it performs a search with
    these genes and returns a list of retrieved PubMed IDs.. 
    """
    articles = ids_by_gene(file, format, pcov, pident)
    gene_vs_ids(articles)
    project_vs_gene(articles)
    id_list = [str(i) for row in articles for i in row[1].split(",")]
    return id_list

def retrieve_ids():
	"""
	This function prompts the user to input one or more PubMed IDs separated by a newline character. 
	If the user clicks the submit button, this function parses the input text area and returns a list of PubMed IDs.
	"""
	pmids = st.text_area("Input PubMed Ids")
	id_list = [keyword.strip() for keyword in pmids.split("\n") if len(keyword.strip()) >= 1]
	return id_list

def retrieve_ids_file():
	"""
	This function prompts the user to upload a text file containing one or more PubMed IDs. If the user uploads a file, 
	this function reads the file and returns a list of PubMed IDs.
	"""
	uploaded_file = st.file_uploader("Choose a file...",['txt'])
	if uploaded_file is not None:
		contents = uploaded_file.read().decode('utf-8')
		ids = [id.strip() for id in contents.split("\n") if len(id.strip()) >= 1]
		return ids


def search_pmids(query, sort_method, max_results):
    """Searches for the given query on PubMed and returns a list of PMIDs found in the search. 

    Arguments:
        query {[string]} -- query for pubmed search
        sort_method {[string]} -- sort method in pubmed database (options: 'relevance','date')
        max_results {[int]} -- maximum number of pubmed ids to obtain in search

    Returns:
        pmids {[list]} -- list of pmids. Length = max_results.
    """
    Entrez.email = generate_random_emails()
    handle = Entrez.esearch(db='pubmed',
                            retmax=max_results,
                            retmode='xml',
                            term=query,
                            sort=sort_method)

    results = Entrez.read(handle)
    pmids = results['IdList']
    return pmids

@backoff.on_exception(backoff.expo, requests.exceptions.RequestException, max_time=20)
def fetchpmid(pmid):
    """ Fetch an article by supplying its pubmed ID

    Arguments:
        pmid {[int or string]} -- pubmed id

    Returns:
        PubMedArticle object
    """

    return fetch.article_by_pmid(pmid)

def remove_html(text):
    """Gets the text content of the BeautifulSoup object without HTML code

    Returns:
        {[string]} -- text without HTML code
    """

    soup = BeautifulSoup(text, 'html.parser')
    return soup.get_text()

def get_abstract(pmid):
    """
    Arguments:
        pmid {[type]} -- pubmed id

    Returns:
        abstract {[string]} -- abstract for given pmid without HTML code
    """
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    try:
        xml_data = Entrez.read(handle)

        if xml_data['PubmedArticle'] == []:
            st.warning(f"Failed to recover information for pmid {pmid}")
            pass
        else:
            try:
                if 'Abstract' in xml_data["PubmedArticle"][0]['MedlineCitation']['Article'].keys():
                    abstract = xml_data["PubmedArticle"][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
                    abstract = "\n".join([remove_html(unidecode(xml_data["PubmedArticle"][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][x])) for x in range(0,len(xml_data["PubmedArticle"][0]['MedlineCitation']['Article']['Abstract']['AbstractText']))])
                    return abstract
                else:
                    st.warning(f"Failed to recover information for pmid {pmid}")
                    pass
            except:
                print(pmid,"abstract of interest")
    except:
        st.warning(f"Failed to recover information for pmid {pmid}")
        pass

def get_metadata(pmid):
    try:
        article = fetchpmid(pmid)
    except:
        st.warning(f"Failed to recover information for pmid {pmid}")
        
    try:
        title = article.title
        if title == None:
            title = ""
    except:
        title = ""
                
    try:
        year = article.year
        if year == None:
            year = ""
    except:
        year = ""
                
    try:
        citation = article.citation
        if citation == None:
            citation = ""
    except:
        citation = ""

    try:
        pmc = article.pmc
        if pmc == None:
            pmc = ""
    except:
        pmc = ""

    try:
        lst_mesh = article.mesh
        if len(lst_mesh) == 0:
            mesh = ""
        else: 
            mesh = ";".join(lst_mesh.keys())
    except:
        mesh = ""
    
    try:
        doi = article.doi
        if doi == None:
            doi = ""
    except:
        doi = ""

    try:
        url = article.url
        if url == None:
            url = ""
    except:
        url = ""
    
    if st.session_state.current_project == None:
        st.error("No project has been chosen")

    return {
        "title": title,
        "citation": citation,
        "year": year,
        "pmc": pmc,
        "mesh": mesh,
        "doi": doi,
        "url": url,
        "analysis": st.session_state.current_project
    }

def update_chromaDB(pmids, chromaDB):
    """
    Updates the ChromaDB collection with new PMIDs and associated metadata and abstracts.

    Parameters:
    pmids (list): List of PMIDs to be added to the collection.
    chromaDB: The collection to be updated.
    """
    #chromaDB = get_create_persist_collection("database",key)
    new_pmids = []
    data = chromaDB.get()
    #st.write(data)
    df = pd.DataFrame.from_dict(data)
    #st.write("previous dataframe")
    #st.dataframe(df)
    for pmid in pmids:
        #st.write(pmid)
        # Check if PMID is already in the collection
        if pmid in df['ids'].values:
            #st.write("exists")
            existing = chromaDB.get(ids=pmid)
            #st.write(existing)
            project = existing["metadatas"][0]["analysis"]
            #st.write(project)
            project_lst = project.split(";")
            if st.session_state.current_project not in project_lst:
                new_project =project+";"+st.session_state.current_project
                existing["metadatas"][0]["analysis"] = new_project
                chromaDB.update(ids=pmid,metadatas=[existing["metadatas"][0]])
            #st.write(new_project)
        else:
            new_pmids.append(pmid)
            abstract = get_abstract(pmid)
            metadata = get_metadata(pmid)
            #st.write("abstract")
            #st.write(abstract)
            #st.write("metadata")
            #st.write(metadata)
            if abstract:
                chromaDB.add(
                    documents = abstract,
                    metadatas = metadata,
                    ids = pmid
                    )
    return new_pmids

def delete_project(project,chromaDB):

    #Delete project from database
    data = chromaDB.get()
    #st.write(data)
    df = pd.DataFrame.from_dict(data)
    df = pd.concat([df.drop(['metadatas'], axis=1), df['metadatas'].apply(pd.Series)], axis=1)
    #st.dataframe(df)
    project_filter = df["analysis"].str.contains(project)
    filtered_df = df[project_filter]
    id_lst = filtered_df["ids"].tolist()
    for pmid in id_lst:
        existing = chromaDB.get(ids=pmid)
        #st.write(existing)
        all_projects = existing["metadatas"][0]["analysis"]
        #st.write(all_projects)
        project_lst = all_projects.split(";")
        #st.write(project_lst)
        #st.write("project")
        #st.write(project)
        project_lst.remove(project)
        new_project = ";".join(project_lst)
        existing["metadatas"][0]["analysis"] = new_project
        #st.write(new_project)
        chromaDB.update(ids=pmid,metadatas=[existing["metadatas"][0]])
        
    #Delete project from project_vs_gene Pickle file
    with open("Pickle/gene_network/project_vs_gene.pkl", "rb") as handle:
        df = pickle.load(handle)
    #st.dataframe(df)
    filtered_df = df[df['project'] != project]
    filtered_df.to_pickle("Pickle/gene_network/project_vs_gene.pkl")

    st.success(f"Project '{project}' has been deleted")




    