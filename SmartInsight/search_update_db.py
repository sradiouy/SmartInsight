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

from utils import generate_random_emails

fetch = PubMedFetcher()

# Create a ChromaDB client with the given settings
@st.cache_resource(show_spinner=False)
def create_client():
    settings = Settings(chroma_db_impl="duckdb+parquet", persist_directory=".chromadb/")
    chroma_client = chromadb.Client(settings)
    return chroma_client

def find_analysis():
    chroma_client = create_client()
    collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
    data = collection_db.get()
    df = pd.DataFrame.from_dict(data)
    df = pd.concat([df.drop(['metadatas'], axis=1), df['metadatas'].apply(pd.Series)], axis=1)
    #st.dataframe(df)
    unique_analysis = list(set(df['analysis'].str.split(';').explode()))
    return unique_analysis 
    
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
    
    if st.session_state.current_analysis == None:
        st.error("No analysis has been chosen")

    return {
        "title": title,
        "citation": citation,
        "year": year,
        "pmc": pmc,
        "mesh": mesh,
        "doi": doi,
        "url": url,
        "analysis": st.session_state.current_analysis
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
            analysis = existing["metadatas"][0]["analysis"]
            #st.write(analysis)
            analysis_lst = analysis.split(";")
            if st.session_state.current_analysis not in analysis_lst:
                new_analysis =analysis+";"+st.session_state.current_analysis
                existing["metadatas"][0]["analysis"] = new_analysis
                chromaDB.update(ids=pmid,metadatas=[existing["metadatas"][0]])
            #st.write(new_analysis)
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


    