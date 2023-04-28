import streamlit as st
import pandas as pd
from metapub import PubMedFetcher
from csv import reader
import backoff
import requests
from Bio import Entrez
from unidecode import unidecode
from bs4 import BeautifulSoup
import pickle

from utils import generate_random_emails

fetch = PubMedFetcher()

@backoff.on_exception(backoff.expo, requests.exceptions.RequestException, max_time=20)
def fetchpmid(pmid):
    """ Fetch an article by supplying its pubmed ID

    Arguments:
        pmid {[int or string]} -- pubmed id

    Returns:
        PubMedArticle object
    """

    return fetch.article_by_pmid(pmid)


def update_data(pmid_df):
    """ Updates 'processed_pmids.csv' file with current run

    Arguments:
    pmid_df {[dataframe]} -- Panda dataframe. Two columns [pubmed id, result]
                             Result: 'True' for correct complete processing
                                     'False' for no information obtained (no abstract, no pmc id)
                                     {PMC ID} for no abstract obtained, but has pmc id

    """

    try:
        df = pd.read_csv('processed_pmids.csv', header=None)
    except:
        df = pd.DataFrame()
    
    updated = pd.concat([df, pmid_df],ignore_index=True)
    updated.to_csv("processed_pmids.csv", index=False, header=False)


def load_processed_data():
    """ Loads previously processed pubmed ids from processed_pmids.csv

    Returns:
        ids_list {[list]} -- list of previously processed pubmed ids
    """

    df = pd.read_csv('processed_pmids.csv', header=None)
    ids_list = df[df.columns[0]].values.tolist()
    return ids_list

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
            pass
        else:
            try:
                if 'Abstract' in xml_data["PubmedArticle"][0]['MedlineCitation']['Article'].keys():
                    abstract = xml_data["PubmedArticle"][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
                    abstract = "\n".join([remove_html(unidecode(xml_data["PubmedArticle"][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][x])) for x in range(0,len(xml_data["PubmedArticle"][0]['MedlineCitation']['Article']['Abstract']['AbstractText']))])
                    return abstract
                else:
                    pass
            except:
                print(pmid,"abstract of interest")
    except:
        pass


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
    return pmids,query


def search_pubmed(pmids,query):
    """Retrieves information from PubMed for each PubMed ID returns a dictionary with the PubMed ID as the key and the article information as the value.
    
    Arguments:
        pmids {[list]} -- list of pmids.

    Returns:
        article_dict {[dict]} -- key: PubMed ID, value: article information
    
    """ 
    try:
        data = load_processed_data()
    except:
        data = []
    lst = []
    article_dict = {}
    

    for pmid in pmids:
        if int(pmid) not in data: 

            try:
                article = fetchpmid(pmid)
            except:
                print(f"Failed to recover pmid {pmid}")
                continue
            
            abstract = get_abstract(pmid)
            pmc = article.pmc
            if not abstract or len(abstract.split()) < 10:
                if pmc != None:
                    st.write(f"Abstract not found or too short for PMID {pmid}, maybe try the PMC: {pmc}")
                    tmp = [pmid] + [pmc]
                    lst.append(tmp)
                else:
                    st.write(f"Abstract not found or too short for PMID {pmid} and PMC not found")
                    tmp = [pmid] + ['False']
                    lst.append(tmp)
                continue

            article_dict[pmid] = {
                'title': None,
                'citation': None,
                'year': None,
                'abstract': None,
                'pmid': None,
                'pmc': None,
                'mesh': None,
                'doi': None,
                'url': None,
                'tag': query
            }
            
            try:
                article_dict[pmid]['title'] = article.title
            except:
                article_dict[pmid]['title'] = None
                
            try:
                article_dict[pmid]['year'] = article.year
            except:
                article_dict[pmid]['year'] = None
                
            try:
                article_dict[pmid]['abstract'] = abstract
            except:
                article_dict[pmid]['abstract'] = None
            
            try:
                article_dict[pmid]['citation'] = article.citation
            except:
                article_dict[pmid]['citation'] = None
            try:
                article_dict[pmid]['pmid'] = article.pmid
            except:
                article_dict[pmid]['pmid'] = None
                
            try:
                article_dict[pmid]['pmc'] = article.pmc
            except:
                article_dict[pmid]['pmc'] = None
                
            try:
                mesh = article.mesh
                if len(mesh) == 0:
                    article_dict[pmid]['mesh'] = None
                else: 
                    article_dict[pmid]['mesh'] = ";".join(mesh.keys())
            except:
                article_dict[pmid]['mesh'] = None
                
            try:
                article_dict[pmid]['doi'] = article.doi
            except:
                article_dict[pmid]['doi'] =None
                
            try:
                article_dict[pmid]['url'] = article.url
            except:
                article_dict[pmid]['url'] = None

            tmp = [pmid] + ['True']
            lst.append(tmp)
            
    df = pd.DataFrame(lst)
    update_data(df)
    return article_dict

def update_database(pmids,data_file,query):
    """Updates database with articles obtained by the given query.

    Arguments:
        query {[string]} -- query for pubmed search
        sort_method {[string]} -- sort method in pubmed database (options: 'relevance','date')
        data_file {[file]} -- Pickle file with information from retrieved articles (title,citation,year,abstract,pmid,pmc,mesh,doi,url,tag)
        max_results {[int]} -- maximum number of pubmed ids to obtain in search

    """

    article_dict = search_pubmed(pmids,query)
    new_pmids = []
    new_pmids = [key for key in article_dict]
    new_df = pd.DataFrame.from_dict(article_dict, orient='index')
    try:
        with open('Pickle/article_db.pkl', 'rb') as handle:
            previous_df = pickle.load(handle)
    except:
        previous_df = pd.DataFrame()
    updated = pd.concat([previous_df, new_df],ignore_index=True)
    updated.to_pickle("Pickle/article_db.pkl")
    return new_pmids