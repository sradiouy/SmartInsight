import streamlit as st
import pandas as pd
from metapub import PubMedFetcher
from csv import reader
import backoff
import requests
import random
import string
from Bio import Entrez
from unidecode import unidecode
from bs4 import BeautifulSoup
import pickle

from utils import convert_df
from load_articles import *
from FetchArticles import *

fetch = PubMedFetcher()

"""
This script provides a user interface that allows users to input PubMed IDs or gene names to retrieve associated articles 
from PubMed. The retrieved information is saved in a Pandas DataFrame object so that it can be queried later.

"""

def open_database():
	"""
	This function loads a pickled Pandas DataFrame object from a file called "article_db.pkl". 
	If this file doesn't exist, it creates an empty dataframe with specified columns (title, citation, year, 
	abstract, pmid, pmc, mesh, doi, url, tag, collection and source) and saves it.
	"""
	try: 
		with open('Pickle/article_db.pkl', 'rb') as handle:
			df = pickle.load(handle)
	except:
		# Define the columns
		columns = ["title", "citation", "year", "abstract", "pmid", "pmc", "mesh", "doi", "url", "tag","collection","source"]
		# Create an empty DataFrame with the specified columns
		df = pd.DataFrame(columns=columns)
		df.to_pickle("Pickle/article_db.pkl")

def search_by_gene(file_input,max,sort,source):
	"""
	This function searches for PubMed articles related to a gene or list of genes provided in a text file. 
	It retrieves article information such as title, citation, year, abstract, pmid, pmc, mesh, doi and url 
	using the metapub and BioPython libraries. It saves the retrieved article information in the pickled DataFrame 
	object. It returns a DataFrame object with the information of each PubMed article which appeared on the search.
	"""
	open_database()
	contents = file_input.read().decode('utf-8')
	genes = [gene.strip() for gene in contents.split("\n") if len(gene.strip()) >= 1]	
	all_pmids = []
	all_new_pmids = []
	with st.spinner("Computing..."):
		for gene in genes:
			pmids,query = search_pmids(gene,sort,max)
			new_pmids = update_database(pmids,'Pickle/article_db.pkl',query,st.session_state.current_collection,source)
			all_pmids.append(pmids)
			all_new_pmids.append(new_pmids)
		pmid_lst = [pmid for sublist in all_pmids for pmid in sublist]
		new_pmid_lst = [pmid for sublist in all_new_pmids for pmid in sublist]
		n_pmids = len(pmid_lst)
		n_new_pmids = len(new_pmid_lst)
		st.write(f"New information retrieved for {n_new_pmids} PMIDs from a total of {n_pmids}")
		with open('Pickle/article_db.pkl', 'rb') as handle:
			df = pickle.load(handle)
		df_filtered = df[df["pmid"].astype(str).isin(pmid_lst)]
		st.dataframe(df_filtered)
		#st.table(df_filtered)
		csv_file = convert_df(df_filtered)
		st.download_button(
			label="Download data as CSV",
			data=csv_file,
			file_name='RetrievedArticles.csv',
			mime='text/csv',
		)
	return df_filtered

def retrieve_ids():
	"""
	This function prompts the user to input one or more PubMed IDs separated by a newline character. 
	If the user clicks the submit button, this function parses the input text area and returns a list of PubMed IDs.
	"""
	pmids = st.text_area("Input PubMed Ids")
	if st.button("Submit"):
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

def search_with_homologues(file,format,pcov,pident):
	"""
	This function searches for homologous genes based on a given gene file (txt or fasta), it performs a search with
	these genes and returns a list of retrieved PubMed IDs.. 
	"""
	articles = ids_by_gene(file, format, pcov, pident)
	#st.write(articles)
	id_list = [str(i) for row in articles for i in row[1].split(",")]
	return id_list

def search_by_id(pmids,source):
	"""
	This function retrieves article information from PubMed such as title, citation, year, abstract, pmid, pmc, mesh, doi and url 
	using the metapub and BioPython libraries, based on a list of PubMed IDs provided as an argument. It saves the retrieved article 
	information in the pickled DataFrame object. It returns a DataFrame object with the information of each PubMed article.
	"""
	open_database()
	#pmids = [str(i) for i in pmids]
	with st.spinner("Computing..."):
		new_pmids = update_database(pmids,'Pickle/article_db.pkl','',st.session_state.current_collection,source)
		n_pmids = len(pmids)
		n_new_pmids = len(new_pmids)
		st.write(f"New information retrieved for {n_new_pmids} PMIDs from a total of {n_pmids}")
		with open('Pickle/article_db.pkl', 'rb') as handle:
			df = pickle.load(handle)
		df_filtered = df[df["pmid"].astype(str).isin(pmids)]
		# Add a download button for the entire table
		csv_file = convert_df(df_filtered)
		st.dataframe(df_filtered)
		#st.table(df_filtered)
		st.download_button(
			label="Download data as CSV",
			data=csv_file,
			file_name='RetrievedArticles.csv',
			mime='text/csv',
		)
	return df_filtered

