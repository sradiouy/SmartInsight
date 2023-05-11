import streamlit as st
from streamlit_extras.stateful_button import button
import openai
import pandas as pd
from pandas import DataFrame
from chromadb.config import Settings
import chromadb
import pickle
import sys
import os
sys.path.append('../')


from load_articles_page import *
from generateEmbeddings import update_embeddings
from clustering_page import choose_n_clusters, k_means
from classifyClusters import *
from updateDB_page import *
from search_update_db import find_analysis

import streamlit as st


hide_default_format = """
       <style>
       footer {visibility: hidden;}
       </style>
       """
	   #MainMenu {visibility: hidden; }
st.markdown(hide_default_format, unsafe_allow_html=True)

# Intialize session state
#articles_df = st.session_state.get('articles_df', DataFrame())
embeddings_df = st.session_state.get('embeddings_df', DataFrame())

if 'articles_df' not in st.session_state:
    st.session_state.articles_df = DataFrame()

if 'n_clusters' not in st.session_state:
    st.session_state.n_clusters = None

if 'vis_dims_PCA' not in st.session_state:
    st.session_state.vis_dims_PCA = None

if 'search_done' not in st.session_state:
	st.session_state.search_done = False

if 'embeddings_done' not in st.session_state:
	st.session_state.embeddings_done = False

if 'cluster_df' not in st.session_state:
	st.session_state.cluster_df = DataFrame()

#Create tooltips for better navigation in the application.

main_help = '''
Update analysis: add new articles to the database of the current analysis \n
Retrieve analysis: download files saved in your analysis \n
Find related topics: cluster articles of your analysis in categories \n
Chatbot: interact with your database!
'''.strip()

search_option_help = '''
Choose input method for PubMed article retrieval
'''.strip()

homologous_file_help = '''
Upload a Fasta file or a Text file. Text files must have gene IDs separated by new lines
'''.strip()


# Title of the application
st.title("SmartInsight")

# Sidebar navigation
st.sidebar.title("Navigation")

#Radio button in the sidebar to navigate between three different pages: "Retrieve Data," "Find related topics," and "Chatbot."
selected_page = st.sidebar.radio("Select a page:", ["Update analysis", "Retrieve analysis", "Find related topics","Chatbot"],help=main_help)

if selected_page == "Update analysis":
	st.header("Update Analysis Database")
	if st.session_state.current_analysis is not None:
		st.subheader(f"Add new articles to the analysis '{st.session_state.current_analysis}'")
	else:
		st.error("Select analysis to update in settings")
	option = st.radio("Select an input method:", ["Input genes", "Input PubMed IDs", "Input keywords"],horizontal=True, help=search_option_help)
	if option == "Input genes":
		search_option = st.radio('Use homologous genes for search?', ['Yes', 'No'],horizontal=True)
		if search_option == "Yes":
			left, right = st.columns(2)   
			with left:     
				pcov=st.slider("Minimum percentage of coverage for an homologous to be accepted:",1,100,30)
			with right:
				pident=st.slider("Minimum percentage of identity for an homologous to be accepted:",1,100,30)
			uploaded_file = st.file_uploader("Choose a file...",help=homologous_file_help)
			file_format = st.selectbox('Pick file format', ['fasta', 'text'])
		if search_option == "No":
			max = st.slider('Pick number of articles to be retrieved for each gene', max_value=500, step=5)
			sort = st.selectbox('Pick sort method', ['Relevance', 'Date'])
			uploaded_file = st.file_uploader("Choose a file...",['txt'])
		if st.button("Submit",key="submit"):
			if uploaded_file is not None:
				if search_option == "Yes":
					ids = search_with_homologues(uploaded_file,file_format,pcov,pident)
					chroma_client = create_client()
					articles_df = search_by_id_db(ids,st.session_state.key,chroma_client)
				elif search_option == "No":
					#st.session_state.articles_df = search_by_gene(uploaded_file,max,sort,"Gene")
					chroma_client = create_client()
					articles_df = search_pubmed(uploaded_file,"genes",max,sort,st.session_state.key,chroma_client)	 
			else:
				st.error("Upload a file")
	elif option == "Input PubMed IDs":
		method = st.radio("Select one:", ["Input text", "Upload file",],horizontal=True)
		if method == "Input text":
			ids = retrieve_ids()
		elif method == "Upload file":
			ids = retrieve_ids_file()
		
		if ids is not None:
			chroma_client = create_client()
			articles_df = search_by_id_db(ids,st.session_state.key,chroma_client)

	elif option == "Input keywords":
		max = st.slider('Pick number of articles to be retrieved for each keyword', max_value=500, step=5)
		sort = st.selectbox('Pick sort method', ['Relevance', 'Date'])
		uploaded_file = st.file_uploader("Choose a file...",['txt'])
		if st.button("Submit",key="submit"):
			chroma_client = create_client()
			articles_df = search_pubmed(uploaded_file,"keywords",max,sort,st.session_state.key,chroma_client)	

elif selected_page == "Retrieve analysis":
	st.header("Retrieve analysis")
	analysis_lst = find_analysis()
	selection = st.multiselect('Choose analysis to download', analysis_lst)
	if button("Submit", key="sm"):
		chroma_client = create_client()
		collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
		df = filter_db_analysis(collection_db,selection)
		df = df.drop(['embeddings'], axis=1)
		st.dataframe(df)
		csv_file = convert_df(df)
		name = ",".join(selection)
		st.write(name)
		st.download_button(
			label="Download analysis as CSV",
			data=csv_file,
			file_name=f'analysis_{name}.csv',
			mime='text/csv',
			)

elif selected_page == "Find related topics":
	st.header("Find related topics")
	#if st.session_state.current_analysis is not None:
	#	st.subheader(f"Fin topics in analysis '{st.session_state.current_analysis}'")
	#else:
	#	st.error("Select analysis to update in settings")
	
	analysis_lst = find_analysis()
	selection = st.multiselect('Choose analysis', analysis_lst)
	#st.write(selection)
	if button("Submit", key="sb"):
		chroma_client = create_client()
		collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
		filtered_df = filter_db_analysis(collection_db,selection)
		#st.dataframe(filtered_df)
		n_clusters, vis_dims_PCA = choose_n_clusters(filtered_df,20,3)
		st.session_state['vis_dims_PCA'] = vis_dims_PCA
		chosen_n_clusters = st.slider(label='Pick number of clusters', min_value=1, max_value=20,value=int(n_clusters))
		st.write("Chosen number is ", chosen_n_clusters)
		st.session_state.n_clusters = chosen_n_clusters
		if button("Done", key="Done"):
			#st.write(f"selected number is {st.session_state.n_clusters}")
			st.session_state['cluster_df'] = k_means(filtered_df,st.session_state['vis_dims_PCA'],st.session_state.n_clusters)
				#st.session_state['cluster_df'] = st.session_state['cluster_df'].merge(st.session_state.articles_df,how='left',on='pmid')
				#st.session_state['cluster_df']['abstract'] = articles_df['abstract']
				#st.dataframe(st.session_state['cluster_df'].Cluster)
			df_result = ClassifyClusters(chosen_n_clusters,10,st.session_state['cluster_df'])
			st.dataframe(df_result)
			name = st.text_input('Save clusters as:')
			if st.button("Save clusters"):
				filename = f"Pickle/{name}.pkl"
				if os.path.isfile(filename):
					st.warning("Name already exists, choose a new one")
				else:
					df_result.to_pickle(filename)
					st.success("File saved")


if selected_page == "Chatbot":
	st.error('This is an error', icon="🚨")
	st.warning('This is a warning', icon="⚠️")

