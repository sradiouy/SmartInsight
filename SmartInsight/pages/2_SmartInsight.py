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

from clusterPage import choose_n_clusters, k_means
from clusterClassification import *
from updateDbPage import *
from updateDb import find_project
from geneNetwork import *
from utils import validate_key
from chatBotPage import chat_page
from genePCA import gene_pca

import streamlit as st


hide_default_format = """
       <style>
       footer {visibility: hidden;}
       </style>
       """
	   #MainMenu {visibility: hidden; }
st.markdown(hide_default_format, unsafe_allow_html=True)

# Intialize session state variables

if 'articles_df' not in st.session_state:
    st.session_state.articles_df = DataFrame()

if 'n_clusters' not in st.session_state:
    st.session_state.n_clusters = None

if 'vis_dims_PCA' not in st.session_state:
    st.session_state.vis_dims_PCA = None

if 'cluster_df' not in st.session_state:
	st.session_state.cluster_df = DataFrame()

if 'project' not in st.session_state:
       st.session_state.project = None

if 'current_project' not in st.session_state:
       st.session_state.current_project = None
elif st.session_state.current_project is not None and len(st.session_state.current_project) == 0:
	st.session_state.current_project = None

if 'cluster_button' not in st.session_state:
	st.session_state.cluster_button = False

if 'key' not in st.session_state:
    st.session_state.key = None

if "generated" not in st.session_state:
    st.session_state.generated = []
if "past" not in st.session_state:
    st.session_state.past = []
if "messages" not in st.session_state:
    st.session_state.messages = []
if "user_text" not in st.session_state:
    st.session_state.user_text = ""

#Create tooltips for better navigation in the application.

main_help = '''
Create gene network: create gene networks by frequent terms and PubMed Ids \n
Find related topics: cluster articles of your project in categories \n
Chatbot: interact with your project or cluster
'''.strip()

analysis_help = '''
Saving the analysis will let you interact with a specific cluster in the chatbot
'''.strip()


# Title of the application
st.title("SmartInsight")

# Sidebar navigation
st.sidebar.title("Navigation")

#Radio button in the sidebar to navigate between pages.
selected_page = st.sidebar.radio("Select a page:", ["Analysis by genes","Find related topics","Chatbot"],help=main_help)

if selected_page == "Find related topics":
	st.header("Find related topics")
	#if st.session_state.current_project is not None:
	#	st.subheader(f"Fin topics in project '{st.session_state.current_project}'")
	#else:
	#	st.error("Select project to update in settings")
	
	if not validate_key():
		st.error("Invalid key, set a new key in 'Home'")
	elif st.session_state.current_project is None:
		st.error("Choose a project to update in 'Home'")
	else:
		st.write(f"**Using project '{st.session_state.current_project}':**")
	
		selection = st.session_state.current_project
		chroma_client = create_client()
		collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
		filtered_df = filter_db_project(collection_db,selection)
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
			df_result = ClassifyClusters(chosen_n_clusters,8,st.session_state['cluster_df'])
			st.dataframe(df_result)
			name = st.text_input('Save analysis as:', help=analysis_help)
			if st.button("Save analysis"):
				filename = f"Pickle/cluster_analysis/{name}.pkl"
				if os.path.isfile(filename):
					st.warning("Name already exists, choose a new one")
				elif len(name)==0:
					st.error("Enter a name to save the analysis")
				else:
					df_result.to_pickle(filename)
					st.success("Analysis saved")


elif selected_page == "Analysis by genes":

	st.header("Create gene network")
	if st.session_state.current_project is not None:
		st.write(f"**Using project '{st.session_state.current_project}':**")
		filtered_df = pd.DataFrame()
		filtered_df = network_options()
		#st.write(filtered_df)
		if filtered_df is not None:
			tab1, tab2, tab3 = st.tabs(["Network by frequent terms", "Network by IDs","Abstract similarity"])
			with tab1:
				with st.expander("Network help"):
					st.write(
					"""
					This interactive network provides a visual representation of gene interactions and their connections to frequent words in scientific articles.
					-	Gene nodes are displayed in blue.
					-	Frequent words nodes are displayed in green.
					-	Thicker edges indicate stronger associations.\n
					Feel free to interact with the network and uncover valuable insights!
					"""
					)
				gene_network_page(filtered_df)
				gene_pca(filtered_df)
			with tab2:
				with st.expander("Network help"):
					st.write(
					"""
					This interactive network provides a visual representation of connections between genes based on shared PubMed IDs.
					-	Gene nodes are displayed in blue.
					-	Edges between gene nodes symbolize shared PubMed IDs.
					-	Stronger edges indicate a higher number of PubMed IDs in common.
					-	Hovering over an edge reveals the associated PubMed IDs.\n
					Feel free to interact with the network and uncover valuable insights!
					"""
					)
				id_network_page(filtered_df)
			with tab3:
				with st.expander("Network help"):
					st.write(
					"""
					This interactive plot provides a three-dimensional view of abstracts using Principal Component Analysis (PCA). 
					-	Each point represents an abstract, color-coded by gene.
					-	Hovering over a point will display the associated PubMed ID.\n
					Feel free to interact with the plot, rotate it, zoom in/out, and explore the abstracts associated with different genes. 
					"""
					)
				gene_pca(filtered_df)
			
	else:
		st.error("Choose a project in 'Home'")

elif selected_page == "Chatbot":
	st.header("ChatBot")
	chat_page()
	#if not validate_key():
	#	st.error("Invalid key, set a new key in 'Home'")
	#else:
	#	chat_page()



