import streamlit as st
from streamlit_extras.stateful_button import button
import openai
import pandas as pd
from pandas import DataFrame


from load_articles_page import *
from generateEmbeddings import update_embeddings
from clustering_page import choose_n_clusters, k_means

import streamlit as st

st.set_page_config(layout="wide", page_title="Application")

hide_default_format = """
       <style>
       footer {visibility: hidden;}
       </style>
       """
	   #MainMenu {visibility: hidden; }
st.markdown(hide_default_format, unsafe_allow_html=True)

# Intialize session state
articles_df = st.session_state.get('articles_df', DataFrame())
embeddings_df = st.session_state.get('embeddings_df', DataFrame())

if 'n_clusters' not in st.session_state:
    st.session_state.n_clusters = None

if 'vis_dims_PCA' not in st.session_state:
    st.session_state.vis_dims_PCA = None

if 'search_done' not in st.session_state:
	st.session_state.search_done = False

if 'embeddings_done' not in st.session_state:
	st.session_state.embeddings_done = False


#Create tooltips for better navigation in the application.

main_help = '''
Retrive data: download article information \n
Find related topics: cluster articles in categories \n
Chatbot: ask all your questions!
'''.strip()

search_option_help = '''
Choose input method for PubMed article retrieval
'''.strip()

homologous_file_help = '''
Upload a Fasta file or a Text file. Text files must have gene IDs separated by new lines
'''.strip()


# Title of the application
st.title("Application")

# Sidebar navigation
st.sidebar.title("Navigation")

#Radio button in the sidebar to navigate between three different pages: "Retrieve Data," "Find related topics," and "Chatbot."
selected_page = st.sidebar.radio("Select a page:", ["Retrieve Data", "Find related topics","Chatbot"],help=main_help)

if selected_page == "Retrieve Data":
	st.header("Pubmed Article Retrieval")
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
					with st.spinner("Computing..."):
						ids = search_with_homologues(uploaded_file,file_format,pcov,pident)
						articles_df = search_by_id(ids)
				elif search_option == "No":
					articles_df = search_by_gene(uploaded_file,max,sort)
			else:
				st.error("Upload a file")
	elif option == "Input PubMed IDs":
		method = st.radio("Select one:", ["Input text", "Upload file",],horizontal=True)
		if method == "Input text":
			ids = retrieve_ids()
		elif method == "Upload file":
			ids = retrieve_ids_file()
		
		if ids is not None:
			articles_df = search_by_id(ids)


	elif option == "Upload file":
		st.subheader("Upload file")
		uploaded_file = st.file_uploader("Choose a file...")
		file_format = st.selectbox('Pick file format', ['fasta', 'text'])
		if st.button("Submit"):
			if uploaded_file is not None:
				with st.spinner("Fetching PMIDs"):
					articles = search_with_homologues(uploaded_file,file_format)
					st.table(articles)
	elif option == "Input keywords":
		articles_df = search_by_gene()


elif selected_page == "Find related topics":
	st.header("Find related topics")
	key = st.sidebar.text_input("Your Open API Key", "sk...")
	openai.api_key = key
	st.subheader("Upload file")
	uploaded_file = st.file_uploader("Choose a file...", type="csv")
	if uploaded_file is not None:
		try:
			df = pd.read_csv(uploaded_file)
			#st.session_state['articles_df'] = df
			if set(['title','citation','year','abstract','pmid','pmc','mesh','doi','url','tag']).issubset(df.columns):
				articles_df = df
				st.session_state['search_done'] = True
			else:
				st.error("Please enter a valid file")
		except:
			st.error("Please enter a valid file")

	if st.session_state['search_done']:
		# Create a button to generate the embeddings:
		st.write("---")
		st.subheader("Generate Embeddings")
		if key == "sk":
			st.error("Please add a valid Open API Key")
		else:
			if button("Generate embeddings", key="Embeddings"):
				with st.spinner("Computing..."):
					embeddings_df = update_embeddings(articles_df)
					#st.table(embeddings_df)
					st.session_state['embeddings_df'] = embeddings_df
					st.session_state['embeddings_done'] = True
					st.success('Embeddings ready')
		
	if st.session_state['embeddings_done']:
		st.write("---")
		st.subheader("Create Clusters")
		if button("Cluster", key="Cluster"):
			n_clusters, vis_dims_PCA = choose_n_clusters(st.session_state['embeddings_df'],20,3)
			st.session_state['vis_dims_PCA'] = vis_dims_PCA
			chosen_n_clusters = st.slider(label='Pick number of clusters', min_value=1, max_value=20,value=int(n_clusters))
			st.write("Chosen number is ", chosen_n_clusters)
			st.session_state.n_clusters = chosen_n_clusters
			if button("Done", key="Done"):
				#st.write(f"selected number is {st.session_state.n_clusters}")
				k_means(st.session_state['embeddings_df'],st.session_state['vis_dims_PCA'],st.session_state.n_clusters)
					
if selected_page == "Chatbot":
	age = st.slider('How old are you?', 0, 130, 25)
	st.write("I'm ", age, 'years old')
	st.error('This is an error', icon="🚨")
	st.warning('This is a warning', icon="⚠️")
	st.info('Info message')
	st.success('Success message')


	
		




