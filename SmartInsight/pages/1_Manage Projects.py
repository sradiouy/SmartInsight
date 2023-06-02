import streamlit as st
from streamlit_extras.stateful_button import button

from chromadb.config import Settings
import chromadb


from updateDb import *
from updateDbPage import *
from utils import validate_key

import sys
sys.path.append('../')

hide_default_format = """
       <style>
       footer {visibility: hidden;}
       </style>
       """
	   #MainMenu {visibility: hidden; }
st.markdown(hide_default_format, unsafe_allow_html=True)


# Initialize session state variables
if 'current_project' not in st.session_state:
       st.session_state.current_project = None
elif st.session_state.current_project is not None and len(st.session_state.current_project) == 0:
	st.session_state.current_project = None

if 'key' not in st.session_state:
    st.session_state.key = None

# Create tooltips for better navigation in the application.
main_help = '''
Update project: add new articles to the current project  \n
Download project: download csv file with information for all articles in the project \n
Delete project: remove existing project \n
'''.strip()

search_option_help = '''
Choose input method for PubMed article retrieval
'''.strip()

homologous_file_help = '''
Upload a Fasta file or a Text file. Text files must have gene IDs separated by new lines
'''.strip()

# Title of the application
st.title("SmartInsight")

selected_page = st.sidebar.radio("Select a page:", ["Update project","Download project","Delete project"],help=main_help)

if selected_page == "Update project":
    st.header("Update Project")

    if not validate_key():
        st.error("Invalid key, set a new key in 'Home'")
    elif st.session_state.current_project is None:
        st.error("Choose a project to update in 'Home'")

    else:
        st.write(f"**Select settings to update project '{st.session_state.current_project}':**")
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
                max = st.slider('Pick number of articles to be retrieved for each gene', min_value = 1, max_value=100, step=1)
                sort = st.selectbox('Pick sort method', ['Relevance', 'Date'])
                uploaded_file = st.file_uploader("Choose a file...",['txt'])
            if st.button(f"Update project '{st.session_state.current_project}'",key="submit"):
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
            
            if st.button(f"Update project '{st.session_state.current_project}'",key="submit"):
                if ids is not None:
                    chroma_client = create_client()
                    articles_df = search_by_id_db(ids,st.session_state.key,chroma_client)
                else:
                    st.error("Enter PubMed IDs")

        elif option == "Input keywords":
            max = st.slider('Pick number of articles to be retrieved for each keyword', min_value=1,max_value=100, step=1)
            sort = st.selectbox('Pick sort method', ['Relevance', 'Date'])
            uploaded_file = st.file_uploader("Choose a file...",['txt'])
            if st.button(f"Update project '{st.session_state.current_project}'",key="submit"):
                chroma_client = create_client()
                articles_df = search_pubmed(uploaded_file,"keywords",max,sort,st.session_state.key,chroma_client)

elif selected_page == "Download project":
    st.header("Download a project")
    project_lst = find_project()
    if not validate_key():
        st.error("Invalid key, set a new key in 'Home'")
    else:
        selection = st.selectbox("**Choose project to download**", project_lst)
        if st.button("Submit"):
            if selection is not None:
                chroma_client = create_client()
                collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
                df = filter_db_project(collection_db,selection)
                df = df.drop(['embeddings','analysis'], axis=1)
                st.dataframe(df)
                csv_file = convert_df(df)
                st.download_button(
                    label=f"Download project '{selection}' as CSV",
                    data=csv_file,
                    file_name=f'project_{selection}.csv',
                    mime='text/csv',
                    )
            else:
                st.error("Choose a project to download")

elif selected_page == "Delete project":
    st.header("Delete a project")
    project_lst = find_project()
    project = st.selectbox("**Choose project to delete**",project_lst)
    if len(project) > 0:
        if st.button(f"Delete project '{project}'"):
            chroma_client = create_client()
            collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
            delete_project(project,collection_db)
        