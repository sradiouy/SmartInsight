import streamlit as st

from utils import convert_df
#rom load_articles import *
from updateDbHomologues import *
from updateDb import *

from chromadb.config import Settings
from chromadb.utils.embedding_functions import OpenAIEmbeddingFunction
import chromadb

from geneVsProject import *


def search_pubmed(file_input,how,max,sort,key,chroma_client):

    # Get or create a collection with the given name and embedding function
    collection_db = get_create_persist_collection(collection_name="database", openai_api_key=key,chroma_client=chroma_client)
    #st.write("previous data from page")
    #prev_data = collection_db.get()
    #st.write(prev_data)
    #st.write("-----------")
    contents = file_input.read().decode('utf-8')
    data = [line.strip() for line in contents.split("\n") if len(line.strip()) >= 1]
    all_pmids = []
    all_new_pmids = []
    articles = []
    #st.write(data)
    with st.spinner("Computing..."):
        for line in data:
            if how == "genes":
                pmids = search_pmids(line,sort,max)
                tmpart = [line] + [pmids]
                articles += [tmpart]
                all_pmids.append(pmids)
            elif how == "keywords":
                query = query_builder(line)
                #st.write(query)
                pmids = search_pmids(query,sort,max)
                all_pmids.append(pmids)  
        
        if how == "genes":
            gene_vs_ids(articles)
            #st.write("articles")
            #st.dataframe(articles)
            project_vs_gene(articles)

        pmid_lst = [pmid for sublist in all_pmids for pmid in sublist]
        unique_pmids = list(set(pmid_lst))
        new_pmids = update_chromaDB(unique_pmids,collection_db)
        
        n_pmids = len(unique_pmids)
        n_new_pmids = len(new_pmids)

        st.success(f"Project {st.session_state.current_project} has been updated")
        st.write(f"New information retrieved for {n_new_pmids} PMIDs from a total of {n_pmids}")
        
        df =[]
        #all_data = collection_db.get()
        #st.write("all_data")
        #st.write(all_data)
        #st.write("filtered data")
        data = collection_db.get(ids = unique_pmids)
        
        df = pd.DataFrame.from_dict(data)
        df = pd.concat([df.drop(['metadatas'], axis=1), df['metadatas'].apply(pd.Series)], axis=1)
        df = df.drop(['embeddings','analysis'], axis=1)
        st.dataframe(df)
        chroma_client.persist()
        
        # Download
        #csv_file = convert_df(df)
        #st.download_button(
		#	label="Download new data as CSV",
		#	data=csv_file,
		#	file_name=f'{st.session_state.current_project}_new_genes.csv',
		#	mime='text/csv',
		#)
        return df
        

def search_by_id_db(pmid_lst,key,chroma_client):
    # Get or create a collection with the given name and embedding function
    unique_pmids = list(set(pmid_lst))
    all_new_pmids = []
    collection_db = get_create_persist_collection(collection_name="database", openai_api_key=key,chroma_client=chroma_client)
    n_pmids = len(unique_pmids)
    with st.spinner("Computing..."):
        new_pmids = update_chromaDB(unique_pmids,collection_db)
    
    n_new_pmids = len(new_pmids)
    st.success(f"Project {st.session_state.current_project} has been updated")
    st.write(f"New information retrieved for {n_new_pmids} PMIDs from a total of {n_pmids}")
    
    #Retrieve df with new data from database
    data = collection_db.get(ids = unique_pmids)
    df = pd.DataFrame.from_dict(data)
    df = pd.concat([df.drop(['metadatas'], axis=1), df['metadatas'].apply(pd.Series)], axis=1)
    df = df.drop(['embeddings','analysis'], axis=1)
    st.dataframe(df)
    chroma_client.persist()

    # Download
    #csv_file = convert_df(df)
    #st.download_button(
	#	label="Download new data as CSV",
	#	data=csv_file,
	#	file_name=f'{st.session_state.current_project}_new_ids.csv',
	#	mime='text/csv',
	#)
    return df



def filter_db_project(collection,selection):
    data = collection.get(include=["embeddings","documents","metadatas"])  
    df = pd.DataFrame.from_dict(data)
    df = pd.concat([df.drop(['metadatas'], axis=1), df['metadatas'].apply(pd.Series)], axis=1)
    # Filter rows that contain the selection (current project)
    filtered_df = df[df['analysis'].str.contains(selection)]
    #st.dataframe(filtered_df)
    return filtered_df
