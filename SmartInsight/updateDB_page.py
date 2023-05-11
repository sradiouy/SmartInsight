import streamlit as st

from utils import convert_df
from load_articles import *
from FetchArticles import *
from search_update_db import *

from chromadb.config import Settings
from chromadb.utils.embedding_functions import OpenAIEmbeddingFunction
import chromadb


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
    #st.write(data)

    with st.spinner("Computing..."):
        for line in data:
            if how == "genes":
                pmids = search_pmids(line,sort,max)
                all_pmids.append(pmids)
            elif how == "keywords":
                query = query_builder(line)
                #st.write(query)
                pmids = search_pmids(query,sort,max)
                all_pmids.append(pmids)                
        pmid_lst = [pmid for sublist in all_pmids for pmid in sublist]
        unique_pmids = list(set(pmid_lst))
        new_pmids = update_chromaDB(unique_pmids,collection_db)
        
        n_pmids = len(unique_pmids)
        n_new_pmids = len(new_pmids)

        st.write(f"New information retrieved for {n_new_pmids} PMIDs from a total of {n_pmids}")
        
        df =[]
        #all_data = collection_db.get()
        #st.write("all_data")
        #st.write(all_data)
        #st.write("filtered data")
        data = collection_db.get(ids = unique_pmids)
        
        df = pd.DataFrame.from_dict(data)
        df = pd.concat([df.drop(['metadatas'], axis=1), df['metadatas'].apply(pd.Series)], axis=1)
        st.dataframe(df)
        chroma_client.persist()
        
        # Download
        csv_file = convert_df(df)
        st.download_button(
			label="Download new data as CSV",
			data=csv_file,
			file_name=f'{st.session_state.current_analysis}_new_genes.csv',
			mime='text/csv',
		)
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
    st.write(f"New information retrieved for {n_new_pmids} PMIDs from a total of {n_pmids}")
    
    #Retrieve df with new data from database
    data = collection_db.get(ids = unique_pmids)
    df = pd.DataFrame.from_dict(data)
    df = pd.concat([df.drop(['metadatas'], axis=1), df['metadatas'].apply(pd.Series)], axis=1)
    st.dataframe(df)
    chroma_client.persist()

    # Download
    csv_file = convert_df(df)
    st.download_button(
		label="Download new data as CSV",
		data=csv_file,
		file_name=f'{st.session_state.current_analysis}_new_ids.csv',
		mime='text/csv',
	)
    
    return df



def filter_db_analysis(collection,analysis_lst):
    data = collection.get(include=["embeddings","documents","metadatas"])  
    df = pd.DataFrame.from_dict(data)
    df = pd.concat([df.drop(['metadatas'], axis=1), df['metadatas'].apply(pd.Series)], axis=1)
    # Filter rows that contain any of the elements in lst
    filtered_df = df[df['analysis'].str.contains('|'.join(analysis_lst))]
    return filtered_df
