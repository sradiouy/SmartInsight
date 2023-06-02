import os
import streamlit as st
import pandas as pd

from updateDb import create_client, get_create_persist_collection
 
def find_analysis():
    folder_path = "Pickle/cluster_analysis/"
    names = []
    for filename in os.listdir(folder_path):
        # Extract the name from the filename
        name = os.path.splitext(filename)[0]
        names.append(name)
    return names

def find_cluster(name):
    
    filename = f"Pickle/cluster_analysis/{name}.pkl"

    # Load the pickle file as a DataFrame
    df = pd.read_pickle(filename)
    topics = df['topic']

    topic_dict = {}

    for topic in topics:
        filtered_df = df[df['topic'] == topic]
        id_lst = [id_str for ids in filtered_df['ids'] for id_str in ids.split(';')]
        topic_dict[topic] = id_lst
        #st.write(id_lst)

    #st.write(topic_dict)
    
    return topic_dict

def query_database(id_lst):
    chroma_client = create_client()
    collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
    data = collection_db.get(ids=id_lst, include = ["documents","embeddings","metadatas"])
    #st.write(data)
    return data

def create_temp_database(data):
    chroma_client = create_client()
    temp_collection = get_create_persist_collection(collection_name="temp", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
    
    documents = data['documents']
    embeddings = data['embeddings']
    metadatas = data['metadatas']
    ids = data['ids']

    # Add articles to the Chroma collection
    temp_collection.upsert(
        documents=documents,
        embeddings=embeddings,
        metadatas=metadatas,
        ids=ids
    )


def query_database_project(project):
    chroma_client = create_client()
    collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
    data = collection_db.get()
    df = pd.DataFrame.from_dict(data)
    df = pd.concat([df.drop(['metadatas'], axis=1), df['metadatas'].apply(pd.Series)], axis=1)

    project_filter = df["analysis"].str.contains(project)
    filtered_df = df[project_filter]
    id_lst = filtered_df["ids"].tolist()
    return id_lst

def query_database_question(text):
    n_return = 5
    chroma_client = create_client()
    temp_collection = get_create_persist_collection(collection_name="temp", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
    if temp_collection.count() < n_return:
        result = temp_collection.query(query_texts=[text],n_results=temp_collection.count())
    else:
        result = temp_collection.query(query_texts=[text],n_results=n_return)
    return result

def create_message(id, document, metadata):
    message_content = f"ID: {id}\nDocument: {document}\nTitle: {metadata['title']}\nCitation: {metadata['citation']}\nYear: {metadata['year']}\npmc: {metadata['pmc']}\nmesh: {metadata['mesh']}\nDOI: {metadata['doi']}\nurl: {metadata['url']}"
    #return {"role": "system", "content": message_content}
    return message_content

def create_context_messages(data):
    id_list = data['ids'][0]
    documents = data['documents'][0]
    metadatas = data['metadatas'][0]

    context_messages = []
    for id, document, metadata in zip(id_list, documents, metadatas):
        #st.write("id")
        #st.write(id)
        #st.write("id list")
        #st.write(id_list)
        context_messages.append(create_message(id, document, metadata))
    #st.write(context_messages)
    #for message in context_messages:
    #    st.write(message['content'])
    chroma_client = create_client()
    chroma_client.delete_collection(name="temp")
    return context_messages


    
