# Example OpenAI Python library request
from IPython.display import display_markdown
import streamlit as st
import openai
import pandas as pd

if 'key' not in st.session_state:
    st.session_state.key = None

openai.api_key = st.session_state.key

@st.cache_data(show_spinner=False)
def ClassifyTopic(abstracts):
    response = openai.ChatCompletion.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "system", "content": "You are a helpful assistant."},
            {"role": "user", "content": "I need your help to summarize a text containing a range of abstracts, define the main topic of them and define keywords. The keywords must be 3, topic must be defined in a maximum of 6 words and summary no more than 50 words.\
                                        I will provide you with all abstracts inside a cluster. I want one summary for all of them, a summary for the cluster NOT for each abstract. Keep the topic and keywords as consise as possible"},
            {"role": "assistant", "content": "Hello! I'm here to help you with your scientific research. To get started, please provide me with the abstracts in the cluster you're working with. Once I have that, I'll be able to provide you with the assistance you need. Thank you!:\n"},
            {"role": "user", "content": "Given the following text provide your answer: 'Horizontal genomics is a new field in prokaryotic biology that is focused on the analysis of DNA sequences in prokaryotic chromosomes that seem to have originated from other prokaryotes or eukaryotes. However, it is equally important to understand the agents that effect DNA movement: plasmids, bacteriophages and transposons.\
                                            Plasmids are key vectors of horizontal gene transfer and essential genetic engineering tools. They code for genes involved in many aspects of microbial biology, including detoxication, virulence, ecological interactions, and antibiotic resistance. "},
            {"role": "assistant", "content": """
                    Topic: Horizontal genomics.\n
                    Summary: The abstracts describe the new field of horizontal genomics in prokaryotic biology, which focuses on analyzing DNA sequences that originated from other prokaryotes or eukaryotes. It highlights the importance of understanding the agents that effect DNA movement, such as plasmids, bacteriophages, and transposons, and recent advancements in mobile genetic element (MGE) genomics.\n
                    Keywords: Horizontal genomics; DNA movement; Mobile genetic element (MGE) genomics."""},
            {"role": "user", "content": "Great! Now please give me the main topic, summary and keywords of the following abstracts: \n" + abstracts +". Please provide your response."}
        ],
             temperature=0)
    return response["choices"][0]["message"]["content"]

def ClassifyClusters(n_clusters,abstract_per_cluster,df_clusters):
    #st.dataframe(df_clusters)
    df_result = pd.DataFrame(columns=['cluster', 'topic', 'summary', 'keywords','ids'])
    for i in range(n_clusters):
        if st.session_state.key is not None:
            try:
                cluster_docs = df_clusters[df_clusters.Cluster == i]
                if len(cluster_docs) < abstract_per_cluster:
                    texts = " ".join(df_clusters[df_clusters.Cluster == i].documents.sample(len(cluster_docs), random_state=42).values)
                else:
                    texts = " ".join(df_clusters[df_clusters.Cluster == i].documents.sample(abstract_per_cluster, random_state=42).values)
                #st.write(texts)
                st.write(f"**Cluster {i}:**")
                response = ClassifyTopic(texts)
                st.write(response)
                response_parts = response.split("\n")
                try:
                    topic = response_parts[0].replace("Topic: ", "")
                except:
                    topic = ""
                try:
                    summary = response_parts[2].replace("Summary: ", "")
                except:
                    summary = ""
                try:
                    keywords = response_parts[4].replace("Keywords: ", "")
                except:
                    keywords = ""
                ids = ";".join(cluster_docs.ids.astype(str))
                # Add the information to the DataFrame
                df_result = df_result.append({'cluster': i,'topic': topic,'summary': summary,'keywords': keywords, 'ids':ids}, ignore_index=True)

            except Exception as e:
                st.warning(f"Could not define topics for Cluster {i}")
                #response = ClassifyTopic(texts)
                #st.write(response)
    return df_result
