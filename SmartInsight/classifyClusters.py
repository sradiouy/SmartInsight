# Example OpenAI Python library request
from IPython.display import display_markdown
import streamlit as st
import openai

openai.api_key = st.session_state.key

@st.cache_data(show_spinner=False)
def ClassifyTopic(abstracts):
    response = openai.ChatCompletion.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "system", "content": "You are a helpful assistant."},
            {"role": "user", "content": "I need your help to summarize a text containing a range of abstracts, define the main topic of them and define keywords. The keywords must be 3, topic must be defined in a maximum of 8 words and summary no more than 50 words."},
            {"role": "assistant", "content": "Hello! I'm here to help you with your scientific research. To get started, please provide me with the list of abstracts you're working with. Once I have that, I'll be able to provide you with the assistance you need. Thank you!:\n"},
            {"role": "user", "content": "Given the following text provide your answer: 'Horizontal genomics is a new field in prokaryotic biology that is focused on the analysis of DNA sequences in prokaryotic chromosomes that seem to have originated from other prokaryotes or eukaryotes. However, it is equally important to understand the agents that effect DNA movement: plasmids, bacteriophages and transposons."},
            {"role": "assistant", "content": """Topic: Horizontal genomics in prokaryotic biology.\n
                    Summary: The abstracts describe the new field of horizontal genomics in prokaryotic biology, which focuses on analyzing DNA sequences that originated from other prokaryotes or eukaryotes. It highlights the importance of understanding the agents that effect DNA movement, such as plasmids, bacteriophages, and transposons, and recent advancements in mobile genetic element (MGE) genomics.\n
                    Keywords: Horizontal genomics; DNA movement; Mobile genetic element (MGE) genomics."""},
            {"role": "user", "content": "Great! Now please give me the main topic, summary and keywords of the following abstracts: \n" + abstracts +". Please provide your response."}
        ],
             temperature=0)
    return response["choices"][0]["message"]["content"]

def ClassifyClusters(n_clusters,abstract_per_cluster,df_clusters):
    for i in range(n_clusters):
        if st.session_state.key is not None:
            try:
                texts = " ".join(
                    df_clusters[df_clusters.Cluster == i].abstract.sample(abstract_per_cluster, random_state=42).values
                )
                #st.write(texts)
                st.write(f"**Cluster {i}:**")
                st.write(ClassifyTopic(texts))
            except Exception as e:
                st.write("____________")
                st.warning(f"Cluster {i}, has less than {abstract_per_cluster} abstracts")
     