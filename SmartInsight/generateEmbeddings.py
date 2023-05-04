import re
import pandas as pd
import openai
import tiktoken
from openai.embeddings_utils import get_embedding
import warnings
warnings.filterwarnings('ignore')
import os
import streamlit as st

def expand_contractions(text, contractions):
    for contraction, expansion in contractions.items():
        text = re.sub(contraction, expansion, text, flags=re.IGNORECASE)
    return text

def preprocess_text(text):
    # Lowercase
    text = text.lower()

    # Expand contractions
    contractions = {
        r"\bI'm\b": "I am",
        r"\bI'll\b": "I will",
        r"\bI'd\b": "I would",
        r"\bI've\b": "I have",
        r"\byou're\b": "you are",
        r"\byou'll\b": "you will",
        r"\byou'd\b": "you would",
        r"\byou've\b": "you have",
        r"\bhe's\b": "he is",
        r"\bhe'll\b": "he will",
        r"\bhe'd\b": "he would",
        r"\bhe've\b": "he have",
        r"\bshe's\b": "she is",
        r"\bshe'll\b": "she will",
        r"\bshe'd\b": "she would",
        r"\bshe've\b": "she have",
        r"\bit's\b": "it is",
        r"\bit'll\b": "it will",
        r"\bit'd\b": "it would",
        r"\bwe're\b": "we are",
        r"\bwe'll\b": "we will",
        r"\bwe'd\b": "we would",
        r"\bwe've\b": "we have",
        r"\bthey're\b": "they are",
        r"\bthey'll\b": "they will",
        r"\bthey'd\b": "they would",
        r"\bthey've\b": "they have",
        r"\bthat's\b": "that is",
        r"\bthat'll\b": "that will",
        r"\bthat'd\b": "that would",
        r"\bthere's\b": "there is",
        r"\bthere'll\b": "there will",
        r"\bthere'd\b": "there would",
        r"\bwhere's\b": "where is",
        r"\bwhere'll\b": "where will",
        r"\bwhere'd\b": "where would",
        r"\bwho's\b": "who is",
        r"\bwho'll\b": "who will",
        r"\bwho'd\b": "who would",
        r"\bwhat's\b": "what is",
        r"\bwhat'll\b": "what will",
        r"\bwhat'd\b": "what would",
        r"\bcan't\b": "cannot",
        r"\bwon't\b": "will not",
        r"\bshan't\b": "shall not",
        r"\bdon't\b": "do not",
        r"\bdoesn't\b": "does not",
        r"\bdidn't\b": "did not",
        r"\bhasn't\b": "has not",
        r"\bhaven't\b": "have not",
        r"\bhadn't\b": "had not",
        r"\bwasn't\b": "was not",
        r"\bwere't\b": "were not",
        r"\bwouldn't\b": "would not",
        r"\bshouldn't\b": "should not",
        r"\bcouldn't\b": "could not",
        r"\bmightn't\b": "might not",
        r"\bmustn't\b": "must not",
        r"\bain't\b": "am not",
        r"\baren't\b": "are not",
        r"\bisn't\b": "is not",
        # Add more contractions as needed
    }
    text = expand_contractions(text, contractions)

    # Remove special characters and URLs
    text = re.sub(r"http\S+", "", text)
    text = re.sub(r"[^a-zA-Z0-9\s]+", " ", text)

    # Replace consecutive spaces with a single space and remove leading/trailing spaces
    text = re.sub(r"\s+", " ", text).strip()

    return text

def num_tokens_from_string(string: str, encoding_name: str) -> int:
    """Returns the number of tokens in a text string."""
    encoding = tiktoken.get_encoding(encoding_name)
    num_tokens = len(encoding.encode(string))
    return num_tokens

@st.cache_data(show_spinner=False)
def update_embeddings(df, max_tokens = 8000, encoding = "cl100k_base", model_name= "text-embedding-ada-002"):
    # Check if embeddings DataFrame exists
    pickle_file = 'Pickle/embeddings_db.pickle'
    if os.path.exists(pickle_file):
        embeddings_df = pd.read_pickle(pickle_file)
        pmids_in_pickle = set(embeddings_df['pmid'].values.tolist())
    else:
        embeddings_df = pd.DataFrame(columns=['pmid', 'embedding'])
        pmids_in_pickle = set()
    
    # Iterate over the rows of the input DataFrame
    count = 0
    for index, row in df.iterrows():
        pmid = row['pmid']
        abstract = row['abstract']
        
        # Check if PMID is already in embeddings DataFrame
        if pmid in pmids_in_pickle:
            count += 1
            continue
        
        # Check if abstract is too long to embed
        num_tokens = num_tokens_from_string(abstract, encoding)
        if num_tokens > max_tokens:
            st.write(f"Abstract with PMID {pmid} is too long to embed (has {num_tokens} tokens). Skipping...")
            continue
        
        # Preprocess abstract
        preprocessed_abstract = preprocess_text(abstract)
        
        # Generate embedding
        embedding = get_embedding(preprocessed_abstract, model_name)
        
        # Add pmid and embedding to embeddings DataFrame
        embeddings_df = embeddings_df.append({'pmid': pmid, 'embedding': embedding}, ignore_index=True)
    # remove duplicates rows from embeddings_df
    embeddings_df.drop_duplicates(subset=['pmid'], inplace=True)
    # Save embeddings DataFrame to pickle file
    embeddings_df.to_pickle(pickle_file)
    #st.write(f"Saved {len(embeddings_df)} embeddings to {pickle_file}")
    return embeddings_df
