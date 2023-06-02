import streamlit as st
import pickle
import pandas as pd
import numpy as np

def pmids_per_gene(df):
    #st.dataframe(df)
    df = pd.DataFrame(df)

    # Create an empty dictionary to store the relationships
    relationships = {}

    # Iterate over each row in the dataframe
    for _, row in df.iterrows():
        gene = row[0]
        pmids = row[1]
        #If string, split PMIDs
        if isinstance(pmids, str):
            pmids = pmids.split(",")
        #st.write(pmids)
        
        # For each PMID, update the relationships dictionary
        for pmid in pmids:
            if pmid not in relationships:
                relationships[pmid] = {}
            relationships[pmid][gene] = 1

    # Convert the relationships dictionary to a dataframe
    result_df = pd.DataFrame.from_dict(relationships, orient='index').fillna(0)

    return result_df


def gene_vs_ids(df):
    try:
        with open("Pickle/gene_network/gene_vs_ids.pkl", "rb") as handle:
            gene_ids_df = pickle.load(handle)
    except:
            gene_ids_df = pmids_per_gene(df)
    
    #st.write("previous")
    #st.write(gene_ids_df)

    new_df = pmids_per_gene(df)
    #st.write("new")
    #st.write(new_df)

    #Combine previous data with new_df
    merged_df = gene_ids_df.combine(new_df, lambda x1, x2: np.logical_or(x1 == 1, x2 == 1).astype(int))
    #st.write("merged")
    #st.write(merged_df)

    # Update Pickle with merged dataframe
    merged_df.to_pickle("Pickle/gene_network/gene_vs_ids.pkl")



