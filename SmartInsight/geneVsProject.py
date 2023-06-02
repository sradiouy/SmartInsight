import streamlit as st
import pickle
import pandas as pd

def genes_per_project(df):
    df = pd.DataFrame(df)
    
    df_exploded = df.explode(1).reset_index(drop=True)

    new_df = pd.DataFrame({
        'project': [st.session_state.current_project] * len(df_exploded),
        'gene': df_exploded[0],
        'pmid': df_exploded[1]
    })
    #st.write("new_df")
    #st.write(new_df)

    new_df = new_df.drop_duplicates()

    return new_df

def project_vs_gene(df):
    try:
        with open("Pickle/gene_network/project_vs_gene.pkl", "rb") as handle:
            project_gene_id_df = pickle.load(handle)
    except:
            project_gene_id_df = genes_per_project(df)
    
    new_df = genes_per_project(df)

    merged_df = pd.concat([project_gene_id_df, new_df]).drop_duplicates().reset_index(drop=True)
    #st.dataframe(merged_df)

    # Update Pickle with merged dataframe
    merged_df.to_pickle("Pickle/gene_network/project_vs_gene.pkl")

