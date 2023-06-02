import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import seaborn as sns

from updateDb import find_project, create_client, get_create_persist_collection

def gene_pca(filtered_df):
    df = pd.DataFrame(filtered_df)
    transformed_df = pd.DataFrame(columns=['Gene', 'ID'])

    # Iterate over each row in the original DataFrame
    for index, row in df.iterrows():
        for column in df.columns:
            if row[column] == 1:
                transformed_df = transformed_df.append({'Gene': column, 'ID': index}, ignore_index=True)
    #st.dataframe(transformed_df)
    unique_ids = transformed_df['ID'].unique().tolist()

    chroma_client = create_client()
    collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
    data = collection_db.get(ids=unique_ids, include=["embeddings","metadatas"])


    embeddings_df = pd.DataFrame({'ID': data['ids'], 'embeddings': data['embeddings']})
    #st.dataframe(embeddings_df)
    df = transformed_df.merge(embeddings_df, on='ID', how='left')
    #st.dataframe(df)
    df = df.dropna(subset=['embeddings'])
    df = df.reset_index(drop=True)
    #st.dataframe(df)
   
    # Stack the embedding values from the DataFrame into a matrix
    matrix = np.vstack(df.embeddings.values)

    # Perform PCA on the matrix with the specified number of components
    pca = PCA(n_components=3)
    vis_dims_PCA = pca.fit_transform(matrix)
    #st.dataframe(vis_dims_PCA)
    # Get the unique genes in the 'Gene' column
    unique_genes = df['Gene'].unique()

    # Generate a color palette with the number of unique genes
    n_colors = len(unique_genes)
    color_palette = sns.color_palette('hls', n_colors).as_hex()

    # Create a dictionary mapping each gene to a color from the palette
    gene_color_map = dict(zip(unique_genes, color_palette))

    # Create a scatter plot trace for each gene
    scatter_traces = []
    for gene in unique_genes:
        # Filter the DataFrame for the current gene
        gene_df = df[df['Gene'] == gene]
        
        # Create a scatter plot trace for the current gene
        scatter_trace = go.Scatter3d(
            x=vis_dims_PCA[gene_df.index, 0],  # Use the corresponding indices in vis_dims_PCA
            y=vis_dims_PCA[gene_df.index, 1],
            z=vis_dims_PCA[gene_df.index, 2],
            mode='markers',
            hovertext=gene_df['ID'],  # Set the hover text to the 'ID' column
            hoverinfo='text',
            marker=dict(
                size=5,
                color=gene_color_map[gene],  # Use the color from the gene_color_map
                opacity=0.8
            ),
            name=gene  # Set the gene name for the legend
        )
        
        scatter_traces.append(scatter_trace)

    # Create the layout
    layout = go.Layout(
        scene=dict(
            xaxis_title='PC1',
            yaxis_title='PC2',
            zaxis_title='PC3'
        ),
        showlegend=True  # Show the legend
    )

    # Create the figure with the scatter plot traces and layout
    fig = go.Figure(data=scatter_traces, layout=layout)

    # Render the plot in Streamlit
    st.plotly_chart(fig)

    metadata_df = pd.DataFrame.from_dict(data)
    metadata_df  = pd.concat([metadata_df .drop(['metadatas'], axis=1), metadata_df ['metadatas'].apply(pd.Series)], axis=1)
    #st.dataframe(metadata_df)
    
    ids = metadata_df['ids']
    urls = metadata_df['url']

    new_df = pd.DataFrame({'PubMed ID': ids, 'URL': urls})
    new_df = new_df.sort_values('PubMed ID')

    markdown_table = new_df.to_markdown(index=False)

    # Add custom CSS to adjust the font size
    css = """
    <style>
    table {
    font-size: 12px;
    }
    </style>
    """

    # Display the Markdown table with custom CSS
    st.markdown(css, unsafe_allow_html=True)
    st.markdown(markdown_table, unsafe_allow_html=True)