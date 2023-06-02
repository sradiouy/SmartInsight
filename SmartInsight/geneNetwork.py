import streamlit as st
import pickle
import pandas as pd

import streamlit.components.v1 as components
import networkx as nx
from pyvis.network import Network

from updateDb import find_project, create_client, get_create_persist_collection
from geneCleanDocs import *


def network_options():
    #Find list of existing genes in gene_vs_ids.pkl 
    with open("Pickle/gene_network/project_vs_gene.pkl", "rb") as handle:
        project_gene_df = pickle.load(handle)

    with open("Pickle/gene_network/gene_vs_ids.pkl", "rb") as handle:
        gene_id_df = pickle.load(handle)
    

    #st.write("project_gene_df")
    #st.write(project_gene_df)

    #st.write("gene_id_df")
    #st.write(gene_id_df)


    #Filter by current project
    filtered_pg_df = project_gene_df[project_gene_df['project'] == st.session_state.current_project]

    #st.write("filtered_pg_df")
    #st.write(filtered_pg_df)

    #Get genes for current project
    genes_current_project = filtered_pg_df['gene'].unique()
    
    #Get ids for current project
    ids_current_project = filtered_pg_df['pmid'].unique()
    ids = [i for i in ids_current_project if i == i]
    #st.write(ids)
    #Filter by genes for current project
    filtered_gi_df = gene_id_df[genes_current_project]
    #st.write(filtered_gi_df)
    # Filter by ids for current project if they exist in the index
    filtered_gi_df = filtered_gi_df.loc[filtered_gi_df.index.intersection(ids)]

    #st.dataframe(filtered_gi_df)

    gene_lst = filtered_gi_df.columns.tolist()
    #st.write(gene_lst)
    gene_lst.sort()

    if gene_lst:
        #User selects genes to visualize in network
        selection = st.multiselect('Select gene(s) to visualize', gene_lst)

        # Set info message on initial site load
        if len(selection) < 2:
            st.write('Choose at least 2 genes to start')

        # Create network graph when user selects >= 1 item
        else:
            #Keep dataframe only with selection
            filtered_df = filtered_gi_df.loc[:, selection]
            #st.dataframe(filtered_df)
            return filtered_df
    else:
        st.error("No data available for current project. Update your database or create a new project")


def find_gene_pmids(filtered_df):
    gene_pmids = {}
    
    # Iterate over the columns
    for column in filtered_df.columns:
        # Get the gene name
        gene = column
        
        # Get the PubMed IDs that have a value of 1 for the gene
        pubmed_ids = filtered_df.index[filtered_df[column] == 1].tolist()
        
        # Store the gene and its PubMed IDs in the dictionary
        gene_pmids[gene] = pubmed_ids

    # Print the gene-PubMed ID mappings
    return gene_pmids
        
def search_collection(gene_pmids):
    chroma_client = create_client()
    collection_db = get_create_persist_collection(collection_name="database", openai_api_key=st.session_state.key ,chroma_client=chroma_client)
    documents_dict = {}
    
    for gene, pubmed_ids in gene_pmids.items():
        df = collection_db.get(ids=pubmed_ids, include=["documents"])
        documents_combined = ' '.join(df['documents'])
        documents_dict[gene] = documents_combined

    # returns dictionary of clean tokens
    tokens_dict = cleanwords(documents_dict,"","",3)

    try:
        # returns df to create network
        counts = countwords(tokens_dict)
        #st.write(counts)
        if not counts.empty:
            return counts
        else:
            exit()
    except:
        st.error("No connection found for selected genes")
        exit()


def create_gene_networks(df_select):

    G = nx.from_pandas_edgelist(df_select, 'gene', 'term', 'weight')

    # Initiate PyVis network object
    gene_net = Network(
                    height='500px',
                    width='100%',
                    bgcolor='#FFFFFF',
                    font_color='black'
                    )


    # Create two separate node groups for genes and terms
    gene_labels = df_select['gene'].unique()
    term_labels = df_select['term'].unique()

    gene_nodes = set(node for node in G.nodes() if str(node) in gene_labels)
    term_nodes = set(node for node in G.nodes() if str(node) in term_labels)

    # Assign colors to the gene and term node groups
    gene_color = '#5DADE2'
    term_color = '#7ecfa7'

    # Add nodes to the network and assign colors based on the node group
    for gene in gene_nodes:
        gene_net.add_node(gene, color=gene_color, size = 20)

    for term in term_nodes:
        gene_net.add_node(term, color=term_color, size = 14)

    # Set the color for all edges
    edge_color = '#AED6F1'

    # Add edges to the network
    for edge in G.edges(data=True):
        gene, term, weight = edge
        gene_net.add_edge(gene, term, value=weight['weight'], color=edge_color,title=str(weight['weight']))

    # Generate network with specific layout settings
    gene_net.repulsion(
                    node_distance=250,
                    central_gravity=0.33,
                    spring_length=80,
                    spring_strength=0.05,
                    damping=0.95
                    )
        

    # Save and read graph as HTML file (on Streamlit Sharing)
    try:
        path = '/tmp'
        gene_net.save_graph(f'{path}/pyvis_graph.html')
        HtmlFile = open(f'{path}/pyvis_graph.html', 'r', encoding='utf-8')
        
    # Save and read graph as HTML file (locally)
    except:
        path = '/html_files'
        gene_net.save_graph(f'{path}/pyvis_graph.html')
        HtmlFile = open(f'{path}/pyvis_graph.html', 'r', encoding='utf-8')
        
    # Load HTML file in HTML component for display on Streamlit page
    components.html(HtmlFile.read(), height=500)

def count_ids(df):
    genes = df.columns.tolist()
    common_pmids_df = pd.DataFrame(columns=['gene1', 'gene2', 'weight'])

    # Find the common PubMed IDs between every pair of genes:
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            gene1 = genes[i]
            gene2 = genes[j]
            
            common_pmids = df[(df[gene1] == 1) & (df[gene2] == 1)].index
            row = {'gene1': gene1, 'gene2': gene2, 'weight': len(common_pmids), 'pmids': common_pmids.tolist()}
            common_pmids_df = common_pmids_df.append(row, ignore_index=True)
    filtered_common_pmids_df = common_pmids_df[common_pmids_df['weight'] != 0]

    if filtered_common_pmids_df.empty:
        st.error("No connection found for selected genes")
        exit()
    else:
        return filtered_common_pmids_df

def create_ids_networks(common_pmids_df):

    G = nx.from_pandas_edgelist(common_pmids_df, 'gene1', 'gene2', 'weight')

    # Initiate PyVis network object
    gene_net = Network(
                    height='500px',
                    width='100%',
                    bgcolor='#FFFFFF',
                    font_color='black'
                    )


    # Create two separate node groups for genes and terms
    gene1_labels = common_pmids_df['gene1'].unique()
    gene2_labels = common_pmids_df['gene2'].unique()

    gene1_nodes = set(node for node in G.nodes() if str(node) in gene1_labels)
    gene2_nodes = set(node for node in G.nodes() if str(node) in gene2_labels)

    # Assign colors to the gene and term node groups
    gene_color = '#5DADE2'

    # Add nodes to the network and assign colors based on the node group
    for gene in gene1_nodes:
        gene_net.add_node(gene, color=gene_color, size = 20)
    
    for gene2 in gene2_nodes:
        gene_net.add_node(gene2, color=gene_color, size = 20)

    # Set the color for all edges
    edge_color = '#AED6F1'

    # Add edges to the network
    for edge in G.edges(data=True):
        gene, term, weight = edge
        edge_data = common_pmids_df[(common_pmids_df['gene1'] == gene) & (common_pmids_df['gene2'] == term)]
        pmids = edge_data['pmids'].tolist()
        gene_net.add_edge(gene, term, value=weight['weight'], color=edge_color, title=f"IDs: {pmids}")


    # Generate network with specific layout settings
    gene_net.repulsion(
                    node_distance=250,
                    central_gravity=0.33,
                    spring_length=80,
                    spring_strength=0.05,
                    damping=0.95
                    )
        

    # Save and read graph as HTML file (on Streamlit Sharing)
    try:
        path = '/tmp'
        gene_net.save_graph(f'{path}/pyvis_graph.html')
        HtmlFile = open(f'{path}/pyvis_graph.html', 'r', encoding='utf-8')
        
    # Save and read graph as HTML file (locally)
    except:
        path = '/html_files'
        gene_net.save_graph(f'{path}/pyvis_graph.html')
        HtmlFile = open(f'{path}/pyvis_graph.html', 'r', encoding='utf-8')
        
    # Load HTML file in HTML component for display on Streamlit page
    components.html(HtmlFile.read(), height=500)

def gene_network_page(filtered_df):
    gene_pmids = find_gene_pmids(filtered_df)
    counts_df = search_collection(gene_pmids)
    create_gene_networks(counts_df)

def id_network_page(filtered_df):
    #st.dataframe(filtered_df)
    common_pmids_df = count_ids(filtered_df)
    #st.write(common_pmids_df)
    create_ids_networks(common_pmids_df)
