import streamlit as st
from clusterEmbeddings import *

if 'sse' not in st.session_state:
    st.session_state.sse = None

if 'silhouette_coefficients' not in st.session_state:
       st.session_state.silhouette_coefficients = None

@st.cache_data(show_spinner=False)
def choose_n_clusters(df,max_clusters=20,n_components=3):
    with st.spinner("Looking for optimal number of clusters"):
        matrix, st.session_state.sse, vis_dims_PCA, st.session_state.silhouette_coefficients = cluster_embeddings(df, max_clusters, n_components)
        st.write("Optimal number of clusters analysis:")
        #col1, col2 = st.columns(2)
        tab1, tab2 = st.tabs(["Elbow Method", "Optimization of the silhouette coefficient"])
        with tab1:
        #with col1:
            plot_sse(max_clusters, st.session_state.sse)
            n_clusters_elbow= get_nCluster(max_clusters,st.session_state.sse)
            st.write(f"Recommended number of clusters is {n_clusters_elbow}")
        with tab2:
        #with col2:
            plot_silhouette(max_clusters,st.session_state.silhouette_coefficients)
            n_clusters_sil = find_optimal_clusters(st.session_state.silhouette_coefficients)
            st.write(f"Recommended number of clusters is {n_clusters_sil}")

        return int(n_clusters_elbow), vis_dims_PCA


@st.cache_data(show_spinner=False)
def k_means(df,vis_dims_PCA,clusters):
    with st.spinner("Computing clusters"):
        clusters_df = fit_kmeans(df,vis_dims_PCA, clusters)
        plot_3d_pca_clusters(clusters_df, vis_dims_PCA)
    return clusters_df







