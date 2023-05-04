import streamlit as st
from clusterEmbeddings import *


@st.cache_data(show_spinner=False)
def choose_n_clusters(df,max_clusters=20,n_components=3):
    with st.spinner("Looking for optimal number of clusters"):
        matrix, sse, vis_dims_PCA, silhouette_coefficients = cluster_embeddings(df, max_clusters, n_components)
        col1, col2 = st.columns(2)
        with col1:
            plot_sse(sse,max_clusters)
            n_clusters_elbow= get_nCluster(max_clusters,sse)
            st.write(f"Recommended number of clusters is {n_clusters_elbow}")
        with col2:
            plot_silhouette(silhouette_coefficients, max_clusters)
            n_clusters_sil = find_optimal_clusters(silhouette_coefficients)
            st.write(f"Recommended number of clusters is {n_clusters_sil}")

        return int(n_clusters_elbow), vis_dims_PCA


@st.cache_data(show_spinner=False)
def k_means(df,vis_dims_PCA,clusters):
    with st.spinner("Computing clusters"):
        clusters_df = fit_kmeans(df,vis_dims_PCA, clusters)
        plot_3d_pca_clusters(clusters_df, vis_dims_PCA)
    return clusters_df







