import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from kneed import KneeLocator
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from IPython.display import display_markdown
import ast
import pickle
import streamlit as st
import mpld3
import streamlit.components.v1 as components

def cluster_embeddings(df, max_clusters, n_components):
    # Stack the embedding values from the DataFrame into a matrix
    matrix = np.vstack(df.embeddings.values)
    
    # Perform PCA on the matrix with the specified number of components
    pca = PCA(n_components=n_components)
    vis_dims_PCA = pca.fit_transform(matrix)

    # Set KMeans clustering parameters
    kmeans_kwargs = {
           "init": "k-means++",
           "n_init": 10,
           "max_iter": 300,
           "random_state": 42}

    # Perform KMeans clustering for a range of cluster numbers and calculate SSE
    sse = []
    for k in range(1, max_clusters + 1):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(vis_dims_PCA)
        sse.append(kmeans.inertia_)
        
    silhouette_coefficients = []
    # Notice you start at 2 clusters for silhouette coefficient
    for k in range(2, max_clusters):
         kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
         kmeans.fit(matrix)
         score = silhouette_score(matrix, kmeans.labels_)
         silhouette_coefficients.append(score)

    # Return the matrix, SSE, and vis_dims_PCA
    return matrix, sse, vis_dims_PCA, silhouette_coefficients

def plot_sse(max_clusters,sse):
    sns.set(style='white', context='notebook', palette='Greens', font='sans-serif', font_scale=1.2)
    fig, ax = plt.subplots()
    ax.plot(range(1, max_clusters+1), sse)
    ax.set_xticks(range(1, max_clusters+1))
    ax.set_xticklabels(range(1, max_clusters+1), fontsize=6)
    ax.set_xlabel("Number of Clusters")
    ax.set_ylabel("SSE")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.set_facecolor('white')
    #plt.show()
    #plt.savefig('sse.png')
    st.write("Elbow method")
    #st.image('sse.png')
    #st.pyplot(fig)
    fig_html = mpld3.fig_to_html(fig)
    components.html(fig_html, height=600)

def plot_silhouette(max_clusters,silhouette_coefficients):
    sns.set(style='white', context='notebook', palette='Greens', font='sans-serif', font_scale=1.2)
    fig, ax = plt.subplots()
    ax.plot(range(2, max_clusters), silhouette_coefficients)
    ax.set_xticks(range(2, max_clusters))
    ax.set_xticklabels(range(2, max_clusters), fontsize=6)
    ax.set_xlabel("Number of Clusters")
    ax.set_ylabel("Silhouette Coefficient")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.set_facecolor('white')
    #plt.show()
    #plt.savefig('silhouette.png')
    st.write("Optimization of the silhouette coefficient")
    #st.image('silhouette.png')
    fig_html = mpld3.fig_to_html(fig)
    components.html(fig_html, height=600)

def get_nCluster(max_clusters, sse):
    kl = KneeLocator(range(1, max_clusters+1), sse, curve="convex", direction="decreasing")
    optimal = kl.elbow
    return optimal

def find_optimal_clusters(silhouette_coefficients):
    max_score = max(silhouette_coefficients)
    optimal_cluster = silhouette_coefficients.index(max_score) + 2
    return optimal_cluster

def fit_kmeans(df,vis_dims3_PCA, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters, init="k-means++", random_state=42, n_init=10)
    kmeans.fit(vis_dims3_PCA)
    labels = kmeans.labels_
    df["Cluster"] = labels
    return df

def plot_3d_pca_clusters(df, vis_dims3_PCA):
    # Set the seaborn style
    sns.set(style='white', context='notebook', palette='deep', font='sans-serif', font_scale=1.2)
    x = [x for x, y, z in vis_dims3_PCA]
    y = [y for x, y, z in vis_dims3_PCA]
    z = [z for x, y, z in vis_dims3_PCA]

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    unique_clusters = sorted(df['Cluster'].unique())  # Ensure clusters are sorted

    legend_elements = []

    for cluster, color in zip(unique_clusters, sns.color_palette("deep", n_colors=len(unique_clusters))):
        xs = np.array(x)[df.Cluster == cluster]
        ys = np.array(y)[df.Cluster == cluster]
        zs = np.array(z)[df.Cluster == cluster]

        scatter = ax.scatter(xs, ys, zs, color=color, alpha=0.7, edgecolors='w', linewidths=0.5, s=60, label=f"Cluster {cluster}")

        legend_elements.append(scatter)

    ax.set_xlabel('PCA Component 1')
    ax.set_ylabel('PCA Component 2')
    ax.set_zlabel('PCA Component 3')
    plt.title("Clusters identified visualized in language 3D using PCA")
    plt.legend(handles=legend_elements, title='Clusters', loc='upper left', bbox_to_anchor=(1.05, 1), fontsize='small', title_fontsize='medium', markerscale=1.2)
    #plt.show()
    plt.savefig('3dPCA.png')
    st.image('3dPCA.png',width=600)


##### Gaussian Mixture Model ####

