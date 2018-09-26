# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 13:30:06 2018

@author: Admin
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from sklearn.neighbors import DistanceMetric
from sklearn.metrics.pairwise import pairwise_distances
from sklearn import cluster
from sklearn.decomposition.pca import PCA
from sklearn.manifold import TSNE
import seaborn as sns
from sklearn import preprocessing

markersizes= [20,200]
markers = ['^', 'o']

def build_images_KMeans(spectra, spectrum_columns, spectra_distances, colors, TSNE_learning_rate=500, TSNE_n_iter=1500, TSNE_learning_rate2=300):
    colors_m = ['red','black']
    
    cols = spectra['marked'].apply(lambda x: colors[x])
    col = spectra['marked']

    plt.subplots(figsize=(18, 6))
    plt.subplot(131)
    plt.title("PCA")
    pca = PCA(n_components=2, random_state=42)
    spectra_2D = pca.fit_transform(spectra[spectrum_columns])
    for i in range(len(spectra_2D)):
        #print(i)
        #print(spectra_2D[i, 0])
        #print(spectra_2D[i, 1])
        #print(cols[i])
        #print(col[i])
        plt.scatter(spectra_2D[i, 0], spectra_2D[i, 1], c=cols[i], alpha=0.5, marker=markers[col[i]], s=markersizes[col[i]])

    plt.subplot(132)
    plt.title("TSNE, Euclidean distance")
    tsne = TSNE(n_components=2, random_state=42, learning_rate=TSNE_learning_rate, n_iter=TSNE_n_iter)
    spectra_2D = tsne.fit_transform(spectra[spectrum_columns])
    for i in range(len(col)):
        plt.scatter(spectra_2D[i, 0], spectra_2D[i, 1], c=cols[i], alpha=0.5, marker=markers[col[i]], s=markersizes[col[i]])

    plt.subplot(133)
    plt.title("TSNE, Chosen distance")
    tsne = TSNE(n_components=2, random_state=42, metric="precomputed", learning_rate=TSNE_learning_rate2)
    spectra_2D = tsne.fit_transform(spectra_distances)
    for i in range(len(col)):
        plt.scatter(spectra_2D[i, 0], spectra_2D[i, 1], c=cols[i], alpha=0.5, marker=markers[col[i]],s=markersizes[col[i]])
    
    # visualization - tsne with chosen distance
    print('Clustering')
    plt.subplots(figsize=(18, 12))

    for n in range(2, 10):
        kmeans = cluster.KMeans(n_clusters=n, random_state=42)
        cluster_labels = kmeans.fit_predict(spectra_distances)

        plt.subplot(3, 3, n-1)
        cols = [colors[l] for l in cluster_labels]
        plt.title("cluster labels ({} clusters)".format(n))
        for i in range(len(col)):
            plt.scatter(spectra_2D[i, 0], spectra_2D[i, 1], c=cols[i], alpha=0.5, marker=markers[col[i]], s=markersizes[col[i]])
    
    plt.show()

    plt.show()
    return spectra_2D

def print_clusters_structure_KMeans(spectra, spectrum_columns, other_names, spectra_distances, n, colors, spectra_2D):
    
    kmeans = cluster.KMeans(n_clusters=n, random_state=42)
    cluster_labels = kmeans.fit_predict(spectra_distances)
    centers = []
    list_spectra_clusters = []
    col = spectra['marked']
    plt.figure(figsize=(18, 6))
    cols = [colors[l] for l in cluster_labels]
    plt.title("cluster labels ({} clusters)".format(n))
    for i in range(len(col)):
        plt.scatter(spectra_2D[i, 0], spectra_2D[i, 1], c=cols[i], alpha=0.5, marker=markers[col[i]], s=markersizes[col[i]])
    
    plt.show()

    spectra['KMeans_clusters_' + str(n)] = cluster_labels
    display(spectra[spectra['marked'] == 1])
    for i in range(n):
        list_spectra_clusters.append(spectra[spectra['KMeans_clusters_' + str(n)] == i][other_names])
        centers.append(spectra.loc[spectra['KMeans_clusters_' + str(n)] == i, spectrum_columns].mean(axis=0).values)
    centers = pd.DataFrame(np.column_stack(centers).T, columns=spectrum_columns)
    del spectra['KMeans_clusters_' + str(n)]

    # centroids
    return list_spectra_clusters, centers
    
def build_images_DBSCAN(spectra, spectrum_columns, spectra_distances, colors, eps_l=[0.01 * i for i in range(12, 0, -2)], TSNE_learning_rate=500, TSNE_n_iter=1500, TSNE_learning_rate2=300):
    markers = ['x', 'o']
    cols = spectra['marked'].apply(lambda x: colors[x])
    col = spectra['marked']

    plt.subplots(figsize=(18, 6))
    plt.subplot(131)
    plt.title("PCA")
    pca = PCA(n_components=2, random_state=42)
    spectra_2D = pca.fit_transform(spectra[spectrum_columns])
    for i in range(len(col)):
        plt.scatter(spectra_2D[i, 0], spectra_2D[i, 1], c=cols[i], alpha=0.5, marker=markers[col[i]], s=markersizes[col[i]])

    plt.subplot(132)
    plt.title("TSNE, Euclidean distance")
    tsne = TSNE(n_components=2, random_state=42, learning_rate=TSNE_learning_rate, n_iter=TSNE_n_iter)
    spectra_2D = tsne.fit_transform(spectra[spectrum_columns])
    for i in range(len(col)):
        plt.scatter(spectra_2D[i, 0], spectra_2D[i, 1], c=cols[i], alpha=0.5, marker=markers[col[i]], s=markersizes[col[i]])

    plt.subplot(133)
    plt.title("TSNE, Chosen distance")
    tsne = TSNE(n_components=2, random_state=42, metric="precomputed", learning_rate=TSNE_learning_rate2)
    spectra_2D = tsne.fit_transform(spectra_distances)
    for i in range(len(col)):
        plt.scatter(spectra_2D[i, 0], spectra_2D[i, 1], c=cols[i], alpha=0.5, marker=markers[col[i]], s=markersizes[col[i]])
    plt.show()
    
    
    # visualization - tsne with chosen distance
    print('Clustering')
    plt.subplots(figsize=(18, 18))

    for i, eps in enumerate(eps_l):
        dbscan = cluster.DBSCAN(eps=eps, min_samples=4)
        cluster_labels = dbscan.fit_predict(spectra_distances)

        plt.subplot(3, 3, i + 1)
        cols = [colors[l] for l in cluster_labels]
        plt.title("cluster labels (eps = {:.2})".format(eps))
        for j in range(len(col)):
            plt.scatter(spectra_2D[j, 0], spectra_2D[j, 1], c=cols[j], alpha=0.5, marker=markers[col[j]], s=markersizes[col[j]])

    plt.show()
    return spectra_2D
    
    
def print_clusters_structure_DBSCAN(spectra, spectrum_columns, other_names, spectra_distances, eps, colors, spectra_2D):
    
    dbscan = cluster.DBSCAN(eps=eps, min_samples=2)
    cluster_labels = dbscan.fit_predict(spectra_distances)
    centers = []
    col = spectra['marked']
    plt.figure(figsize=(18, 6))
    cols = [colors[l] for l in cluster_labels]
    col_name = 'DBSCAN_clusters_eps={}'.format(eps)
    spectra[col_name] = cluster_labels
    plt.title(col_name)
    for i in range(len(col)):
        plt.scatter(spectra_2D[i, 0], spectra_2D[i, 1], c=cols[i], alpha=0.5, marker=markers[col[i]], s=markersizes[col[i]])
    plt.show()
    list_spectra_clusters = []
    display(spectra[spectra['marked'] == 1])
    
    for i in set(cluster_labels):
        list_spectra_clusters.append(spectra[spectra[col_name] == i][other_names])
        centers.append(spectra.loc[spectra[col_name] == i, spectrum_columns].mean(axis=0).values)
    centers = pd.DataFrame(np.column_stack(centers).T, columns=spectrum_columns)
    
    del spectra['DBSCAN_clusters_eps={}'.format(eps)]

    # centroids
    return list_spectra_clusters, centers
    # plt.figure(figsize=(20, 10))
    # sns.heatmap(centers, vmin=0, vmax=1)
    # plt.show()

def show_all_spectra(spectra_columns, spectra):
    plt.clf()
    fig, ax = plt.subplots(figsize=(16,6))
    marked_spectra = spectra
    
    for i in range(0,len(marked_spectra)):
        table = marked_spectra[spectra_columns].iloc[i]
        values = table.values

        spectra_points_number = len(spectra_columns)
        plt.xticks(())
        #if i == 0:
        #    plt.setp([plt], title='Normalized columns')
        if i == 3:
            #plt.set_xlabel('optics - gev - tev')
            plt.xticks(range(spectra_points_number), spectra_columns,rotation=45)
        ax.plot(range(len(values)), marked_spectra[spectra_columns].iloc[i].values, label=marked_spectra['gev_1FGL_Name'].values[i])
        legend = ax.legend(loc = 'upper right')   
    plt.show()
