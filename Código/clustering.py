from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans, Birch, BirchDTW
from tslearn.clustering import TimeSeriesKMeans
from skfuzzy import cmeans
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

## Reading files ##

def read_region(roi):
    original = pd.read_csv("D:\\Projects\\GEE\\Time Series\\LS{}TS.csv".format(roi))
    original = original.dropna(how='all', axis=1)
    return original

def load_region(original, cc_value):    
    filtered = filter_cloudcover(original, cc_value)
    different_pixels = int(len(filtered.index) / filtered['Date'].nunique())
    
    pivot = filtered.pivot(index = 'Pixel', columns = 'Date', values = 'NDVI')
    pivot = pivot.dropna(how = 'any', axis = 1)
    pixels = filtered[:different_pixels]
    pixels = pixels[['Pixel']]
    return pivot, pixels

def filter_cloudcover(df, value):
    return df[df['CLOUD_COVER_REGION'] <= value]

###################################################

## Auxiliar Methods ##

def add_computed_typed(method, assignments, original):
    dataFrame = pd.DataFrame({'Pixel': [i for i in range(len(assignments))], 'cluster':assignments})
    dic = dataFrame.set_index('Pixel').T.to_dict('index').get('cluster')
    original['{} Type'.format(method)] = original['Pixel'].map(dic)
    return original
    
def add_computed_typed_pixels(method, assignments, pixels):
    pixels['{}'.format(method)] = assignments
    return pixels

def add_metrics_typed(method, metrics, assignments, pivot):
    if not(assignments is None) : 
        silhouette = silhouette_score(pivot, assignments, sample_size = 10000)
        metrics.loc[['Silhouette'], method] = silhouette
    return metrics

##################################################

## Best Parameter Calculation ##
    
def plot_KMeans_elbow(pivot):
    maxK = 2
    maxSil = -1
    for i in range(2, 20):
        kmeans = KMeans(n_clusters = i, random_state = 42)
        labels = kmeans.fit_predict(pivot)
        silhouette = silhouette_score(pivot, labels, sample_size = 10000)
        print("For n_clusters =", i,
          "The average silhouette_score is :", silhouette,
          "SSE is :", kmeans.inertia_)
        if silhouette > maxSil:
            maxSil = silhouette
            maxK = i
    return maxK
   
def plot_KMeansDTW_elbow(pivot):
    maxK = 2
    maxSil = -1
    for i in range(2, 20):
        kmeans = TimeSeriesKMeans(n_clusters=i, metric="dtw", random_state=42)
        labels = kmeans.fit_predict(pivot)
        silhouette = silhouette_score(pivot, labels, sample_size = 10000)
        print("For n_clusters =", i,
          "The average silhouette_score is :", silhouette,
          "SSE is :", kmeans.inertia_)
        if silhouette > maxSil:
            maxSil = silhouette
            maxK = i
    return maxK

def plot_birch_elbow(pivot):
    maxK = 2
    maxSil = -1
    for i in range(50, 1001, 25):
        birch = Birch(branching_factor=i, n_clusters=5, threshold=0.5, compute_labels=True)
        labels = birch.fit_predict(pivot)
        silhouette = silhouette_score(pivot, labels, sample_size = 10000)
        print("For n_clusters =", i,
          "The average silhouette_score is :", silhouette)
        if silhouette > maxSil:
            maxSil = silhouette
            maxK = i
    return maxK

def plot_birchDTW_elbow(pivot):
    maxK = 2
    maxSil = -1
    for i in range(50, 1001, 25):
        birch = BirchDTW(branching_factor=i, n_clusters=5, threshold=0.5, compute_labels=True)
        labels = birch.fit_predict(pivot)
        silhouette = silhouette_score(pivot, labels, sample_size = 10000)
        print("For n_clusters =", i,
          "The average silhouette_score is :", silhouette)
        if silhouette > maxSil:
            maxSil = silhouette
            maxK = i
    return maxK

def plot_fuzzyCMeans_elbow(pivot):
    maxK = 2
    maxSil = -1
    for i in range(2, 20):
        cntr, u, u0, d, jm, p, fpc = cmeans(pivot, i, 2, error=0.005, maxiter=1000, init=None)     
        silhouette = silhouette_score(pivot, np.argmax(u, axis = 0), sample_size = 10000)
        print("For n_clusters =", i,
          "The average silhouette_score is :", silhouette,
          "FPC is :", fpc)
        if silhouette > maxSil:
            maxSil = silhouette
            maxK = i
    return maxK 
    
def plot_fuzzyCMeansDTW_elbow(pivot):
    maxK = 2
    maxSil = -1
    for i in range(2, 20):
        cntr, u, u0, d, jm, p, fpc = cmeans(pivot, i, 2, metric = 'dtw', error=0.005, maxiter=1000, init=None)     
        silhouette = silhouette_score(pivot, np.argmax(u, axis = 0), sample_size = 10000)
        print("For n_clusters =", i,
          "The average silhouette_score is :", silhouette,
          "FPC is :", fpc)
        if silhouette > maxSil:
            maxSil = silhouette
            maxK = i
    return maxK

#################################################################

## Clustering ## 
     
def cl_KMeans(pivot, original, pixels, metrics, k, cc):
    kmeans = KMeans(n_clusters = k, random_state = 42).fit(pivot)   
    method = 'K_{0}_{1}'.format(k, cc*10)
    assignments = kmeans.labels_
    prediction = kmeans.predict(pivot)
    metrics.loc[['SSE'],method] = kmeans.inertia_
    return add_computed_typed(method, assignments, original), add_computed_typed_pixels(method, assignments, pixels), add_metrics_typed(method, metrics, prediction, pivot)
        
def cl_KMeansDTW(pivot, original, pixels, metrics, k, cc):
    kmeans = TimeSeriesKMeans(n_clusters = k, metric="dtw", random_state = 42)
    method = 'KDTW_{0}_{1}'.format(k, cc*10)
    newpivot = pivot.values
    assignments = kmeans.fit_predict(newpivot)  
    return add_computed_typed(method, assignments, original), add_computed_typed_pixels(method, assignments, pixels), add_metrics_typed(method, metrics, assignments, pivot)

def cl_birch(pivot, original, pixels, metrics, bf, k, cc):
    birch = Birch(branching_factor= bf, n_clusters=k, threshold=0.5, compute_labels=True).fit(pivot)
    method = 'B_{0}_{1}'.format(k, cc*10)
    assignments = birch.labels_
    prediction = birch.predict(pivot)
    return add_computed_typed(method, assignments, original) , add_computed_typed_pixels(method, assignments, pixels), add_metrics_typed(method, metrics, prediction, pivot)
      
def cl_birchDTW(pivot, original, pixels, metrics, bf, k, cc):
    birch = BirchDTW(branching_factor= bf, n_clusters=k, threshold=0.5, compute_labels=True).fit(pivot)
    method = 'B_{0}_{1}'.format(k, cc*10)
    assignments = birch.labels_
    prediction = birch.predict(pivot)
    return add_computed_typed(method, assignments, original) , add_computed_typed_pixels(method, assignments, pixels), add_metrics_typed(method, metrics, prediction, pivot)

def cl_fuzzyCMeans(pivot, original, pixels, metrics, k, cc):
    cntr, u, u0, d, jm, p, fpc = cmeans(pivot.T, k, 2, error=0.005, maxiter=1000, init=None)
    method = 'F_{0}_{1}'.format(k, cc*10)
    assignments = np.argmax(u, axis = 0) 
    metrics.loc[['FPC'],method] = fpc 
    return add_computed_typed(method, assignments, original), add_computed_typed_pixels(method, assignments, pixels), add_metrics_typed(method, metrics, assignments, pivot)

def cl_fuzzyCMeansDTW(pivot, original, pixels, metrics, k, cc):
    cntr, u, u0, d, jm, p, fpc = cmeans(pivot.T, k, 2, metric = 'dtw', error=0.005, maxiter=1000, init=None)
    method = 'F_{0}_{1}'.format(k, cc*10)
    assignments = np.argmax(u, axis = 0) 
    metrics.loc[['FPC'],method] = fpc 
    return add_computed_typed(method, assignments, original), add_computed_typed_pixels(method, assignments, pixels), add_metrics_typed(method, metrics, assignments, pivot)

###########################################################

cc_values = [0, 2.5, 5, 10]
regions = ['Santar', 'SCTSCL', 'BenRib', 'BGL']   
metrics = pd.DataFrame({'Metric': ['Silhouette', 'FPC', 'SSE']})
metrics.set_index("Metric", inplace=True)
for roi in regions:
    print('Region =', roi)
    original = read_region(roi)
    result = original
    pixels_all = pd.DataFrame()
    for cc_value in cc_values:
        print('Cloud Cover = ', cc_value)
        pivot, pixels = load_region(original, cc_value)
        
        #best_bf= plot_birch_elbow(pivot)
        #best_k = plot_KMeans_elbow(pivot)
        #best_k = plot_KMeansDTW_elbow(pivot)
        #best_k = plot_fuzzyCMeans_elbow(pivot.T)
        #best_k = plot_fuzzyCMeansDTW_elbow(pivot.T)
        
        for k in range(3,8,2):
            print('K =', k)
            result, pixels, metrics = cl_KMeans(pivot, result, pixels, metrics, k, cc_value)
            result, pixels, metrics = cl_KMeansDTW(pivot, original, pixels, metrics, k, cc_value)  
            result, pixels, metrics = cl_fuzzyCMeans(pivot, result, pixels, metrics, k, cc_value) 
            result, pixels, metrics = cl_fuzzyCMeansDTW(pivot, result, pixels, metrics, k, cc_value) 
            result, pixels, metrics = cl_birch(pivot, result, pixels, metrics, 75, k, cc_value)
            result, pixels, metrics = cl_birchDTW(pivot, result, pixels, metrics, 75, k, cc_value)            
        pixels_all = pd.concat([pixels_all, pixels], axis = 1)
        
    #All class Attributions and Extras
    result.to_csv('Results_{}.csv'.format(roi))
    
    #All class Attributions
    pixels_all.to_csv('Pixel_{}.csv'.format(roi))
    
    #All class Attributions Metrics
    metrics.to_csv('Metrics_{}.csv'.format(roi))