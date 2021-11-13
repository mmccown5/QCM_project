import numpy as np
import pandas as pd
import logging
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import mean_squared_error

from scipy.sparse import csc_matrix

#import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

def method(input_train_mod1, input_train_mod2, input_test_mod1):
    '''Basic user-implemented method'''
    
    
    logging.info('Performing dimensionality reduction on modality 1 values...')    
    input_mod1 = ad.concat(
        {"train": input_train_mod1, "test": input_test_mod1},
        axis=0,
        join="outer",
        label="group",
        fill_value=0,
        index_unique="-"
    )
    
    embedder_mod1 = TruncatedSVD(n_components=50)
    mod1_pca = embedder_mod1.fit_transform(input_mod1.X)

    logging.info('Performing dimensionality reduction on modality 2 values...')
    embedder_mod2 = TruncatedSVD(n_components=50)
    mod2_pca = embedder_mod2.fit_transform(input_train_mod2.X)

    # split dimred back up
    X_train = mod1_pca[input_mod1.obs['group'] == 'train']
    X_test = mod1_pca[input_mod1.obs['group'] == 'test']
    y_train = mod2_pca

    assert len(X_train) + len(X_test) == len(mod1_pca)

    logging.info('Running K nearest neigbors...')

    # Train the model on the PCA reduced modality 1 and 2 data
    neigh = KNeighborsRegressor(n_neighbors=2)
    neigh.fit(X_train, mod2_pca)

    # Project the predictions back to the modality 2 feature space
    y_pred = neigh.predict(X_test)
    y_pred = y_pred @ embedder_mod2.components_

    # Store as sparse matrix to be efficient. Note that this might require
    # different classifiers/embedders before-hand. Not every class is able
    # to support such data structures.
    y_pred = csc_matrix(y_pred)

    adata = ad.AnnData(
        X=y_pred,
        obs=input_test_mod1.obs,
        var=input_train_mod2.var,
        uns={
            'dataset_id': input_train_mod1.uns['dataset_id'],
            'method_id': 'basic_beans'
        },
    )

    return adata