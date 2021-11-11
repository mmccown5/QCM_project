import numpy as np
import pandas as pd
import logging
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

from scipy.sparse import csc_matrix

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

def method(input_train_mod1, input_train_mod2, input_test_mod1):
    '''Basic user-implemented method'''
    logging.info('TO DO: Calculate something useful...')
    
    ###
    
    pred = baseline_mean(input_train_mod1, input_train_mod2, input_test_mod1)
    
    y_pred = pred.X
    #X should be csc_matrix
    x2 = csc_matrix(y_pred)
    pred.X = x2
        
    pred.uns["method"] = "basic_beans"
    pred.uns['dataset_id'] = input_train_mod1.uns['dataset_id']
    
    
    
    return pred

def starter_method(input_train_mod1, input_train_mod2, input_test_mod1):
    '''Method from starter kit'''
    # Do PCA on the input data
    logging.info('Performing dimensionality reduction on modality 1 values...')
    embedder_mod1 = TruncatedSVD(n_components=50)
    mod1_pca = embedder_mod1.fit_transform(input_train.X)

    logging.info('Performing dimensionality reduction on modality 2 values...')
    embedder_mod2 = TruncatedSVD(n_components=50)
    mod2_pca = embedder_mod2.fit_transform(input_train_mod2.X)

    # split dimred back up
    X_train = mod1_pca[input_train.obs['group'] == 'train']
    X_test = mod1_pca[input_train.obs['group'] == 'test']
    y_train = mod2_pca

    assert len(X_train) + len(X_test) == len(mod1_pca)

    # Get all responses of the training data set to fit the
    # KNN regressor later on.
    #
    # Make sure to use `toarray()` because the output might
    # be sparse and `KNeighborsRegressor` cannot handle it.

    logging.info('Running Linear regression...')

    reg = LinearRegression()

    # Train the model on the PCA reduced modality 1 and 2 data
    reg.fit(X_train, y_train)
    y_pred = reg.predict(X_test)

    # Project the predictions back to the modality 2 feature space
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
            'method_id': 'starter_kit'
        },
    )
    
    return adata


def baseline_mean(input_train_mod1, input_train_mod2, input_test_mod1):
    '''Dummy method that predicts mean(input_train_mod2) for all cells'''
    logging.info('Calculate mean of the training data modality 2...')
    y_pred = np.repeat(input_train_mod2.X.mean(axis=0).reshape(-1,1).T, input_test_mod1.shape[0], axis=0)
    
    # Prepare the ouput data object
    pred_test_mod2 = ad.AnnData(
        X=y_pred,
        obs=input_test_mod1.obs,
        var=input_train_mod2.var,
    )
    
    pred_test_mod2.uns["method"] = "mean"
    pred_test_mod2.uns['dataset_id'] = input_train_mod1.uns['dataset_id']

    return pred_test_mod2


def baseline_linear(input_train_mod1, input_train_mod2, input_test_mod1):
    '''Baseline method training a linear regressor on the input data'''
    input_mod1 = ad.concat(
        {"train": input_train_mod1, "test": input_test_mod1},
        axis=0,
        join="outer",
        label="group",
        fill_value=0,
        index_unique="-", 
    )
    
    # Do PCA on the input data
    logging.info('Performing dimensionality reduction on modality 1 values...')
    embedder_mod1 = TruncatedSVD(n_components=50)
    mod1_pca = embedder_mod1.fit_transform(input_mod1.X)
    

    
    # split dimred mod 1 back up for training
    X_train = mod1_pca[input_mod1.obs['group'] == 'train']
    X_test = mod1_pca[input_mod1.obs['group'] == 'test']
    y_train = input_train_mod2.X.toarray()
    
    assert len(X_train) + len(X_test) == len(mod1_pca)
    
    logging.info('Running Linear regression...')
    
    reg = LinearRegression()
    
    # Train the model on the PCA reduced modality 1 and 2 data
    reg.fit(X_train, y_train)
    y_pred = reg.predict(X_test)
    
    # Project the predictions back to the modality 2 feature space
    
    pred_test_mod2 = ad.AnnData(
        X = y_pred,
        obs = input_test_mod1.obs,
        var = input_train_mod2.var,
    
    )
    
    # Add the name of the method to the result
    pred_test_mod2.uns["method"] = "linear"
    pred_test_mod2.uns['dataset_id'] = input_train_mod1.uns['dataset_id']
    
    return pred_test_mod2
