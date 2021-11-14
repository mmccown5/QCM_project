# Dependencies:
# pip: scikit-learn, anndata, scanpy
#
# Modified from the Python starter kit for the NeurIPS 2021 Single-Cell Competition.
# Parts with `TODO` are supposed to be changed by you.
#
# More documentation:
#
# https://viash.io/docs/creating_components/python/

import logging
import anndata as ad
import sys

from scipy.sparse import csc_matrix

from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import KNeighborsRegressor
from statistics import sqrt
from math import ceil

#from beans import method # unit test wasn't finding it
def method(input_train_mod1, input_train_mod2, input_test_mod1, k=0):
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
    
    #by default, calculate k thusly:
    if k==0:
        N = X_train.shape[0]
        k = ceil(sqrt(N))

    logging.info('Running K nearest neigbors...')

    # Train the model on the SVD reduced modality 1 and 2 data
    neigh = KNeighborsRegressor(n_neighbors=k)
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


###
#  Return to script
logging.basicConfig(level=logging.INFO)

## VIASH START
# Anything within this block will be removed by `viash` and will be
# replaced with the parameters as specified in your config.vsh.yaml.
meta = { 'resources_dir': '.' }

par = {
    'input_train_mod1': 'sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod1.h5ad',
    'input_train_mod2': 'sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod2.h5ad',
    'input_test_mod1': 'sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod1.h5ad',
    'distance_method': 'minkowski',
    'output': 'output.h5ad',
    'n_pcs': 50,
}
## VIASH END

# TODO: change this to the name of your method
method_id = "basic_beans"

sys.path.append(meta['resources_dir'])

logging.info('Reading `h5ad` files...')
input_train_mod1 = ad.read_h5ad(par['input_train_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])
input_test_mod1 = ad.read_h5ad(par['input_test_mod1'])

# TODO: implement own method
adata = method(input_train_mod1, input_train_mod2, input_test_mod1)

adata.uns["method_id"] = method_id

logging.info('Storing annotated data...')
adata
adata.write_h5ad(par['output'], compression = "gzip")



