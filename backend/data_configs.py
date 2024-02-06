# Dataset configuration
# Each dataset is a dictionary with the following keys:
# - filepath (str): path to the h5ad file
# - name (str): name of the dataset
# - hd_metric (str): metric to use for the high-dimensional space [euclidean, angular]
# - description (str): description of the dataset
# - hd_data_key (str): key to access the high-dimensional data in adata.obsm (use "X" for adata.X)
# - metadata_features (list[str], optional): list of metadata features to include in the dataset


dataset_configs = {
    "Mammoth": {
        "filepath": "../data/mammoth/mammoth.h5ad",
        "name": "Mammoth",
        "hd_metric": "euclidean",
        "description": "Mammoth dataset from Wang et al. (2021)",
        "hd_data_key": "X",
    },
    "Mouse CD45neg": {
        "filepath": "../data/mouse_fibro/CD45_quality.h5ad",
        "name": "Mouse CD45neg",
        "hd_metric": "angular",
        "description": "Mouse cells, CD45neg (from Guilliams M. et al. 2022)",
        "hd_data_key": "X_totalVI",
        "metadata_features": ["annot", "louvain", "type", "sampleName"],
    },
    "Mouse Fibroblasts": {
        "filepath": "../data/mouse_fibro/mouseFibroblasts_20231218_154655.h5ad",
        "name": "Mouse Fibroblasts",
        "hd_metric": "angular",
        "description": "Mouse cells, Fibroblasts subset (from Guilliams M. et al. 2022)",
        "hd_data_key": "X_totalVI",
        "metadata_features": ["annot", "louvain", "type", "sampleName"],
    },
    "Human Immune": {
        "filepath": "../data/immune/immune_with_embeddings_new.h5ad",
        "name": "Human Immune",
        "hd_metric": "euclidean",
        "description": "Human Immune, Luecken et al. (2022)",
        "hd_data_key": "X_pca",
        "metadata_features": [
            "celltype",
            "batch",
            "chemistry",
            "species",
            "study",
            "tissue",
        ],
    },
    # default dataset that should always work (data is included in the repo)
    "GaussLine": {
        "filepath": "../data/gauss_line/gauss_line.h5ad",
        "name": "GaussLine",
        "hd_metric": "euclidean",
        "description": "Gaussian clusters shifted along a line from Böhm et al. (2022)",
        "hd_data_key": "X",
    },
}
