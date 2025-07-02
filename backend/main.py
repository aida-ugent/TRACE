import argparse
import os
import sys
from contextlib import asynccontextmanager
from typing import List

import continuous_palettes
import neighbors
import numpy as np
import pandas as pd
import uvicorn
from dataset import Dataset
from fastapi import Body, FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pandas.api.types import CategoricalDtype
from pydantic import BaseModel
from utils import get_available_datasets
import time
API_PORT = 8000
dataset = None

if os.path.isfile("./data_configs.yaml"):
    dataset_configs = get_available_datasets("./data_configs.yaml")
else:
    print(
        "Error: data_configs.yaml not found. Did you forget to create your own "
        "data config using the data_configs.yaml.template?"
    )
    sys.exit(1)


@asynccontextmanager
async def lifespan(app: FastAPI):
    global dataset
    print("Starting server")
    dataset = Dataset(name="tmp", hd_data=np.zeros((10, 10)))
    yield
    if dataset is not None:
        dataset.cleanup()
    print("Shutdown server")


app = FastAPI(title="TRACE: Interactive data visualization", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://localhost:8000",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class NeighborRequest(BaseModel):
    k: int
    points: List[int]
    hd_metric: str


class PointSelection(BaseModel):
    points: List[int]
    selection_name: str


class clusterSelection(BaseModel):
    selectionA: List[int]
    selectionB: List[int]


class FeatureList(BaseModel):
    feature_list: List[str]


@app.get("/backend/datasetOptions")
async def getDatasetNames():
    return {"result": list(dataset_configs.keys())}


@app.get("/backend/loadDataset")
def loadDataset(datasetName: str):
    """
    Load a dataset by name.

    Args:
        datasetName (str): The name of the dataset to load.

    Returns:
        dict: A dictionary containing the result, hd_metric, and dataset_info
            to be displayed in the dashboard.
    """
    global dataset

    if dataset is not None and dataset.name != datasetName:
        dataset.cleanup()

        if datasetName in dataset_configs.keys():
            datasetNameConfig = dataset_configs[datasetName].pop("name", datasetName)
            try:
                dataset = Dataset(
                    name=datasetNameConfig, **dataset_configs[datasetName]
                )
                print(f"Loaded dataset {datasetName}")
            except FileNotFoundError as e:
                print(f"Error loading dataset {datasetName}: {e}")
                print(f"Now loading GaussianLine")
                datasetName = "GaussLine"
                dataset = Dataset(name=datasetName, **dataset_configs["GaussLine"])
        else:
            raise HTTPException(
                status_code=500, detail=f"Dataset {datasetName} not found"
            )
    if dataset.adata.n_vars > 4000:
        print(
            f"Warning: dataset has {dataset.adata.n_vars} features. Limiting to 4000."
        )

    point_color_options = {
        "quality": dataset.get_quality_features(),
        "metadata": dataset.get_metadata_features() + ["HD distances", "none"],
        "features": dataset.adata.var.index.tolist()[
            0 : 4000 if dataset.adata.n_vars > 4000 else dataset.adata.n_vars
        ],
    }
    return {
        "dataset_name": datasetName,
        "metric_options": dataset.get_hd_metric_options(),
        "max_neighbors": dataset.get_precomputed_neighbors_maxK(),
        "hd_metric": dataset.hd_metric,
        "dataset_info": dataset.description,
        "embedding_options": dataset.get_embedding_options(),
        "point_color_options": point_color_options,
        "exclus": dataset.get_exclus_results(),
    }


@app.get("/backend/metadataNames")
async def getMetadataNames():
    return {"result": dataset.get_metadata_features()}


@app.get("/backend/featureNames")
async def getFeatureNames():
    return {"result": dataset.adata.var.index.tolist()}


@app.get("/backend/landmarkPoints")
async def getLandmarkPoints():
    return {"result": dataset.adata.uns["landmark_indices"].tolist()}


@app.get("/backend/featureValues/{fname}")
async def getFeatureValues(fname: str):
    if fname in dataset.adata.var.index:
        expression = dataset.adata.obs_vector(fname).tolist()
        return {"result": expression}
    else:
        raise ValueError(f"feature {fname} not in variable names of dataset")


@app.post("/backend/metadataFeatures")
async def getMetadataFeatures(featureList: FeatureList = Body(...)):
    """
    Get feature values for a list of features.

    Args:
        feature_list (List[str]): A list of feature names.

    Returns:
        dict: A dictionary containing the feature values for each feature in the list.
    """
    features = featureList.feature_list
    result = {}
    for feature in features:
        # feature is data dimension
        if feature in dataset.adata.var.index:
            expression = dataset.adata.obs_vector(feature).tolist()
            result[feature] = expression
        # feature is metadata
        elif feature in dataset.adata.obs_keys():
            result[feature] = dataset.adata.obs[feature].tolist()
        else:
            raise ValueError(
                f"feature {feature} not in variable names or data dimensions."
            )
    return result


@app.get("/backend/embedding")
async def getEmbedding(embName: str):
    """Retrieve embedding from anndata.obsm

    Args:
        embName (str): The name of the embedding.

    Returns:
        dict: A dictionary containing the x and y coordinates of the embedding.
    """
    if embName in dataset.adata.obsm_keys():
        result = {
            "x": dataset.adata.obsm[embName][:, 0].tolist(),
            "y": dataset.adata.obsm[embName][:, 1].tolist(),
        }
    else:
        print(f"Embedding {embName} not found in obsm_keys {dataset.adata.obsm_keys()}")
        result = {
            "x": np.zeros((dataset.adata.n_obs,), dtype=int).tolist(),
            "y": np.zeros((dataset.adata.n_obs,), dtype=int).tolist(),
        }
    return result


@app.get("/backend/pointColor/{fname}")
async def getPointColors(
    fname: str, embeddingName: str, selectedPoint: int = None, hdMetric: str = None
):
    """Retrieve values for metadata or feature, normalize,
    and compute colormaps.

    Args:
        fname (str): name of the feature to be found in the anndata.obs_keys,
            anndata.var.index, or anndata.obsm[embeddingName].colnames()
        embeddingName (str): name of embedding when retrieving
            embedding specific values (such as quality scores).
            Defaults to None.

    Returns:
        dict: Dictionary with the following keys:
            "values": list,
            "colorMap": dict with keys 'ticks' and 'colors' which are two lists of values
            "type": continuous | categorical,
    """
    scale_continuous = True
    range = None
    fgroup = None
    categorical_dtypes = ["object", "category", "bool"]

    if fname == "HD distances" and (selectedPoint is not None and hdMetric is not None):
        # compute HD distances for selected point
        fvalues = dataset.get_HD_landmark_distances(selectedPoint, hd_metric=hdMetric)
        ftype = "continuous"
        fgroup = "metadata"
        colors = continuous_palettes.palettes["viridis"]

    # metadata
    elif fname in dataset.adata.obs_keys():
        fvalues = dataset.adata.obs[fname]
        fdtype = dataset.adata.obs[fname].dtype.name
        fgroup = "metadata"

        if fdtype in categorical_dtypes or (
            "int" in fdtype and len(dataset.adata.obs[fname].unique()) < 50
        ):
            ftype = "categorical"
            # do sorting by frequency only for categories
            if fdtype in categorical_dtypes:
                value_counts = list(
                    fvalues.value_counts(ascending=False, sort=True).keys()
                )
            else:
                value_counts = sorted(pd.unique(fvalues).tolist())
            max_colors = 19

            if len(value_counts) > max_colors:
                # show the first 20 categories in different colors, the rest as grey
                other_cat_name = f"other ({len(value_counts) - max_colors})"
                categories = value_counts[:max_colors] + [other_cat_name]
                cat_dtype = CategoricalDtype(categories=categories)
                fvalues = fvalues.astype(cat_dtype)
                fvalues[~fvalues.isin(value_counts[:max_colors])] = other_cat_name
                encoded_fvalues = pd.Series(fvalues.values.tolist()).astype(cat_dtype)
                encoded_fvalues = encoded_fvalues.cat.codes
                colorMap = dataset.get_category_colors(fname, categories=categories)
                colorMap["colors"][-1] = "#777777"
            else:
                cat_dtype = CategoricalDtype(categories=value_counts)
                encoded_fvalues = pd.Series(fvalues.values.tolist()).astype(cat_dtype)
                encoded_fvalues = encoded_fvalues.cat.codes
                colorMap = dataset.get_category_colors(fname, categories=value_counts)
                range = [0, 1]
        else:
            ftype = "continuous"
            colors = dataset.get_continuous_colors(fname)
            if colors is None:
                colors = continuous_palettes.palettes["viridis"]

    # features
    elif fname in dataset.adata.var.index:
        if "norm_genes" in dataset.adata.layers:
            layer = "norm_genes"
        else:
            layer = "X"
        fvalues = dataset.adata.obs_vector(fname, layer=layer)
        colors = continuous_palettes.palettes["viridis"]
        ftype = "continuous"
        fgroup = "features"

    # quality score for specific dataset
    elif (
        embeddingName in dataset.adata.uns_keys()
        and fname in dataset.adata.uns[embeddingName]["quality"]
    ):
        fvalues = dataset.adata.uns[embeddingName]["quality"][fname]

        # if correlation scale from [-1, 1] to [0,1]
        if "corr" in fname:
            range = [-1, 1]
            colors = continuous_palettes.palettes["diverging_purple_green"]
        else:
            # other quality features are within [0,1] (neighborhood preservation)
            range = [0, 1]
            colors = continuous_palettes.palettes["viridis"]

        # only scale if outside of [0,1] range
        if range[0] >= 0 and range[1] <= 1:
            scale_continuous = False

        ftype = "continuous"
        fgroup = "quality"
        print(f"quality score {fname} has mean {fvalues.mean()}")

    # fallback
    else:
        print(f"{fname} is not in obs_names {dataset.adata.obs_keys()}")
        fvalues = encoded_fvalues = np.zeros((dataset.adata.n_obs,), dtype=int)
        colorMap = {"ticks": ["none"], "colors": ["#444444"]}
        ftype = "categorical"
        fgroup = "metadata"

    # encode values
    if ftype == "continuous":
        if isinstance(fvalues, pd.Series):
            fvalues = fvalues.values

        # check for NaN values
        if np.isnan(np.sum(fvalues)):
            isNaN = np.isnan(fvalues)
            NaNmean = np.nanmean(fvalues)
            fvalues[isNaN] = NaNmean
            print(
                f"{fname} contains {sum(isNaN)} NaN values. Replacing with mean {NaNmean}."
            )
        if range is None:
            range = [float(fvalues.min()), float(fvalues.max())]
        ticks = np.linspace(
            start=range[0], stop=range[1], num=len(colors), endpoint=True
        )
        ticks = [f"{x:.2}" for x in ticks]
        colorMap = {"ticks": ticks, "colors": colors}
        if scale_continuous:
            if range[1] - range[0] == 0:
                encoded_fvalues = np.zeros((dataset.adata.n_obs,), dtype=int)
            else:
                encoded_fvalues = (fvalues - range[0]) / (range[1] - range[0])
        else:
            encoded_fvalues = fvalues

    res = {
        "values": fvalues.tolist(),
        "encoded_values": encoded_fvalues.tolist(),
        "colorMap": colorMap,
        "type": ftype,
        "group": fgroup,
    }
    return res


@app.post("/backend/precomputeAllNeighbors")
async def precomputeAllNeighbors(maxK: int, hd_metric: str):
    if dataset is not None and (
        dataset.adata.uns["hd_neighbors"] is None
        or hd_metric not in dataset.adata.uns["hd_neighbors"].keys()
        or maxK > dataset.adata.uns["hd_neighbors"][hd_metric].shape[1]
    ):
        dataset.precompute_HD_neighbors(maxK, hd_metric)
    return {"result": "success"}


@app.post("/backend/intrusions")
async def getIntrusions(item: NeighborRequest = Body(...)):
    """
    Get intrusions based on the given NeighborRequest. Intrusions are defined
    as points that are in the selection item.points, but not in the k nearest
    neighbors of these points.

    Args:
        item (NeighborRequest): The NeighborRequest object containing the parameters for intrusion detection.

    Returns:
        dict: A dictionary containing the computed results.
            - "result": A list of computed intrusions.
            - "binary": A list representing the binary representation of the intrusions.
    """
    neighbors = dataset.get_HD_neighbors(
        k=item.k, metric=item.hd_metric, indices=np.asarray(item.points)
    )
    neighbors = neighbors.flatten()
    binary = np.zeros((dataset.adata.n_obs,), dtype=int)
    binary[item.points] = 1
    binary[neighbors] = 0
    return {"result": binary.nonzero()[0].tolist(), "binary": binary.tolist()}


@app.post("/backend/savePointSelection")
async def savePointSelection(item: PointSelection = Body(...)):
    """
    Save a point selection to the dataset.
    """
    dataset.save_user_annotation(item.points, item.selection_name)
    return {"result": "success"}


@app.post("/backend/computeHDNeighbors")
async def computeHDNeighbors(item: NeighborRequest):
    """
    Compute high-dimensional neighbors based on the given parameters.

    Args:
        item (NeighborRequest): The request object containing the parameters.

    Returns:
        dict: A dictionary containing the computed results.
            - "result": A list of computed neighbors.
            - "binary": A list representing the binary representation of the neighbors.
    """
    neighbors = dataset.get_HD_neighbors(
        k=item.k, metric=item.hd_metric, indices=np.asarray(item.points)
    )
    neighbors = neighbors.flatten()
    binary = np.zeros((dataset.adata.n_obs,), dtype=int)
    binary[neighbors] = 1
    return {"result": binary.nonzero()[0].tolist(), "binary": binary.tolist()}


@app.post("/backend/computeReverseHDNeighbors")
async def computeReverseHDNeighbors(item: NeighborRequest):
    """
    Compute reverse HD neighbors of a point selection. Reverse HD neighbors are
    defined as points that have at least one of the points in the selection as
    one of their k nearest neighbors.

    Args:
        item (NeighborRequest): The NeighborRequest object containing the parameters.

    Returns:
        dict: A dictionary containing the indices and a binary encoding of the
            reverse neighbors of the selection.
    """
    neighbors = dataset.get_HD_neighbors(
        k=item.k, metric=item.hd_metric, indices=np.arange(dataset.adata.n_obs)
    )
    query_points = np.asarray(item.points)

    # having the second array as the smaller one is faster
    if len(query_points) > item.k:
        bool_result = np.apply_along_axis(
            lambda a: np.isin(query_points, a, assume_unique=True).any(),
            1,
            neighbors,
        )
    else:
        bool_result = np.apply_along_axis(
            lambda a: np.isin(a, query_points, assume_unique=True).any(),
            1,
            neighbors,
        )

    binary = bool_result.astype(int)
    return {"result": binary.nonzero()[0].tolist(), "binary": binary.tolist()}


@app.get("/backend/getUnstablePoints")
async def getUnstablePoints(
    embNameA: str, embNameB: str, maxFraction: float, k: int = 50
):
    """Computing points that are moving a lot between two embeddings."""
    if embNameA in dataset.adata.obsm_keys() and embNameB in dataset.adata.obsm_keys():
        return {
            "result": neighbors.unstable_points(
                dataset.adata.obsm[embNameA],
                dataset.adata.obsm[embNameB],
                maxFraction,
                k,
            )
        }
    else:
        return {"result": []}


@app.post("/backend/explainCluster")
async def explainCluster(item: PointSelection = Body(...)):
    indices = item.points
    method = item.selection_name
    start = time.time()
    features, higher_mean = dataset.explain_cluster(indices, method=method)
    print(
        f"Explaining cluster {method} with {len(indices)} points took {time.time() - start:.2f} seconds"
    )
    return {"features": features.tolist(), "higher_mean": higher_mean.tolist()}


@app.post("/backend/compareClusters")
async def compareClusters(item: clusterSelection = Body(...)):
    features, higher_mean = dataset.compareClusters(item.selectionA, item.selectionB)
    return {"features": features.tolist(), "higher_mean": higher_mean.tolist()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=API_PORT)
    parser.add_argument("--host", type=str, default="127.0.0.1")
    args = parser.parse_args()

    uvicorn.run("main:app", host=args.host, port=args.port, log_level="info")
