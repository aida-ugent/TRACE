from typing import List
from fastapi import FastAPI, Body, HTTPException
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
import numpy as np
from utils import get_available_datasets
from dataset import Dataset
from contextlib import asynccontextmanager
import pandas as pd
from pandas.api.types import CategoricalDtype
import continuous_palettes
import neighbors
import uvicorn
import argparse

API_PORT = 8000
dataset = None
dataset_configs = get_available_datasets("./data_configs.yaml")


@asynccontextmanager
async def lifespan(app: FastAPI):
    global dataset
    print("Starting server")
    dataset = Dataset(**dataset_configs[list(dataset_configs.keys())[0]])
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
            try:
                dataset = Dataset(**dataset_configs[datasetName])
                print(f"Loaded dataset {datasetName}")
            except FileNotFoundError as e:
                print(f"Error loading dataset {datasetName}: {e}")
                print(f"Now loading GaussianLine")
                dataset = Dataset(**dataset_configs["GaussLine"])
                datasetName = "GaussLine"
        else:
            raise HTTPException(
                status_code=500, detail=f"Dataset {datasetName} not found"
            )

    point_color_options = {
        "metadata": dataset.get_metadata_features() + ["HD distances", "none"],
        "quality": dataset.get_quality_features(),
        "features": dataset.adata.var.index.tolist(),
    }
    return {
        "dataset_name": datasetName,
        "metric_options": dataset.get_hd_metric_options(),
        "hd_metric": dataset.hd_metric,
        "dataset_info": dataset.description,
        "embedding_options": dataset.get_embedding_options(),
        "point_color_options": point_color_options,
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
async def getPointColors(fname: str, embeddingName: str, selectedPoint: int = None):
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
            "colorMap": dict with labels as keys and colors as values
            "type": continuous | categorical,
    """
    scale_continuous = True
    range = None
    fgroup = None

    if fname == "HD distances" and selectedPoint is not None:
        # compute HD distances for selected point
        fvalues = dataset.get_HD_landmark_distances(selectedPoint)
        ftype = "continuous"
        fgroup = "metadata"
        colors = continuous_palettes.palettes["viridis"]

    # metadata
    elif fname in dataset.adata.obs_keys():
        fvalues = dataset.adata.obs[fname]
        fgroup = "metadata"

        if (
            dataset.adata.obs[fname].dtype.name == "category"
            or dataset.adata.obs[fname].dtype.name == "object"
        ):
            ftype = "categorical"
            value_counts = list(fvalues.value_counts(ascending=False, sort=True).keys())

            if len(value_counts) > 30:
                print(f"{fname} contains too many unique values ({len(value_counts)})")
                fvalues = encoded_fvalues = np.zeros((dataset.adata.n_obs,), dtype=int)
                colors = {"too many values": "#444444"}
                ftype = "categorical"
                fgroup = "metadata"
            else:
                cat_dtype = CategoricalDtype(categories=value_counts)
                encoded_fvalues = pd.Series(fvalues.values.tolist()).astype(cat_dtype)
                encoded_fvalues = encoded_fvalues.cat.codes
                colors = dataset.get_category_colors(fname, categories=value_counts)
                range = [0, 1]
        else:
            ftype = "continuous"
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
        else:
            # other quality features are within [0,1] (neighborhood preservation)
            range = [0, 1]

        # only scale if outside of [0,1] range
        if range[0] >= 0 and range[1] <= 1:
            scale_continuous = False

        colors = continuous_palettes.palettes["viridis"]
        ftype = "continuous"
        fgroup = "quality"
        print(f"quality score {fname} has mean {fvalues.mean()}")

    # fallback
    else:
        print(f"{fname} is not in obs_names {dataset.adata.obs_keys()}")
        fvalues = encoded_fvalues = np.zeros((dataset.adata.n_obs,), dtype=int)
        colors = {"none": "#444444"}
        ftype = "categorical"
        fgroup = "metadata"

    # encode values
    if ftype == "continuous":
        # check for NaN values
        if np.isnan(np.sum(fvalues)):
            print(f"{fname} contains NaN values. Replacing with mean.")
            fvalues = fvalues.fillna(fvalues.mean())
        if range is None:
            range = [float(fvalues.min()), float(fvalues.max())]
        colorticks = np.linspace(
            start=range[0], stop=range[1], num=len(colors), endpoint=True
        )
        colorticks = [f"{x:.2}" for x in colorticks]
        colors = dict(zip(colorticks, colors))
        if scale_continuous:
            encoded_fvalues = (fvalues - range[0]) / (range[1] - range[0])
        else:
            encoded_fvalues = fvalues

    res = {
        "values": fvalues.tolist(),
        "encoded_values": encoded_fvalues.tolist(),
        "colorMap": colors,
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=API_PORT)
    parser.add_argument("--host", type=str, default="127.0.0.1")
    args = parser.parse_args()

    uvicorn.run("main:app", host=args.host, port=args.port, log_level="info")
