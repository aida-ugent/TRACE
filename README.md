# Pattern or Artifact? Interactively Exploring Embedding Quality with TRACE

TRACE<sup>[1](#trace)</sup> supports you in analyzing **global and local quality 🕵🏽‍♀️** of two-dimensional embeddings, based on [Regl-scatterplot](https://github.com/flekschas/regl-scatterplot)<sup>[2](#regl_citation)</sup> .


## Installation

<details>
<summary><b>OPTION 1: Using Docker 🐋</b></summary>

Make sure you have [Docker Compose](https://docs.docker.com/compose/install/) installed. Then build the container that includes the backend and frontend.

First you might need to install [docker-compose](https://docs.docker.com/compose/install/).

```bash
docker-compose build
docker-compose up
```
This will mount the /frontend, /backend, and /data directories into the repective containers. 

Open [http://localhost:3000](http://localhost:3000) with your browser to see the result.

</details>

<details>
<summary><b>OPTION 2: Without Docker</b></summary>

#### Required packages
**Backend**: Install the required python packages for the backend, tested with Python 3.11 from `backend/pip_requirements.txt` or `backend/conda_requirements.yaml`. 

**Frontend**: Install the packages in `frontend/package.json` using e.g. `npm install`.

First, start the backend within the right python evironment:
```bash
conda activate backend_env/
python main.py
# or
python -m uvicorn main:app --reload
```

Then start the frontend development server:
```bash
npm run dev
```
Open [http://localhost:3000](http://localhost:3000) with your browser to see the result.

</details>

## Data Preparation

The easiest way to load your data into TRACE is using the `Dataset` class to add embeddings and compute quality measures. This will create the necessary Anndata structure under the hood. Examples can be found in the notebooks of each dataset folder. 

```python
trace_data = Dataset(
    hd_data=data,
    name="Gaussian Line",
    verbose=True,
    hd_metric="euclidean",
)
```

<details>
<summary>How is the the Anndata object structured?</summary>

The TRACE backend can load data structured in the [Anndata](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) format. It includes the following fields:

* `adata.X` high-dimensional data
* [optional] `adata.obs`: dataframe with metadata e.g. cluster labels
* `adata.obsm` low-dimensional embeddings, one entry for each embedding, e.g. `adata.obsm["t-SNE (exag. 5)"]` for a t-SNE embedding. 
* `adata.uns` unstructured data:

    * `adata.uns["methods"]`: a dictionary that structures all available embeddings into groups (exactly one level with keys and a list as values such as in the example). This defines which embeddings can be selected in the interface. For example one could group according to DR methods and and list all corresponding **two-dimensional** embedding keys in adata.obsm:
        ```json
        {
            "t-SNE": ["t-SNE (exag. 5)", "t-SNE (exag. 1)"],
            "UMAP": ["UMAP 20", "UMAP 100"]
        }
        ```
    * [optional] `adata.uns["neighbors"]`: an _nxk_ array of the k-nearest high-dimensional neighbors of each point
    * [optional] `adata.uns["t-SNE (exag. 5)"]`: dictionary with additional data for each embedding, such as **quality** scores or **parameters** used to obtain the embedding. For example:
        ```json
        {
            "quality": {"qnx@50": [...], "qnx@200": [...]},
            "parameters": {"perplexity": 100, "exaggeration": 5, "epochs": 750}
        }
        ```
    * [optional] 🌈 You can add custom colors for metadata features by adding a list of HEX values to `trace_data.adata.uns["featureName_colors"]`. For categorical features, the number of colors should match the number of categories. The colors for continuous features will be mapped to the [min, max] range of the feature values.

</details>

### 1. Adding 2-dimensional embeddings
After preprocessing your data and computing a range of 2-dimensional embeddings using your favorite DR method, add the data and the embeddings to the data object:

```python
# Repeat for each embedding
trace_data.add_embedding(
    name= "tSNE (perplexity 30)",
    embedding = tsne_emb,
    category="tSNE",
)
```


### 2. Computing High-Dimensional Neighbors and Quality Measures

To provide snappy interactions in TRACE, the HD neighbors and a range of quality measures need to be precomputed. We use [ANNOY](https://github.com/spotify/annoy) to obtain the approximate neighbors and provide implementations of the following **quality measures** to be visualized via point colors in TRACE:

* **neighborhood preservation** measures the fraction of k high-dimensional neighbors that are preserved in the low-dimensional embedding. 
* **landmark distance correlation**: Sampling landmark points using a random or kmeanss++ (supports only Euclidean distance) from the high-dimensional data. We then compute the pairwise distances between all landmarks in high-dimensional space and each embedding and the rank correlation of their distance vectors. Points that are not landmark points are colored according to their nearest landmark point in the embedding. 
* **random triplet accuracy** quantifies the ratio of random triplets (i,j,k), where relative order of j and k with respect to i in the high-dimensional space is preserved in the embedding. 
* **point stability** measures how much the distances between each point and a random sample of other points vary across all embeddings. If a point has a very different global or local position in the embeddings, the stability will be low. 

To compute all available quality measures:
```python
trace_data.compute_quality(filename="./gauss_line.h5ad", hd_metric="euclidean")
trace_data.print_quality()
```

<details>
<summary>How can I chose the parameters of the quality measures?</summary>
Instead of calling the ```compute_quality``` function, you can also call each function separately.

```python
trace_data.precompute_HD_neighbors(maxK=200)
trace_data.compute_neighborhood_preservation(
    neighborhood_sizes=[200, 100, 50]
)
trace_data.compute_global_distance_correlation(
    max_landmarks=1000, LD_landmark_neighbors=True,
    hd_metric="euclidean", sampling_method="random",
)
trace_data.compute_random_triplet_accuracy(
    num_triplets=10
)
trace_data.compute_point_stability(num_samples=50)

# alighn the embeddings such that point movement is minimized
trace_data.align_embeddings(reference_embedding="PCA")
trace_data.save_adata(filename="./gauss_line.h5ad")

```
</details>

### 3. Add Dataset Configuration

To include a dataset in the dashboard you need to extend the configuration in [data_configs.yaml](./backend/data_configs.yaml). For the Gaussian Line dataset the configuration would be:
```json
"GaussLine": {
    "filepath": "../data/gauss_line/gauss_line.h5ad",
    "name": "GaussLine",
    "description": "Gaussian clusters shifted along a line from Böhm et al. (2022)",
}
```

## Example Datasets

### Gaussian Line 🟢 🟠 🟣
A small example dataset that is included in the repository. 

### Mammoth 🦣
This dataset from Wang et al. can be downloaded from their [PaCMAP](https://github.com/YingfanWang/PaCMAP/blob/master/data/mammoth_3d_50k.json) repository. It then needs to be processed using the `mammoth.ipynb` notebook. 

### Single-Cell Mouse Data 🐁
The processed dataset of gene expressions from [Guilliams et al.](https://pubmed.ncbi.nlm.nih.gov/35021063/) is not available online, please reach out if you are interested. A raw version is available under [GSE192742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192742).


***

<a name="trace">[1]</a> TRACE stands for Two-dimensional representation Analysis and Comparison Engine<br />
<a name="regl_citation">[2]</a> Lekschas, Fritz. "Regl-Scatterplot: A Scalable Interactive JavaScript-based Scatter Plot Library." Journal of Open Source Software (2023)

[⬆️ Back to top](#pattern-or-artefact-interactively-exploring-embedding-quality-with-trace)
