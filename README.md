# Interactively Exploring Embedding Quality with TRACE

TRACE<sup>[1](#trace)</sup> supports you in analyzing **global and local quality 🕵🏽‍♀️** of two-dimensional embeddings, based on [Regl-scatterplot](https://github.com/flekschas/regl-scatterplot)<sup>[2](#regl_citation)</sup> .


## Installation

### Option 1: Using Docker 🐋

```bash
docker-compose build
docker-compose up
```
This will mount the /frontend, /backend, and /data directories into the repective containers.
The /frontend/.next folder will be recreated every time the frontend is started. One might have to delete this folder before running the app without docker (as it is owned by 'root'). 

Open [http://localhost:3000](http://localhost:3000) with your browser to see the result.

### Option 2: Without Docker

#### Required packages
**Backend**: Install the required python packages for the backend, tested with Python 3.11 from `backend/pip_requirements.txt` or `backend/conda_requirements.yaml`. 

**Frontend**: Install the packages in `frontend/package.json` using e.g. `npm install`.


#### Usage

💡 Did you use Docker before to run TRACE?<br />
Make sure your user has write access to /frontend/.next or delete this folder. 

First, start the backend within the right python evironment:
```bash
conda activate backend_env/
python main.py
```

Then start the frontend development server:
```bash
npm run dev
# or
yarn dev
# or
pnpm dev
```
Open [http://localhost:3000](http://localhost:3000) with your browser to see the result.

## Data

The datasets have to be prepared as an [Anndata](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) objects including the following fields:

* `adata.X` high-dimensional data
* [optional] `adata.obs`: dataframe with metadata e.g. cluster labels
* `adata.obsm` low-dimensional embeddings, one entry for each embedding, e.g. `adata.obsm["t-SNE (exag. 5)"]` for a t-SNE embedding. 
* `adata.uns` unstructured data:

    * `adata.uns["methods"]`: a dictionary with DR methods as keys and a list of all **two-dimensional** embedding keys as value. For example:
        ```json
        {
            "t-SNE": ["t-SNE (exag. 5)", "t-SNE (exag. 1)"]
            "UMAP": ["UMAP 20", "UMAP 100"]
        }
        ```
    * [optional] `adata.uns["neighbors"]`: an _nxk_ array of the k-nearest high-dimensional neighbors of each point
    * [optional] `adata.uns["t-SNE (exag. 5)"]`: dictionary with additional data for each embedding, such as **quality** scores or **parameters** used to obtain the embedding. For example:
        ```json
        {
            "quality": {"qnx@50": [...], "qnx@200": [...]}
            "parameters": {"perplexity": 100, "exaggeration": 5, "epochs": 750}
        }
        ```

To include a dataset in the dashboard you need to extend the configuration in [data_configs.yaml](./backend/data_configs.yaml). For the Gaussian Line dataset the configuration would be:
```json
"GaussLine": {
    "filepath": "../data/gauss_line/gauss_line.h5ad",
    "name": "GaussLine",
    "hd_metric": "euclidean",
    "description": "Gaussian clusters shifted along a line from Böhm et al. (2022)",
    "hd_data_key": "X",
}
```

## Quality Measures and High-Dimensional Neighbors

When the backend is started it will compute high-dimensional neighbors and quality measures when they are not yet part of the anndata object. Specifically, it will compute
* the 200 nearest high-dimensional neighbors using the specified *hd_metric* in the configuration. 
* **neighborhood preservation** scores for `k = 50` and `k = 200` for each embedding. The values are computed as the size of the intersection between the sets of high-dimensional and low-dimensional neighbors, normalized by `k`.
* **landmark distance correlation**: Sampling `n = 1000` landmark points using kmeanss++ (supports only Euclidean distance) from the high-dimensional data. We then compute the pairwise distances between all landmarks in high-dimensional space and each embedding and the rank correlation of their distance vectors. Points that are not landmark points are colored according to their nearest landmark point in the embedding. 

The anndata object including the new quality scores will be saved as a new file and can be used by updating the filepath in the config.

### Gaussian Line 🟢 🟠 🟣
A small example dataset that is included in the repository. UMAP and t-SNE are already precomputed, but the HD neighbors and quality measures are computed when the dashboard is loaded.

### Mammoth 🦣
This dataset is from Wang et al. and can be downloaded from their [PaCMAP](https://github.com/YingfanWang/PaCMAP/blob/master/data/mammoth_3d_50k.json) repository. It then needs to be processed using the `mammoth.ipynb` notebook. 

### Single-Cell Mouse Data 🐁
The processed dataset of gene expressions from [Guilliams et al.](https://pubmed.ncbi.nlm.nih.gov/35021063/) is not available online, please reach out if you are interested. A raw version is available under [GSE192742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192742).


***

<a name="trace">[1]</a> TRACE stands for Two-dimensional representation Analysis and Comparison Engine<br />
<a name="regl_citation">[2]</a> Lekschas, Fritz. "Regl-Scatterplot: A Scalable Interactive JavaScript-based Scatter Plot Library." Journal of Open Source Software (2023)