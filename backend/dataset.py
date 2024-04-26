import numpy as np
import pandas as pd
import os
import neighbors
import time
import glasbey
import anndata as ad
import centroid_correlation
import utils
import json
from sklearn.neighbors import NearestNeighbors
from numba import set_num_threads, get_num_threads
from alignment import umeyama_alignment
import evaluation.triplet_accuracy as triplet_accuracy
from evaluation.stability import stability_across_embeddings
from collections.abc import Mapping


class Dataset:
    """
    The Dataset class is used to store and manage the data for a single dataset.

    The quality scores are stored in adata.uns['embedding_name']['quality'] as a dictionary.
    The dictionary adata.uns['methods'] contains the available embeddings
    as keys and the number of different embeddings from that method as values.
    """

    def __init__(
        self,
        name: str,
        hd_data: np.ndarray = None,
        adata: ad.AnnData = None,
        filepath: str = None,
        hd_data_key: str = "X",
        hd_metric: str = None,
        description: str = None,
        verbose: bool = False,
        n_jobs: int = 8,
        **kwargs,
    ):
        """
        Initialize a new dataset.

        Args:
            name (str): The name of the dataset.
            hd_data (np.ndarray, optional): The high-dimensional data array. Defaults to None.
            adata (ad.AnnData, optional): The AnnData object. Defaults to None.
            filepath (str, optional): The path to the h5ad file. Defaults to None.
            hd_data_key (str, optional): The key for the high-dimensional data array in adata.obsm or "X".
            hd_metric (str, optional): The metric to use for high-dimensional space. Defaults to None.
            description (str, optional): A description of the dataset. Defaults to None.
        """
        self.name = name
        self.hd_metric = hd_metric
        self.description = description
        self.hd_data_key = hd_data_key
        self.filepath = filepath
        self.verbose = verbose
        self.n_jobs = n_jobs
        set_num_threads(n_jobs)

        self.hd_annoy_filepath = dict()

        if hd_data is not None:
            self.adata = ad.AnnData(X=np.asarray(hd_data))
        elif adata is not None:
            self.__load_anndata(adata, hd_data_key)
        elif filepath is not None:
            self.__load_anndata(filepath, hd_data_key)
            self.filepath = filepath
        else:
            raise ValueError("No data provided.")

        self.__check_andata()

    def __load_anndata(self, adata: ad.AnnData | str, hd_data_key: str):
        """
        Load an AnnData object into the dataset.

        Args:
            adata (ad.AnnData | str): The AnnData object or path to the h5ad file.
            hd_data_key (str): The key for the high-dimensional data array in adata.obsm or "X".
        """
        self.hd_data_key = hd_data_key

        if isinstance(adata, str):
            if os.path.isfile(adata):
                self.filepath = adata
                self.adata = ad.read_h5ad(adata, backed="r")
            else:
                raise FileNotFoundError(f"File {adata} not found.")
        elif isinstance(adata, ad.AnnData):
            self.adata = adata
            self.hd_data_key = hd_data_key

        self.__check_andata()
        if self.verbose:
            print(
                f"\nLoaded {self.name} data.\n"
                + f"Embeddings {self.get_embedding_names()}.\n"
                + f"Quality features {self.get_quality_features()}."
            )

    def __check_andata(self):
        """
        Check if the AnnData object is valid.
        - embeddings are normalized between [-1, 1]
        - 'hd_neighbors' is a dictionary in adata.uns
        - every embedding has a 'quality' dictionary in adata.uns[embedding_name]
        """
        # check HD data key
        if self.hd_data_key not in self.adata.obsm_keys() and self.hd_data_key != "X":
            raise ValueError(f"HD data key {self.hd_data_key} not found in adata.obsm.")

        # keep list of embeddings in adata.uns['methods]
        if "methods" not in self.adata.uns_keys():
            self.adata.uns["methods"] = dict()

            # add existing embeddings to adata.uns['methods']
            for name in self.adata.obsm_keys():
                emb_dim = self.adata.obsm[name].shape[1]
                if emb_dim == 2:
                    self.adata.uns["methods"]["other"] = self.adata.uns["methods"].get(
                        "other", []
                    ) + [name]
                elif "pca" in name or "PCA" in name:
                    # add first two dimensions as extra embedding
                    self.adata.obsm[name + "_2D"] = self.adata.obsm[name][:, :2]
                    self.adata.uns["methods"]["other"] = self.adata.uns["methods"].get(
                        "other", []
                    ) + [name + "_2D"]

        # check that all embeddings are in obsm
        for name in self.get_embedding_names():
            if name not in self.adata.obsm_keys():
                raise ValueError(f"Embedding {name} not found in adata.obsm.")

        # normalization
        for name in self.get_embedding_names():
            self.adata.obsm[name] = utils.normalizeEmbedding(self.adata.obsm[name])

        if "hd_neighbors" not in self.adata.uns_keys():
            self.adata.uns["hd_neighbors"] = dict()

        # ensure that we have a dictionary in adata.uns for each embedding
        for name in self.get_embedding_names():
            if name not in self.adata.uns_keys():
                self.adata.uns[name] = {"quality": {}}
            elif "quality" not in self.adata.uns[name].keys():
                self.adata.uns[name]["quality"] = {}

        if self.hd_metric is None and len(self.adata.uns["hd_neighbors"].keys()) > 0:
            self.hd_metric = list(self.adata.uns["hd_neighbors"].keys())[0]

    def add_embedding(
        self,
        embedding: np.ndarray,
        name: str,
        category: str = None,
        meta_info: dict = None,
    ):
        """
        Add a new embedding to the dataset.

        Parameters:
            embedding (np.ndarray): The embedding to be added.
            name (str): The name of the embedding.
            category (str, optional): The category of the embedding (e.g. DR method). Defaults to None.
            meta_info (dict, optional): Additional metadata for the embedding. Defaults to None.
                Could be precomputed quality scores or parameter information. Example:
                {"parameters": {"perplexity": 30, "learning_rate": 200, "n_iter": 1000}}
        """
        if name in self.adata.uns.keys():
            raise ValueError(f"Embedding {name} already exists.")

        if category is None:
            self.adata.uns["methods"]["other"] = self.adata.uns["methods"].get(
                "other", []
            ) + [name]
        elif category in self.adata.uns["methods"]:
            self.adata.uns["methods"][category].append(name)
        else:
            self.adata.uns["methods"][category] = [name]

        self.adata.obsm[name] = utils.normalizeEmbedding(embedding)

        if meta_info is not None:
            # check if quality scores are point-wise
            if "quality" in meta_info.keys():
                for key, value in meta_info["quality"].items():
                    if len(value) != self.adata.n_obs:
                        raise ValueError(
                            f"Quality score {key} should have a value for each point."
                        )
            else:
                meta_info["quality"] = {}
            self.adata.uns[name] = meta_info
        else:
            self.adata.uns[name] = {"quality": {}}

    def add_metadata(self, metadata: pd.DataFrame | Mapping[str, np.ndarray]):
        """
        Add metadata to the dataset.

        Args:
            metadata (pd.DataFrame | Mapping[str, np.ndarray]): The metadata to be added.
                If a pd.DataFrame is provided, it should have a value for each point in the dataset.
                If a Mapping[str, np.ndarray] is provided, each key-value pair represents a metadata column,
                where the key is the column name and the value is an array of values for each point.

        Raises:
            ValueError: If the length of the metadata does not match the number of points in the dataset.
        """
        if isinstance(metadata, pd.DataFrame):
            if len(metadata) != self.adata.n_obs:
                raise ValueError(f"Metadata should have a value for each point.")
            self.adata.obs = pd.concat([self.adata.obs, metadata], axis=1)
        else:
            for key, value in metadata.items():
                if len(value) != self.adata.n_obs:
                    raise ValueError(
                        f"Metadata {key} should have a value for each point."
                    )
                self.adata.obs[key] = value

    def rename_quality_column(self, old_name: str, new_name: str):
        """
        Rename a quality column in the dataset.

        Parameters:
            old_name (str): The old name of the quality column.
            new_name (str): The new name of the quality column.
        """
        for name in self.get_embedding_names():
            if old_name in self.adata.uns[name]["quality"]:
                self.adata.uns[name]["quality"][new_name] = self.adata.uns[name][
                    "quality"
                ].pop(old_name)

    def get_quality_features(self):
        """
        List of quality features that can be used to color points.
        """
        features = set()
        for embedding in self.get_embedding_names():
            features = features.union(self.adata.uns[embedding]["quality"].keys())
        return list(features)

    def get_metadata_features(self):
        """
        List of column names that can be used to color points.
        """
        return self.adata.obs_keys()

    def get_HD_data(self):
        """Return the ndarray to use for HD neighbor computation.
        e.g. obsm["X_pca"]
        """
        if self.hd_data_key in self.adata.obsm_keys():
            return self.adata.obsm[self.hd_data_key]
        else:
            return np.asarray(self.adata.X)

    def save_adata(self, filename=None):
        """
        Save the AnnData object to a file.

        Parameters:
            filename (str): Filename to save the h5ad file.
                If None, create a new name based on the dataset name and current timestamp
        """

        if filename is None:
            if self.filepath is None:
                filename = self.name.replace(" ", "_") + ".h5ad"

                if os.path.isfile(filename):
                    filename = (
                        self.name.replace(" ", "_")
                        + f"_{time.strftime('%Y%m%d_%H%M%S')}.h5ad"
                    )
            else:
                filename = (
                    os.path.splitext(os.path.basename(self.filepath))[0]
                    + f"_{time.strftime('%Y%m%d_%H%M%S')}.h5ad"
                )
                filename = os.path.join(os.path.dirname(self.filepath), filename)

        if self.verbose:
            print(f"Saving dataset to {filename}")

        self.adata.write(
            filename,
            compression="gzip",
        )

    def __get_annoy_index(self, metric: str):
        """
        Builds an Annoy index for the dataset using the specified metric.

        Parameters:
            metric (str): The metric to be used for building the index.
        """
        if self.hd_annoy_filepath.get(metric, None) is None:
            self.hd_annoy_filepath[metric] = neighbors.build_annoy_index(
                self.get_HD_data(),
                metric=metric,
                filepath=os.path.join(
                    (
                        os.path.dirname(self.filepath)
                        if self.filepath is not None
                        else "./"
                    ),
                    self.name.replace(" ", "_") + "_" + metric + "_annoy.ann",
                ),
            )

        return self.hd_annoy_filepath[metric]

    def cleanup(self):
        """
        Clean up the dataset by removing any existing HD_annoy files.

        This method iterates over the HD_annoy filepaths and removes the corresponding files if they exist.
        """
        for fpath in self.hd_annoy_filepath.values():
            if fpath is not None and os.path.isfile(fpath):
                os.remove(fpath)

    def get_hd_metric_options(self):
        if "hd_neighbors" in self.adata.uns.keys():
            return list(self.adata.uns["hd_neighbors"].keys())
        else:
            return []

    def get_precomputed_neighbors_maxK(self):
        if "hd_neighbors" in self.adata.uns.keys():
            return min([v.shape[1] for v in self.adata.uns["hd_neighbors"].values()])
        else:
            return 0

    def get_embedding_options(self):
        """
        Returns a list with a dictionary for each method. This is the format required by the grouped-option select.

        Returns:
            list: Containing dictionaries for each method.
                  Example: {label: tSNE,
                            options: [{value: tSNE_perplexity_30, label: perplexity 30, group: tSNE},
                                      {value: tSNE_perplexity_50, label: perplexity 50, group: tSNE}]
        """
        options = []
        for method, embeddings in self.adata.uns["methods"].items():
            options.append(
                {
                    "label": method,
                    "options": [
                        {"value": name, "label": name, "group": method}
                        for name in embeddings
                    ],
                }
            )
        return options

    def get_pointcolor_options(self):
        options = []

        # Metadata
        options.append(
            {
                "label": "Metadata",
                "options": [
                    {"value": name, "label": name}
                    for name in self.get_metadata_features() + ["HD distances", "none"]
                ],
            }
        )

        # Quality
        options.append(
            {
                "label": "Quality",
                "options": [
                    {"value": name, "label": name}
                    for name in self.get_quality_features()
                ],
            }
        )

        # Data Features
        options.append(
            {
                "label": "Features",
                "options": [
                    {"value": name, "label": name}
                    for name in self.adata.var.index.tolist()
                ],
            }
        )
        return options

    def get_embedding_names(self):
        """
        Get the names of all embeddings in the dataset.

        Returns:
            list: A list of embedding names to be used to index into adata.obsm.
        """
        names = []
        for embeddingNames in self.adata.uns["methods"].values():
            names.extend(embeddingNames)
        return names

    def save_user_annotation(self, points: list[int], selection_name: str):
        """
        Save the user annotations for the dataset.

        Args:
            points (list[int]): List with point indices.
            selection_name (str): Description of the selection.
        """
        fname = os.path.join(
            (os.path.dirname(self.filepath) if self.filepath is not None else "./"),
            self.name.replace(" ", "_")
            + f"_user_annotations_{time.strftime('%Y%m%d')}.json",
        )

        if "user_annotations" not in self.adata.uns:
            self.adata.uns["user_annotations"] = {}

        if selection_name in self.adata.uns["user_annotations"]:
            # add current time as stinrg to make it unique
            selection_name = selection_name.replace(" ", "_")
            selection_name += f"_{time.strftime('%H%M%S')}"

        self.adata.uns["user_annotations"][selection_name] = points
        json.dump(self.adata.uns["user_annotations"], open(fname, "w"))
        # need to change the ownership of the file to the user (when running from docker)
        os.chown(fname,
                 os.stat(os.path.dirname(self.filepath)).st_uid, 
                 os.stat(os.path.dirname(self.filepath)).st_gid)
        print("Wrote user annotations to ", fname)

    def get_category_colors(self, fname, categories: list):
        """
        Get the colors associated with each category in the specified feature.

        Args:
            fname (str): The name of the feature.

        Returns:
            dict: A dictionary with "ticks" and "colors" keys.
        """
        if fname + "_colors" in self.adata.uns_keys():
            colors = self.adata.uns[fname + "_colors"][: len(categories)]
        else:
            colors = glasbey.create_palette(
                palette_size=len(categories), grid_size=32, colorblind_safe=True
            )
        return {"ticks": categories, "colors": list(colors)}

    def align_embeddings(self, reference_embedding: str):
        if reference_embedding not in self.get_embedding_names():
            raise ValueError(f"Reference embedding {reference_embedding} not found.")

        for name in self.get_embedding_names():
            if name != reference_embedding:
                self.adata.obsm[name] = umeyama_alignment(
                    self.adata.obsm[name], self.adata.obsm[reference_embedding]
                )
                # make sure the embedding is still normalized
                self.adata.obsm[name] = utils.normalizeEmbedding(self.adata.obsm[name])

    #############################
    ## Quality Score Functions ##
    #############################

    def compute_quality(self, filename: str = None):
        """
        Compute all quality measures for the dataset.
        """
        # this will take a while, the user should be updated
        self.verbose = True
        print(f"Quality measures for {self.name}...")
        print(f"NUMBA uses {get_num_threads()} threads")
        start_time = time.time()

        if self.adata.n_obs > 50:
            num_samples = int(
                10 + 20 * (np.log10(self.adata.n_obs) - 3)
                if self.adata.n_obs > 1000
                else 10
            )
            neighborhood_sizes = [5 * num_samples, num_samples]
        else:
            num_samples = int(self.adata.n_obs / 2)
            neighborhood_sizes = [num_samples]

        self.precompute_HD_neighbors(maxK=max(neighborhood_sizes))
        self.compute_neighborhood_preservation(neighborhood_sizes=neighborhood_sizes)
        self.compute_global_distance_correlation()

        self.compute_random_triplet_accuracy(num_triplets=num_samples)
        self.compute_point_stability(num_samples=num_samples)

        self.align_embeddings(reference_embedding=self.get_embedding_names()[0])
        print(f"Done in {time.time()-start_time:.2f}s.")
        self.save_adata(filename=filename)

    def print_quality(self):
        """
        Print the quality of the dataset for each method.

        This method iterates over the methods in the dataset and prints the quality
        information for each method. It calculates the mean value and the percentage
        of zeros for each embedding in the dataset.

        Returns:
            None
        """
        avg_quality_df = pd.DataFrame(
            columns=self.get_quality_features(), index=self.get_embedding_names()
        )

        for name in self.get_embedding_names():
            for quality_key in self.adata.uns[name]["quality"].keys():
                avg_quality_df.loc[name, quality_key] = self.adata.uns[name]["quality"][
                    quality_key
                ].mean()
        print(avg_quality_df)

    def precompute_HD_neighbors(self, maxK: int, metric: str = None, exact=False):
        """
        Precomputes the high-dimensional (HD) neighbors for each data point.

        Args:
            maxK (int): The maximum number of neighbors to compute.
            metric (str): The distance metric to use for computing neighbors.
            exact (bool, optional): Whether to compute exact neighbors or use an approximate algorithm. Defaults to False.
        """

        if metric is None:
            metric = self.hd_metric
        if (
            metric in self.adata.uns["hd_neighbors"]
            and self.adata.uns["hd_neighbors"][metric].shape[1] < maxK
        ) or (metric not in self.adata.uns["hd_neighbors"]):
            if self.verbose:
                print(f"Computing all HD neighbors until k = {maxK}")
            start_time = time.time()
            self.adata.uns["hd_neighbors"][metric] = neighbors.get_nearest_neighbors(
                data=self.get_HD_data(),
                indices=range(self.adata.n_obs),
                k=maxK,
                metric=metric,
                filepath=self.__get_annoy_index(metric),
                exact=exact,
                n_jobs=self.n_jobs,
            )
            if self.verbose:
                print(f"Done ({time.time()-start_time:.2f}s).")

        if os.path.isfile(self.__get_annoy_index(metric)):
            os.remove(self.__get_annoy_index(metric))

    def get_HD_neighbors(self, k: int, metric: str, indices: np.ndarray):
        """
        Get the high-dimensional neighbors for the given indices.

        Parameters:
            k (int): The number of neighbors to retrieve.
            metric (str): The distance metric to use for neighbor search.
            indices (np.ndarray): The indices of the data points for which to find neighbors.

        Returns:
            np.ndarray: An array of shape (len(indices), k) with the high-dimensional neighbors
            for the given indices.
        """
        if (
            metric in self.adata.uns["hd_neighbors"]
            and self.adata.uns["hd_neighbors"][metric].shape[1] >= k
        ):
            return self.adata.uns["hd_neighbors"][metric][indices, :k]
        elif indices.shape[0] > self.adata.n_obs / 2:
            # precompute all neighbors
            self.precompute_HD_neighbors(maxK=k, metric=metric)
            return self.adata.uns["hd_neighbors"][metric][indices, :k]
        else:
            # compute neighbors on the fly
            return neighbors.get_nearest_neighbors(
                data=self.get_HD_data(),
                indices=indices,
                k=k,
                metric=metric,
                filepath=self.__get_annoy_index(metric),
                n_jobs=self.n_jobs,
            )

    def compute_neighborhood_preservation(
        self,
        neighborhood_sizes: list[int],
        hd_metric: str = None,
        ld_metric: str = "euclidean",
    ):
        """
        Computes the quality scores for embeddings based on neighborhood preservation.
        It stores the quality scores in the adata.uns dictionary of every embedding.

        Args:
            k (int): The number of nearest neighbors to consider.
            hd_metric (str): The metric to use for high-dimensional space.
            ld_metric (str, optional): The metric to use for low-dimensional space. Defaults to "euclidean".
        """
        if hd_metric is None:
            hd_metric = self.hd_metric
        if hd_metric not in self.adata.uns["hd_neighbors"] or self.adata.uns[
            "hd_neighbors"
        ][hd_metric].shape[1] < max(neighborhood_sizes):
            if self.verbose:
                print("Computing HD neighbors before calculating quality")
            self.precompute_HD_neighbors(
                maxK=max(neighborhood_sizes), hd_metric=hd_metric
            )
        if self.verbose:
            print("Computing neighborhood preservation scores...", end="")
        start_time = time.time()
        for name in self.get_embedding_names():
            if any(
                [
                    f"neighborhood preservation k={k}"
                    not in self.adata.uns[name]["quality"]
                    for k in neighborhood_sizes
                ]
            ):
                if self.verbose:
                    print(f"{name}, ", end="")
                preservation = neighbors.neighborhood_preservation_multi(
                    hd_neighbors=self.adata.uns["hd_neighbors"][hd_metric],
                    embedding=self.adata.obsm[name],
                    neighborhood_sizes=neighborhood_sizes,
                    ld_metric=ld_metric,
                )
                for i, k in enumerate(neighborhood_sizes):
                    self.adata.uns[name]["quality"][
                        f"neighborhood preservation k={k}"
                    ] = preservation[:, i]
        if self.verbose:
            print(f"Done ({time.time()-start_time:.2f}s).")

    def compute_global_distance_correlation(
        self,
        hd_metric: str = None,
        max_landmarks: int = 10000,
        LD_landmark_neighbors: bool = True,
        sampling_method: str = "auto",
    ):
        """
        Computes the quality scores for embeddings based on global distance correlation.
        It selects a subset of landmarks and computes the distance correlation between the
        low-dimensional and high-dimensional distances. Each point will be assigned the
        correlation score of the closest landmark.

        Args:
            hd_metric (str, optional): Metric to compute HD distances. Defaults to None.
            max_landmarks (int, optional): The maxium number of landmarks. Defaults to 10000.
            LD_landmark_neighbors (bool, optional): If the closest landmark is
                defined in LD space (True) or HD space (False). Defaults to True.
            sampling_method (str, optional): How to sample the landmarks. Default is "auto" which uses "kmeans++"
                for data up to 10000 points and "random" otherwise. Options ["kmeans++", "random", "auto"]
        """
        start_time = time.time()
        if hd_metric is None:
            hd_metric = self.hd_metric

        sampling_choices = ["kmeans++", "random", "auto"]
        if sampling_method not in sampling_choices:
            raise ValueError(
                f"Sampling method {sampling_method} not found. \
                             Options are {sampling_choices}."
            )
        if sampling_method == "auto":
            if self.get_HD_data().shape[1] <= 10000:
                sampling_method = "kmeans++"
            else:
                sampling_method = "random"

        num_samples = 100 if self.adata.n_obs < 1000 else int(self.adata.n_obs**0.7)
        num_samples = min(num_samples, max_landmarks)

        # get embeddings for which the quality score is not yet computed
        embedding_names = [
            name
            for name in self.get_embedding_names()
            if "landmark distance correlation" not in self.adata.uns[name]["quality"]
        ]

        if len(embedding_names) == 0:
            if self.verbose:
                print(
                    "Skipping landmark distance correlation. It is already available for all embeddings."
                )
            return

        if self.verbose:
            print("Computing landmark distance correlation...", end="", flush=True)

        if "landmark_indices" not in self.adata.uns_keys():
            if sampling_method == "kmeans++" and self.get_HD_data().shape[1] > 10000:
                print(
                    f"Sampling {num_samples} landmarks using kmeans++ for data of shape "
                    + f"{self.get_HD_data().shape} will take a while. Consider using sampling_method='random' instead."
                )
            if sampling_method == "kmeans++":
                landmark_indices = centroid_correlation.sample_landmarks(
                    data=self.get_HD_data(),
                    num_samples=num_samples,
                )
            else:
                landmark_indices = np.random.choice(
                    self.adata.n_obs, size=num_samples, replace=False
                )
            self.adata.uns["landmark_indices"] = landmark_indices

        hd_landmark_distances = centroid_correlation.compute_pairwise_distance(
            X=self.get_HD_data()[self.adata.uns["landmark_indices"]],
            Y=None,
            metric=hd_metric,
        )

        for name in embedding_names:
            if self.verbose:
                print(f"{name}, ", end="", flush=True)
            self.adata.uns[name]["quality"]["landmark distance correlation"] = (
                centroid_correlation.compute_landmark_correlation(
                    ld_data=self.adata.obsm[name],
                    hd_data=self.get_HD_data(),
                    landmark_indices=self.adata.uns["landmark_indices"],
                    hd_landmark_distances=hd_landmark_distances,
                    LD_landmark_neighbors=LD_landmark_neighbors,
                    hd_metric=hd_metric,
                    ld_metric="euclidean",
                )
            )
        if self.verbose:
            print(f"Done ({time.time()-start_time:.2f}s).")

    def get_HD_landmark_distances(self, selectedPoint: int):
        """
        Get the distances of the closest landmark of the selected point to all other landmarks.

        Parameters:
            selectedPoint (int): The index of the selected point.

        Returns:
            np.ndarray: An array of shape (n_points,) with the prolongated distance.
        """

        if ["closestLandmarks"] not in self.adata.uns_keys():
            # nearest neighbor index for HD data
            hd_nbr_index = NearestNeighbors(n_neighbors=2, algorithm="auto").fit(
                self.get_HD_data()[self.adata.uns["landmark_indices"], :]
            )
            _, neighbors = hd_nbr_index.kneighbors(self.get_HD_data(), n_neighbors=1)
            neighbors = neighbors.flatten()
            # translate neighbor indices to landmark indices
            # neighbors = self.adata.uns["landmark_indices"][neighbors]
            self.adata.uns["closestLandmarks"] = neighbors

        closestLandmark = self.adata.uns["landmark_indices"][
            self.adata.uns["closestLandmarks"][selectedPoint]
        ]

        distance_to_landmarks = centroid_correlation.compute_pairwise_distance(
            X=self.get_HD_data()[self.adata.uns["landmark_indices"]],
            Y=self.get_HD_data()[closestLandmark].reshape(1, -1),
            metric=self.hd_metric,
        ).flatten()

        # prolongate distances to all points
        distances = np.take(
            distance_to_landmarks, self.adata.uns["closestLandmarks"], axis=0
        )
        return distances

    def compute_point_stability(
        self, embedding_keys: list[str] = None, num_samples: int = 50, alpha: float = 20
    ):
        """
        Compute the point stability for the given embeddings.

        Args:
            embedding_keys (list[str], optional): Subset of embeddings to
                compute the point stability for. Defaults to None.
            num_samples (int, optional): Number of samples to use for the distances
                variance computation. Defaults to 50.
            alpha (float, optional): Exponent for the stability calculation.
                Values between 5 and 50 are typically sufficient. Defaults to 20.
        """
        if self.verbose:
            print(
                f"Computing stability scores with num_samples={num_samples} and alpha={alpha}...",
                flush=True,
            )

        if embedding_keys is None:
            embedding_keys = self.get_embedding_names()
        else:
            # check if all embeddings are available
            for key in embedding_keys:
                if key not in self.get_embedding_names():
                    raise ValueError(f"Embedding {key} not found.")

        embeddings = [self.adata.obsm[key] for key in embedding_keys]

        start_time = time.time()
        stability = stability_across_embeddings(
            embeddings=embeddings, num_samples=num_samples, alpha=alpha
        )
        if self.verbose:
            print(f"Done ({time.time()-start_time:.2f}s).")

        # add stability to metadata
        self.add_metadata({"point stability": stability})

    def delete_quality_scores(self, quality_scores: list[str]):
        """
        Delete quality scores from the dataset.

        Args:
            quality_scores (list[str]): The quality scores to be deleted.
        """
        for name in self.get_embedding_names():
            for quality_score in quality_scores:
                if quality_score in self.adata.uns[name]["quality"]:
                    del self.adata.uns[name]["quality"][quality_score]

    def compute_random_triplet_accuracy(
        self, num_triplets=10, hd_metric=None, same_triplets=True
    ):
        if hd_metric is None:
            hd_metric = self.hd_metric

        if self.verbose:
            print(
                f"Computing random triplet accuracy with {num_triplets} per point...",
                flush=True,
            )
        start_time = time.time()

        if same_triplets:
            if self.verbose:
                print(f"\tGetting HD triplets and their distances...", flush=True)
            labels, triplets = triplet_accuracy.compute_hd_triplets(
                X=self.get_HD_data(),
                num_triplets=num_triplets,
                hd_metric=hd_metric,
            )
            for name in self.get_embedding_names():
                if "random triplet accuracy" not in self.adata.uns[name]["quality"]:
                    if self.verbose:
                        print(f"{name}, ", end="", flush=True)
                    self.adata.uns[name]["quality"]["random triplet accuracy"] = (
                        triplet_accuracy.get_triplet_accuracy(
                            Y=self.adata.obsm[name],
                            hd_labels=labels,
                            triplets=triplets,
                        )
                    )
        else:
            for name in self.get_embedding_names():
                if "random triplet accuracy" not in self.adata.uns[name]["quality"]:
                    if self.verbose:
                        print(f"{name}, ", end="", flush=True)
                    self.adata.uns[name]["quality"]["random triplet accuracy"] = (
                        triplet_accuracy.random_triplet_eval(
                            X=self.get_HD_data(),
                            Y=self.adata.obsm[name],
                            num_triplets=num_triplets,
                            hd_metric=hd_metric,
                        )
                    )
        if self.verbose:
            print(f"Done ({time.time()-start_time:.2f}s).")
