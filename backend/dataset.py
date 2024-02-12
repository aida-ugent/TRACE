import numpy as np
import pandas as pd
import os
import neighbors
import time
import glasbey
import anndata as ad
import centroid_correlation
import utils
from sklearn.neighbors import NearestNeighbors


class Dataset:
    """
    Embeddings are stored in adata.obsm[methodName_methodIndex]
    and should be normalized to [-1,1].

    The quality scores are stored in adata.uns[methodName_methodIndex][quality_score_name]

    The dictionary adata.uns['methods'] contains the available embeddings
    as keys and the number of different embeddings from that method as values.

    """

    def __init__(
        self,
        filepath: str,
        name: str,
        hd_metric: str,
        description: str,
        hd_data_key: str,
        metadata_features: list = None,
    ):
        self.filepath = filepath
        self.name = name
        self.hd_metric = hd_metric
        self.description = description
        self.hd_data_key = hd_data_key
        self.metadata_features = metadata_features
        self.hd_annoy_filepath = dict()

        if not os.path.isfile(self.filepath):
            raise FileNotFoundError(f"File {self.filepath} not found.")
        self.adata = ad.read_h5ad(self.filepath)
        self.check_andata()

        if self.adata.uns["hd_neighbors"].get(self.hd_metric, None) is None:
            print("Computing global HD neighbors")
            self.precompute_HD_neighbors(maxK=200, metric=self.hd_metric)

        self.rename_quality_column("quality qnx@50", "neighborhood preservation k=50")
        self.rename_quality_column("quality qnx@200", "neighborhood preservation k=200")
        self.rename_quality_column(
            "quality landmark corr", "landmark distance correlation"
        )

        # make sure the local quality is computed for all embeddings
        old_quality_keys = self.get_quality_features()
        self.compute_landmark_correlation(hd_metric=self.hd_metric, recompute=False)
        self.compute_quality(k=50, hd_metric=self.hd_metric)
        self.compute_quality(k=200, hd_metric=self.hd_metric)
        if self.get_quality_features() != old_quality_keys:
            print("Quality features changed. Saving dataset.")
            self.save_adata(overwrite=False)

        print(
            f"Loaded {self.name} data.\n"
            + f"Embeddings {self.get_embedding_names()}.\n"
            + f"Quality features {self.get_quality_features()}."
        )

    def check_andata(self):
        """
        Check if the AnnData object is valid.
        - embeddings are normalized between [-1, 1]
        - 'hd_neighbors' is a dictionary in adata.uns
        - every embedding has a 'quality' dictionary in adata.uns[embedding_name]
        """
        # check HD data key
        if self.hd_data_key not in self.adata.obsm_keys() and self.hd_data_key != "X":
            raise ValueError(f"HD data key {self.hd_data_key} not found in adata.obsm.")

        # check that all embeddings are in obsm
        for name in self.get_embedding_names():
            if name not in self.adata.obsm_keys():
                raise ValueError(f"Embedding {name} not found in adata.obsm.")

        # normalization
        for name in self.get_embedding_names():
            if np.max(np.abs(self.adata.obsm[name])) > 1.0:
                print(f"Normalizing embedding {name} to [-1,1].")
                self.adata.obsm[name] = utils.normalizeEmbedding(self.adata.obsm[name])

        if "hd_neighbors" not in self.adata.uns_keys():
            self.adata.uns["hd_neighbors"] = dict()

        # ensure that we have a dictionary in adata.uns for each embedding
        for name in self.get_embedding_names():
            if name not in self.adata.uns_keys():
                self.adata.uns[name] = {"quality": {}}
            elif "quality" not in self.adata.uns[name].keys():
                self.adata.uns[name]["quality"] = {}

    def compute_landmark_correlation(self, hd_metric: str, recompute=False):
        if (
            "landmark distance correlation" not in self.get_quality_features()
            or recompute
        ):
            print("Computing landmark correlation")
            self.sample_HD_landmarks(hd_metric)

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
        if self.metadata_features is None:
            return self.adata.obs_keys()
        else:
            return self.metadata_features

    def get_HD_data(self):
        """Return the ndarray to use for HD neighbor computation.
        e.g. obsm["X_pca"]
        """
        if self.hd_data_key in self.adata.obsm_keys():
            return self.adata.obsm[self.hd_data_key]
        else:
            return self.adata.X

    def save_adata(self, overwrite: bool = False):
        """
        Save the AnnData object to a file.

        Parameters:
            overwrite (bool): If True, overwrite the existing file. If False, create a new file with a timestamp in the filename.
        """
        if overwrite:
            self.adata.write(self.filepath, compression="gzip")
            print(f"Saved dataset to {self.filepath}")
        else:
            new_filename = (
                os.path.splitext(os.path.basename(self.filepath))[0]
                + f"_{time.strftime('%Y%m%d_%H%M%S')}.h5ad"
            )
            self.adata.write(
                os.path.join(os.path.dirname(self.filepath), new_filename),
                compression="gzip",
            )
            print(f"Saved dataset to {new_filename}")

    def get_annoy_index(self, metric: str):
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
                    os.path.dirname(self.filepath),
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

    def get_category_colors(self, fname, categories: list):
        """
        Get the colors associated with each category in the specified feature.

        Args:
            fname (str): The name of the feature.

        Returns:
            dict: A dictionary mapping each category to its corresponding color.
        """
        if fname + "_colors" in self.adata.uns_keys():
            colors = self.adata.uns[fname + "_colors"][: len(categories)]
        else:
            colors = glasbey.create_palette(
                palette_size=len(categories), grid_size=32, colorblind_safe=True
            )
        return dict(zip(categories, colors))

    def print_quality(self):
        """
        Print the quality of the dataset for each method.

        This method iterates over the methods in the dataset and prints the quality
        information for each method. It calculates the mean value and the percentage
        of zeros for each embedding in the dataset.

        Returns:
            None
        """
        for name in self.get_embedding_names():
            print(f"Quality for {name}:")
            for quality_key in self.adata.uns[name]["quality"].keys():
                print(
                    f" {quality_key}: {self.adata.uns[name]['quality'][quality_key].mean():.3f}, "
                    f"{100*np.sum(self.adata.uns[name]['quality'][quality_key] == 0.0)/self.adata.n_obs:.3f}% zero"
                )

    def precompute_HD_neighbors(self, maxK: int, metric: str, exact=False):
        """
        Precomputes the high-dimensional (HD) neighbors for each data point.

        Args:
            maxK (int): The maximum number of neighbors to compute.
            metric (str): The distance metric to use for computing neighbors.
            exact (bool, optional): Whether to compute exact neighbors or use an approximate algorithm. Defaults to False.
        """
        if (
            metric in self.adata.uns["hd_neighbors"]
            and self.adata.uns["hd_neighbors"][metric].shape[1] < maxK
        ) or (metric not in self.adata.uns["hd_neighbors"]):
            print(f"Computing all HD neighbors until k = {maxK}")
            start = time.time()
            self.adata.uns["hd_neighbors"][metric] = neighbors.get_nearest_neighbors(
                data=self.get_HD_data(),
                indices=range(self.adata.n_obs),
                k=maxK,
                metric=metric,
                filepath=self.get_annoy_index(metric),
                exact=exact,
            )
            print(f"Computing neighbors in {time.time() - start:.2f} seconds")

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
                filepath=self.get_annoy_index(metric),
            )

    def compute_quality(self, k: int, hd_metric: str, ld_metric: str = "euclidean"):
        """
        Computes the quality scores for embeddings based on neighborhood preservation.
        It stores the quality scores in the adata.uns dictionary of every embedding.

        Args:
            k (int): The number of nearest neighbors to consider.
            hd_metric (str): The metric to use for high-dimensional space.
            ld_metric (str, optional): The metric to use for low-dimensional space. Defaults to "euclidean".
        """
        if (
            hd_metric not in self.adata.uns["hd_neighbors"]
            or self.adata.uns["hd_neighbors"][hd_metric].shape[1] < k
        ):
            print("Computing HD neighbors before calculating quality")
            self.precompute_HD_neighbors(maxK=k, hd_metric=hd_metric)

        embedding_names = []
        embeddings = []
        quality_key = f"neighborhood preservation k={k}"
        for name in self.get_embedding_names():
            if quality_key not in self.adata.uns[name]["quality"]:
                embedding_names.append(name)
                embeddings.append(self.adata.obsm[name])

        if len(embeddings) > 0:
            print(f"Computing quality for {embedding_names}")
            quality_scores = neighbors.neighborhood_preservation(
                hd_data=self.get_HD_data(),
                ld_data_arr=embeddings,
                k=k,
                hd_metric=hd_metric,
                ld_metric=ld_metric,
                hd_annoy_filepath=self.get_annoy_index(hd_metric),
                hd_neighbors=self.adata.uns["hd_neighbors"][hd_metric],
            )

            for name, quality_result in zip(embedding_names, quality_scores):
                self.adata.uns[name]["quality"][quality_key] = quality_result

    def compute_centroid_distances(self, hd_metric: str, label_key: str):
        """
        Compute the correlation between centroid distances in HD and LD space.

        Parameters:
            hd_metric (str): The metric used to compute distances between HD data points.
            label_key (str): Key for cluster labels in adata.obs.
        """
        hd_centroid_distances = centroid_correlation.compute_centroid_distances(
            data=self.get_HD_data(),
            labels=self.adata.obs[label_key],
            metric=hd_metric,
        )

        # store distances in adata.obs
        self.adata.obs = pd.merge(
            self.adata.obs,
            centroid_correlation.get_hd_centroid_distance_df(
                labels=self.adata.obs[label_key],
                label_key=label_key,
                hd_centroid_distances=hd_centroid_distances,
            ),
            on=label_key,
            how="left",
        )

        # compute correlations to ld centroids
        embedding_names = self.get_embedding_names()
        for name in embedding_names:
            ld_centroid_distances = centroid_correlation.compute_centroid_distances(
                data=self.adata.obsm[name],
                labels=self.adata.obs[label_key],
                metric="euclidean",
            )
            corr = centroid_correlation.distance_correlation(
                hd_centroid_distances, ld_centroid_distances
            )
            centroid_correlations = pd.merge(
                pd.DataFrame({"cluster": self.adata.obs[label_key]}),
                pd.DataFrame(
                    {
                        "cluster": self.adata.obs[label_key].unique(),
                        "centroid_correlation": corr,
                    }
                ),
                on="cluster",
                how="left",
            )["centroid_correlation"].values
            self.adata.uns[name]["quality"][
                "quality centroid_corr"
            ] = centroid_correlations

    def sample_HD_landmarks(
        self,
        hd_metric: str,
        max_samples: int = 1000,
        LD_landmark_neighbors: bool = True,
    ):
        """
        Sample HD landmarks and compute distance correlations.

        Args:
            hd_metric (str): The metric to use for computing pairwise distances in HD space.
            max_samples (int, optional): The maximum number of landmarks to sample. Defaults to 1000.
            LD_landmark_neighbors (bool, optional): Flag indicating whether to find the nearest
                landmark neighbor in LD space or HD space. Defaults to True.
        """
        landmark_indices = centroid_correlation.sample_landmarks(
            data=self.get_HD_data(),
            max_samples=max_samples,
        )
        print(f"Sampled {len(landmark_indices)} landmarks")
        self.adata.uns["landmark_indices"] = landmark_indices

        hd_landmark_distances = centroid_correlation.compute_pairwise_distance(
            X=self.get_HD_data()[landmark_indices], Y=None, metric=hd_metric
        )
        embedding_names = self.get_embedding_names()

        for name in embedding_names:
            ld_landmark_distances = centroid_correlation.compute_pairwise_distance(
                X=self.adata.obsm[name][landmark_indices], Y=None, metric="euclidean"
            )

            corr = centroid_correlation.distance_correlation(
                hd_landmark_distances, ld_landmark_distances
            )

            if LD_landmark_neighbors:
                # find the nearest landmark neighbor in LD space for all points
                # alternative: find the nearest landmark neighbor in HD space for all points.
                landmark_nbr_index = NearestNeighbors(
                    n_neighbors=2, algorithm="auto", metric="euclidean"
                ).fit(self.adata.obsm[name][landmark_indices])
                _, neighbors = landmark_nbr_index.kneighbors(
                    self.adata.obsm[name], n_neighbors=1
                )
            else:
                # nearest neighbor index for HD Data
                landmark_nbr_index = NearestNeighbors(
                    n_neighbors=2,
                    algorithm="auto",
                    metric="cosine" if hd_metric == "angular" else hd_metric,
                ).fit(self.get_HD_data()[landmark_indices, :])
                _, neighbors = landmark_nbr_index.kneighbors(
                    self.get_HD_data(), n_neighbors=1
                )

            neighbors = neighbors.flatten()

            # place the remaining points at the location of their nearest neighbor
            prolongated_correlations = np.take(corr, neighbors, axis=0)
            print(
                f"Landmark correlation for embedding {name}: {prolongated_correlations.mean()}"
            )
            self.adata.uns[name]["quality"][
                "landmark distance correlation"
            ] = prolongated_correlations

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
