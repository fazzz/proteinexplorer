"""Ensemble clustering (spec section "Ensemble Clustering").

Clusters a set of conformations (either several structures already in the
project, or several MODEL records within one multi-model file such as an
NMR ensemble) by pairwise RMSD, same two-tier pattern as ChemExplorer's
and BioExplorer's clustering commands:

- greedy(): pure-Python, no extra dependency (CD-HIT-style incremental
  clustering -- each conformation joins the first existing cluster whose
  seed is within `threshold` RMSD, or starts a new cluster otherwise).
- hierarchical(): agglomerative clustering via scipy.cluster.hierarchy,
  requires the optional `cluster` extra (`pip install -e ".[cluster]"`).

Both report a medoid (the member with the lowest average RMSD to the rest
of its cluster) as the cluster representative.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from proteinexplorer import geometry as geom


class ClusterError(ValueError):
    pass


def pairwise_rmsd_matrix(atom_groups: list[list], fit: bool = True) -> np.ndarray:
    """All-vs-all RMSD matrix. Every group must have the same atom count
    and a consistent atom order (e.g. all "protein and atom CA" from
    structures/models of the same numbered protein)."""
    n = len(atom_groups)
    if n < 2:
        raise ClusterError("Need at least 2 conformations to cluster")
    counts = {len(g) for g in atom_groups}
    if len(counts) != 1:
        raise ClusterError(
            f"All conformations must select the same number of atoms to compare "
            f"(got sizes {sorted(counts)}); use a --selection that matches across all of them."
        )

    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            value = geom.rmsd(atom_groups[i], atom_groups[j], fit=fit)
            matrix[i, j] = matrix[j, i] = value
    return matrix


def _medoid(indices: list[int], matrix: np.ndarray) -> int:
    if len(indices) == 1:
        return indices[0]
    best, best_avg = indices[0], float("inf")
    for i in indices:
        avg = float(np.mean([matrix[i, j] for j in indices if j != i]))
        if avg < best_avg:
            best, best_avg = i, avg
    return best


@dataclass
class Cluster:
    id: int
    member_indices: list[int]
    member_labels: list[str]
    representative_label: str
    representative_index: int


@dataclass
class ClusterResult:
    method: str
    clusters: list[Cluster]
    matrix: np.ndarray
    labels: list[str]


def greedy(labels: list[str], matrix: np.ndarray, threshold: float = 2.0) -> ClusterResult:
    """CD-HIT-style incremental clustering: process conformations in
    order, join the first existing cluster whose seed (first member) is
    within `threshold` RMSD, else start a new cluster."""
    n = len(labels)
    seeds: list[int] = []
    members: list[list[int]] = []

    for i in range(n):
        joined = False
        for c, seed in enumerate(seeds):
            if matrix[i, seed] <= threshold:
                members[c].append(i)
                joined = True
                break
        if not joined:
            seeds.append(i)
            members.append([i])

    clusters = []
    for cluster_id, indices in enumerate(members, start=1):
        rep = _medoid(indices, matrix)
        clusters.append(
            Cluster(
                id=cluster_id, member_indices=indices,
                member_labels=[labels[i] for i in indices],
                representative_label=labels[rep], representative_index=rep,
            )
        )
    return ClusterResult(method="greedy", clusters=clusters, matrix=matrix, labels=labels)


class ClusterExtraNotAvailableError(RuntimeError):
    pass


def hierarchical(
    labels: list[str],
    matrix: np.ndarray,
    n_clusters: int | None = None,
    distance_threshold: float | None = None,
    linkage_method: str = "average",
) -> ClusterResult:
    """Agglomerative clustering via scipy.cluster.hierarchy. Exactly one
    of n_clusters / distance_threshold must be given."""
    try:
        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import squareform
    except ImportError as exc:
        raise ClusterExtraNotAvailableError(
            "Hierarchical clustering needs scipy. Install it with "
            "`pip install -e \".[cluster]\"` (or `pip install scipy`), "
            "or use --method greedy which has no extra dependency."
        ) from exc

    if (n_clusters is None) == (distance_threshold is None):
        raise ClusterError("Provide exactly one of n_clusters or distance_threshold")

    n = len(labels)
    condensed = squareform(matrix, checks=False)
    Z = linkage(condensed, method=linkage_method)

    if n_clusters is not None:
        assignments = fcluster(Z, t=n_clusters, criterion="maxclust")
    else:
        assignments = fcluster(Z, t=distance_threshold, criterion="distance")

    groups: dict[int, list[int]] = {}
    for i, cluster_num in enumerate(assignments):
        groups.setdefault(int(cluster_num), []).append(i)

    clusters = []
    for cluster_id, (_, indices) in enumerate(sorted(groups.items()), start=1):
        rep = _medoid(indices, matrix)
        clusters.append(
            Cluster(
                id=cluster_id, member_indices=indices,
                member_labels=[labels[i] for i in indices],
                representative_label=labels[rep], representative_index=rep,
            )
        )
    return ClusterResult(method="hierarchical", clusters=clusters, matrix=matrix, labels=labels)
