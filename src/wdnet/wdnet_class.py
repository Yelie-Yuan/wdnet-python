from typing import List, Tuple, Dict, Optional, Union
from pandas import DataFrame
import numpy as np
import warnings
from ._utils import node_strength_py


def weighted_cor(x, y, w):
    # Calculate the weighted means
    mean_x = np.average(x, weights=w)
    mean_y = np.average(y, weights=w)

    # Calculate the weighted covariance and variances
    cov_xy = np.sum(w * (x - mean_x) * (y - mean_y)) / np.sum(w)
    var_x = np.sum(w * (x - mean_x) ** 2) / np.sum(w)
    var_y = np.sum(w * (y - mean_y) ** 2) / np.sum(w)
    if var_x == 0 or var_y == 0 or np.isnan(var_x) or np.isnan(var_y):
        return None

    return cov_xy / np.sqrt(var_x * var_y)


def is_WDNet(netwk):
    check_wdnet_inputs(
        edgelist=netwk.edgelist,
        edgeweight=netwk.edge_attr["weight"],
        directed=netwk.directed,
        edge_attr=netwk.edge_attr,
        node_attr=netwk.node_attr,
    )
    outs, ins = node_strength_py(netwk.edgelist, netwk.edge_attr["weight"].to_numpy())
    if netwk.directed:
        difference = np.sum(np.abs(outs - netwk.node_attr["outs"].to_numpy())) + np.sum(
            np.abs(ins - netwk.node_attr["ins"])
        )
    else:
        difference = np.sum(np.abs(outs + ins - netwk.node_attr["s"].to_numpy()))

    if difference > 1e-5:
        raise ValueError(
            "The node strengths (or degrees) in 'node_attr' are not consistent with 'edgelist' and 'edgeweight'."
        )
    return True


def check_wdnet_inputs(edgelist, edgeweight, directed, edge_attr, node_attr):
    # Check for directed to be bool
    if not isinstance(directed, bool):
        raise ValueError("The 'directed' attribute must be of type bool.")

    # Check for edgelist to have at least one element
    if len(edgelist) < 1:
        raise ValueError("'edgelist' must contain at least one element.")

    # Check if all elements in each tuple in edgelist are nonnegative integers
    if not all(
        len(edge) == 2 and all(isinstance(i, int) and i >= 0 for i in edge)
        for edge in edgelist
    ):
        raise ValueError("Each edge must be represented by two nonnegative integers.")

    # Check if the lengths of edgelist and edgeweight are the same
    if edgeweight is not None:
        if len(edgelist) != len(edgeweight):
            raise ValueError("The lengths of edgelist and edgeweight must be the same.")
        if not all(x > 0 for x in edgeweight):
            raise ValueError("All elements in 'edgeweight' must be greater than zero.")

    # Check if edge_attr is a DataFrame or None
    if edge_attr is not None:
        if not isinstance(edge_attr, DataFrame):
            raise ValueError(
                "The 'edge_attr' attribute must be of type DataFrame or None."
            )
        if edge_attr.shape[0] != len(edgelist):
            raise ValueError(
                "The number of rows in 'edge_attr' must be the same as the length of 'edgelist'."
            )

    # Check if node_attr is a DataFrame or None
    if node_attr is not None:
        if not isinstance(node_attr, DataFrame):
            raise ValueError(
                "The 'node_attr' attribute must be of type DataFrame or None."
            )
        max_node = max(
            x for edge in edgelist for x in edge
        )  # maximum node value in edgelist
        if node_attr.shape[0] != max_node + 1:
            raise ValueError(
                "The number of rows in 'node_attr' should be (maximum value in 'edgelist') + 1."
            )

    return True


class WDNet:
    """
    Weighted Directed Network (WDNet) class for representing networks.

    Attributes:
        edgelist (Union[List[Tuple[int]], np.ndarray]): List of edges represented as tuples (u, v) or a 2-column numpy array.
        directed (bool): Indicates whether the network is directed or undirected. Defaults to True.
        edgeweight (Union[List[float], np.ndarray]): List or numpy array of weights associated with each edge. If None, all edges are assumed to have weight 1. Defaults to None.
        edge_attr (DataFrame): Edge attributes stored as a DataFrame. Defaults to None.
        node_attr (DataFrame): Node attributes stored as a DataFrame. Defaults to None.
        weighted (bool): Indicates whether the network is weighted or not.
    """

    def __init__(
        self,
        edgelist: Union[List[Tuple[int]], np.ndarray] = [(1, 2)],
        directed: bool = True,
        edgeweight: Optional[Union[List[float], np.ndarray]] = None,
        edge_attr: Optional[DataFrame] = None,
        node_attr: Optional[DataFrame] = None,
    ):
        if isinstance(edgelist, np.ndarray):
            edgelist = edgelist.tolist()
        if isinstance(edgeweight, np.ndarray):
            edgeweight = edgeweight.tolist()
        # Check for required attributes
        check_wdnet_inputs(
            edgelist=edgelist,
            edgeweight=edgeweight,
            directed=directed,
            edge_attr=edge_attr,
            node_attr=node_attr,
        )
        # Update attributes: directed, weighted, edgeweight and edgelist
        self.directed = directed
        if edgeweight is None:
            edgeweight = [1] * len(edgelist)
            self.weighted = False
        else:
            self.weighted = all(x == 1 for x in edgeweight)
        self.edgelist = edgelist

        # Combine edgeweight into edge_attr
        self.edge_attr = DataFrame(index=range(len(edgelist)))
        self.edge_attr["weight"] = edgeweight

        if edge_attr is not None:
            if "weight" in edge_attr.columns:
                warnings.warn(
                    "'edgeweight' has overwritten the 'weight' column in 'edge_attr'."
                )
            self.edge_attr = edge_attr.combine_first(self.edge_attr)

        # Compute node strengths (or degrees) and store them in node_attr
        outs, ins = node_strength_py(edgelist, edgeweight)
        self.node_attr = DataFrame(index=range(len(outs)))
        if self.directed:
            self.node_attr["outs"] = outs
            self.node_attr["ins"] = ins
            if node_attr is not None:
                if "outs" in node_attr.columns:
                    warnings.warn(
                        "Out-strength has overwritten the 'outs' column in 'node_attr'."
                    )
                    node_attr = node_attr.drop(columns=["outs"])
                if "ins" in node_attr.columns:
                    warnings.warn(
                        "In-strength has overwritten the 'ins' column in 'node_attr'."
                    )
                    node_attr = node_attr.drop(columns=["ins"])
                self.node_attr = self.node_attr.combine(node_attr)
        else:
            self.node_attr["s"] = outs + ins
            if node_attr is not None:
                if "s" in node_attr.columns:
                    warnings.warn(
                        "Strength has overwritten the 's' column in 'node_attr'."
                    )
                    node_attr = node_attr.drop(columns=["s"])
                self.node_attr = self.node_attr.combine(node_attr)

    def is_WDNet(self):
        # Check for required attributes
        if not all(hasattr(self, attr) for attr in ["edgelist", "edgeweight"]):
            return False

        # Validate that the lengths of 'edgelist' and 'edgeweight' are the same
        if len(self.edgelist) != len(self.edgeweight):
            return False

        return True

    def assortcoef(self) -> Union[List[float], float]:
        """
        Calculate the assortativity coefficient of the network.

        Returns:
        -------
        Union[Dict[float], float]
            A dictionary of floats for directed networks and a single float for undirected networks.

        Examples:
        --------
        >>> wdnet = WDNet(edgelist=[(1, 2), (2, 3), (3, 4)], directed=False)
        >>> wdnet.assortcoef()
        """
        src, tar = map(list, zip(*self.edgelist))
        if self.directed:
            return {
                "outout": weighted_cor(
                    self.node_attr.loc[src, "outs"].to_numpy(),
                    self.node_attr.loc[tar, "outs"].to_numpy(),
                    self.edge_attr["weight"].to_numpy(),
                ),
                "outin": weighted_cor(
                    self.node_attr.loc[src, "outs"].to_numpy(),
                    self.node_attr.loc[tar, "ins"].to_numpy(),
                    self.edge_attr["weight"].to_numpy(),
                ),
                "inout": weighted_cor(
                    self.node_attr.loc[src, "ins"].to_numpy(),
                    self.node_attr.loc[tar, "outs"].to_numpy(),
                    self.edge_attr["weight"].to_numpy(),
                ),
                "inin": weighted_cor(
                    self.node_attr.loc[src, "ins"].to_numpy(),
                    self.node_attr.loc[tar, "ins"].to_numpy(),
                    self.edge_attr["weight"].to_numpy(),
                ),
            }
        else:
            return weighted_cor(
                self.node_attr.loc[src + tar, "s"].to_numpy(),
                self.node_attr.loc[tar + src, "s"].to_numpy(),
                np.concatenate(
                    [
                        self.edge_attr["weight"].to_numpy(),
                        self.edge_attr["weight"].to_numpy(),
                    ]
                ),
            )

    def centrality(self):
        # Method implementation
        pass

    def to_undirected(self):
        if not self.directed:
            return self

        self.directed = False
        self.node_attr["s"] = (
            self.node_attr["outs"].to_numpy() + self.node_attr["ins"].to_numpy()
        )
        self.node_attr = self.node_attr.drop(columns=["outs", "ins"])

        is_WDNet(self)

        return self

    def to_unweighted(self):
        self.weighted = False
        self.edge_attr["weight"] = [1] * len(self.edgelist)
        outs, ins = node_strength_py(self.edgelist, self.edge_attr["weight"].to_numpy())
        if self.directed:
            self.node_attr["outs"] = outs
            self.node_attr["ins"] = ins
        else:
            self.node_attr["s"] = outs + ins
        
        is_WDNet(self)

        return self

    def to_igraph(self):
        # Convert to an igraph.Graph object
        pass

    @classmethod
    def from_igraph(self, igraph_obj):
        # Convert from an igraph.Graph object
        pass

    @classmethod
    def from_adj(self, adj, directed=False, weighted=True):
        # Convert from an adjacency matrix
        pass

    def to_edgelist(self, file=None):
        # Convert to an edgelist
        pass

    def to_adj(self, file=None):
        adj = self.to_adj()

        if file is None:
            print(adj)
        else:
            # Save to file
            # np.save(file, adj, delimiter=",")
            pass

    def __str__(self):
        return f"WDNet Object: directed={self.directed}, weighted={self.weighted}"