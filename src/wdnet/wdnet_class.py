from typing import List, Tuple, Dict, Optional, Union
from pandas import DataFrame
from igraph import Graph
from copy import deepcopy
import numpy as np

# import warnings
from .utils_cy import node_strength_cy


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
    # Check for required attributes
    check_wdnet_inputs(
        edgelist=netwk.edgelist,
        edgeweight=netwk.edge_attr["weight"],
        directed=netwk.directed,
        edge_attr=netwk.edge_attr,
        node_attr=netwk.node_attr,
    )
    # Check if the node strengths (or degrees) are consistent with edgelist and edgeweight
    outs, ins = node_strength_cy(netwk.edgelist, netwk.edge_attr["weight"].to_numpy())
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

    if not all(len(edge) == 2 for edge in edgelist):
        raise ValueError("All edges in 'edgelist' must have length 2.")
    tmp = set(x for edge in edgelist for x in edge)
    nnode = len(tmp)
    # Nodes must be consecutive integers starting from 0
    if not nnode == np.max(list(tmp)) + 1 or not all(isinstance(i, int) for i in tmp):
        raise ValueError("Nodes must be consecutive integers starting from 0.")

    # Check if node_attr is a DataFrame or None
    if node_attr is not None:
        if not isinstance(node_attr, DataFrame):
            raise ValueError(
                "The 'node_attr' attribute must be of type DataFrame or None."
            )
        if node_attr.shape[0] != nnode:
            raise ValueError(
                "The number of rows in 'node_attr' should be (maximum value in 'edgelist') + 1."
            )

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
            self.weighted = not all(x == 1 for x in edgeweight)
        self.edgelist = edgelist

        # Combine edgeweight into edge_attr
        self.edge_attr = DataFrame(index=range(len(edgelist)))
        self.edge_attr["weight"] = edgeweight

        if edge_attr is not None:
            if "weight" in edge_attr.columns:
                print(
                    "Edge weight has overwritten 'weight' in 'edge_attr'."
                )
            self.edge_attr = edge_attr.combine_first(self.edge_attr)

        # Compute node strengths (or degrees) and store them in node_attr
        outs, ins = node_strength_cy(edgelist, edgeweight)
        self.node_attr = DataFrame(index=range(len(outs)))
        if self.directed:
            self.node_attr["outs"] = outs
            self.node_attr["ins"] = ins
            if node_attr is not None:
                if "outs" in node_attr.columns:
                    print(
                        "Out-strength has overwritten 'outs' in 'node_attr'."
                    )
                    node_attr = node_attr.drop(columns=["outs"])
                if "ins" in node_attr.columns:
                    print(
                        "In-strength has overwritten 'ins' in 'node_attr'."
                    )
                    node_attr = node_attr.drop(columns=["ins"])
                self.node_attr = self.node_attr.combine_first(node_attr)
        else:
            self.node_attr["s"] = outs + ins
            if node_attr is not None:
                if "s" in node_attr.columns:
                    print("Strength has overwritten 's' column in 'node_attr'.")
                    node_attr = node_attr.drop(columns=["s"])
                self.node_attr = self.node_attr.combine_first(node_attr)

    def is_WDNet(self) -> bool:
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

    def to_directed(self) -> "WDNet":
        """
        Convert to a directed network.

        Returns:
        -------
        WDNet
        """
        if self.directed:
            print("The network is already directed.")
            return self.copy()
        
        netwk = WDNet(
            edgelist=self.edgelist,
            directed=True,
            edgeweight=self.edge_attr["weight"],
            edge_attr=self.edge_attr.drop(columns=["weight"]),
            node_attr=self.node_attr.drop(columns=["s"]),
        )
        is_WDNet(netwk)
        return netwk

    def to_undirected(self) -> "WDNet":
        """
        Convert to an undirected network.

        Returns:
        -------
        WDNet
        """
        if not self.directed:
            print("The network is already undirected.")
            return self.copy()

        netwk = WDNet(
            edgelist=self.edgelist,
            directed=False,
            edgeweight=self.edge_attr["weight"],
            edge_attr=self.edge_attr.drop(columns=["weight"]),
            node_attr=self.node_attr.drop(columns=["outs", "ins"]),
        )
        is_WDNet(netwk)
        return netwk

    def to_unweighted(self) -> "WDNet":
        """
        Convert to an unweighted network.

        Returns:
        -------
        WDNet
        """
        if not self.weighted:
            print("The network is already unweighted.")
            return self.copy()

        netwk = WDNet(
            edgelist=self.edgelist,
            directed=self.directed,
            edgeweight=None,
            edge_attr=self.edge_attr.drop(columns=["weight"]),
            node_attr=self.node_attr.drop(columns=["outs", "ins"] if self.directed else ["s"]),
        )
        is_WDNet(netwk)
        return netwk

    def to_igraph(self) -> Graph:
        """
        Convert to an igraph.Graph object.

        Returns:
        -------
        igraph.Graph
        """
        g = Graph.TupleList(
            edges=self.edgelist,
            directed=self.directed,
        )
        for attr in self.edge_attr.columns:
            g.es[attr] = self.edge_attr[attr].tolist()
        for attr in self.node_attr.columns:
            g.vs[attr] = self.node_attr[attr].tolist()
        return g

    @classmethod
    def from_igraph(self, g) -> "WDNet":
        """
        Convert from an igraph.Graph object.

        Parameters:
        ----------
        g (igraph.Graph): igraph.Graph object.
        """
        # Convert from an igraph.Graph object
        if not isinstance(g, Graph):
            raise ValueError("The input must be an igraph.Graph object.")
        edgelist = g.get_edgelist()
        directed = g.is_directed()
        if "weight" in g.es.attributes():
            edgeweight = g.es["weight"]
        else:
            edgeweight = None
        edge_attr = DataFrame(index=range(len(edgelist)))
        for attr in g.es.attributes():
            edge_attr[attr] = g.es[attr]
        node_attr = DataFrame(index=range(g.vcount()))
        for attr in g.vs.attributes():
            node_attr[attr] = g.vs[attr]
        return WDNet(
            edgelist=edgelist,
            directed=directed,
            edgeweight=edgeweight,
            edge_attr=edge_attr,
            node_attr=node_attr,
        )

    @classmethod
    def from_adjacency(self, adj, directed=True, weighted=True) -> "WDNet":
        """
        Convert from an adjacency matrix represented as a n by n numpy array,
        where n is the number of nodes.

        Parameters:
        ----------
        adj (np.ndarray): Adjacency matrix represented as a n by n numpy array.
        directed (bool): Indicates whether the network is directed or undirected. Defaults to True.
        weighted (bool): Indicates whether the network is weighted or not. Defaults to True. If
        False, all non-zero elements in the adjacency matrix are considered to be the number of edges
        between the corresponding nodes.
        """
        if not adj.shape[0] == adj.shape[1]:
            raise ValueError("The adjacency matrix must be square.")
        if not directed and not np.allclose(adj, adj.T):
            raise ValueError(
                "The adjacency matrix must be symmetric for undirected networks."
            )
        if not weighted and not np.allclose(adj, adj.astype(int)):
            raise ValueError(
                "The adjacency matrix must be integer-valued for unweighted networks."
            )
        if not all(x >= 0 for x in adj.flatten()):
            raise ValueError(
                "All elements in the adjacency matrix must be non-negative."
            )
        if not directed:
            adj = np.triu(adj)
        edgelist = []
        edgeweight = []
        for i in range(adj.shape[0]):
            for j in range(adj.shape[1]):
                if adj[i, j] != 0:
                    if weighted:
                        edgeweight.append(adj[i, j])
                        edgelist.append((i, j))
                    else:
                        edgeweight.append([1] * adj[i, j])
                        edgelist.append([(i, j)] * adj[i, j])
        return WDNet(
            edgelist=edgelist,
            directed=directed,
            edgeweight=edgeweight,
        )

    def to_edgelist(self, file=None):
        if file is None:
            return self.edgelist, self.edge_attr["weight"].tolist()
        else:
            edgelist_w = [(*edge, self.edge_attr["weight"][i]) for i, edge in enumerate(self.edgelist)]
            if self.weighted:
                np.savetxt(file, edgelist_w, delimiter=",", fmt=('%d', '%d', '%f'))
            else:
                np.savetxt(file, edgelist_w, delimiter=",", fmt=('%d', '%d', '%d'))
        

    def to_adjacency(self, file=None):
        nnode = self.node_attr.shape[0]
        adj = np.zeros(shape=(nnode, nnode))
        for i, edge in enumerate(self.edgelist):
            adj[edge[0], edge[1]] = self.edge_attr["weight"][i]
        if not self.directed:
            tmp_diag = np.diag(adj).copy()
            np.fill_diagonal(adj, 0)
            adj = adj + adj.T
            np.fill_diagonal(adj, tmp_diag)
        if file is None:
            return adj
        else:
            np.savetxt(file, adj, delimiter=",")

    def copy(self) -> "WDNet":
        """
        Return a copy of the network.

        Returns:
        -------
        WDNet
        """
        return deepcopy(self)

    def __str__(self):
        """
        Return a summary of the network, including the number of nodes and edges.
        """
        return (f"WDNet object: directed={self.directed}, "
                f"weighted={self.weighted}, "
                f"number of nodes={self.node_attr.shape[0]}, "
                f"number of edges={len(self.edgelist)}.")
