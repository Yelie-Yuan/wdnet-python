from typing import List, Tuple, Dict, Optional, Union
from pandas import DataFrame
from igraph import Graph
from copy import deepcopy
import numpy as np

# import warnings
from ._utils import node_strength_cy


def weighted_cor(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> Union[float, None]:
    """
    Compute the weighted correlation between two one-dimensional arrays.

    Parameters:
    ----------
    x: np.ndarray
        One-dimensional array.
    y: np.ndarray
        One-dimensional array.
    w: np.ndarray
        One-dimensional array of weights.

    Returns:
    -------
    Union[float, None]:
        Weighted correlation between x and y or None if the variances are zero.
    """
    # Ensure input arrays are numpy arrays with dtype float64 for higher precision
    x, y, w = map(lambda arr: np.asarray(arr, dtype=np.float64), (x, y, w))

    # Calculate the weighted means
    sum_w = np.sum(w)
    mean_x = np.sum(w * x) / sum_w
    mean_y = np.sum(w * y) / sum_w

    # Centralize the arrays
    x_centralized = x - mean_x
    y_centralized = y - mean_y

    # Calculate the weighted covariance and variances
    cov_xy = np.sum(w * x_centralized * y_centralized) / sum_w
    var_x = np.sum(w * x_centralized**2) / sum_w
    var_y = np.sum(w * y_centralized**2) / sum_w

    # Check for zero variance or NaN values
    if var_x == 0 or var_y == 0 or np.isnan(var_x) or np.isnan(var_y):
        return None

    return cov_xy / np.sqrt(var_x * var_y)


def is_WDNet(netwk: "WDNet") -> bool:
    # Check for required attributes
    nnode, nedge = check_wdnet_inputs(
        edgelist=netwk.edgelist,
        edgeweight=netwk.edge_attr["weight"].to_numpy(),
        directed=netwk.directed,
        edge_attr=netwk.edge_attr,
        node_attr=netwk.node_attr,
    )
    # Check if the node strengths (or degrees) are consistent with edgelist and edgeweight
    outs, ins = node_strength_cy(
        edgelist=netwk.edgelist,
        edgeweight=netwk.edge_attr["weight"].to_numpy(dtype=np.float64),
        nnode=nnode,
        nedge=nedge,
    )
    if netwk.directed:
        assert np.allclose(outs, netwk.node_attr["outs"].to_numpy(dtype=np.float64))
        assert np.allclose(ins, netwk.node_attr["ins"].to_numpy(dtype=np.float64))
    else:
        assert np.allclose(outs + ins, netwk.node_attr["s"].to_numpy(dtype=np.float64))
    return True


def check_wdnet_inputs(
    edgelist: np.ndarray,
    edgeweight: Optional[np.ndarray],
    directed: bool,
    edge_attr: Optional[DataFrame],
    node_attr: Optional[DataFrame],
) -> Tuple[int, int]:
    # Check for directed to be bool
    if not isinstance(directed, bool):
        raise ValueError("The 'directed' attribute must be of type bool.")
    # Check for edgelist to have at least one element
    if edgelist.shape[0] == 0:
        raise ValueError("'edgelist' must contain at least one edge.")

    if edgelist.shape[1] != 2:
        raise ValueError("'edgelist' must be a 2-column array.")
    tmp = set(edgelist.flatten())
    nnode = len(tmp)
    nedge = edgelist.shape[0]
    # Nodes must be consecutive integers starting from 0
    if not (nnode == max(tmp) + 1) or not (edgelist.dtype in [np.int_, np.int64]):
        raise ValueError("Nodes must be consecutive integers starting from 0.")
    # Check if the lengths of edgelist and edgeweight are the same
    if edgeweight is not None:
        if edgeweight.shape[0] != nedge:
            raise ValueError(
                "Number of elements in 'edgeweight' must be the same as 'edgelist'."
            )
        if not all(edgeweight > 0):
            raise ValueError("All elements in 'edgeweight' must be positive.")
    # Check if node_attr is a DataFrame or None
    if node_attr is not None:
        if not isinstance(node_attr, DataFrame):
            raise ValueError("'node_attr' must be of type DataFrame or None.")
        if node_attr.shape[0] != nnode:
            raise ValueError(
                "Number of rows in 'node_attr' should be the same as number of nodes."
            )

    # Check if edge_attr is a DataFrame or None
    if edge_attr is not None:
        if not isinstance(edge_attr, DataFrame):
            raise ValueError('"edge_attr" must be of type DataFrame or None.')
        if edge_attr.shape[0] != nedge:
            raise ValueError(
                "Number of rows in 'edge_attr' must be the same as number of edges."
            )
    return nnode, nedge


class WDNet:
    """
    Weighted Directed Network (WDNet) class.

    Attributes:
    ----------
    edgelist: np.ndarray
        List of edges represented as tuples (u, v) or a 2-column numpy
        array.
    directed: bool
        Indicates whether the network is directed or undirected.
        Defaults to True.
    edgeweight: np.ndarray
        List or numpy array of weights associated with each edge. If
        None, all edges are assumed to have weight 1. Defaults to
        None.
    edge_attr: Optional[DataFrame]
        Edge attributes stored as a DataFrame. Defaults to None.
    node_attr: Optional[DataFrame]
        Node attributes stored as a DataFrame. Defaults to None.
    """

    def __init__(
        self,
        edgelist: np.ndarray = np.array([[0, 1]]),
        directed: bool = True,
        edgeweight: Optional[np.ndarray] = None,
        edge_attr: Optional[DataFrame] = None,
        node_attr: Optional[DataFrame] = None,
    ):
        # Check for required attributes
        nnode, nedge = check_wdnet_inputs(
            edgelist=edgelist,
            edgeweight=edgeweight,
            directed=directed,
            edge_attr=edge_attr,
            node_attr=node_attr,
        )
        # Update attributes: directed, weighted, edgeweight and edgelist
        if edgeweight is None:
            edgeweight = np.ones(shape=(nedge,), dtype=np.int_)
            weighted = False
        else:
            if np.all(edgeweight == 1):
                weighted = False
                edgeweight = edgeweight.astype(np.int_)
            else:
                weighted = True
                edgeweight = edgeweight.astype(np.float64)
        edgelist = edgelist.astype(np.int_)
        self.edgelist = edgelist
        self.directed = directed
        self.weighted = weighted
        del weighted, edgelist, directed

        # Combine edgeweight into edge_attr
        self.edge_attr = DataFrame(index=range(nedge))
        self.edge_attr["weight"] = edgeweight
        if edge_attr is not None:
            if "weight" in edge_attr.columns:
                print("Edge weight has overwritten 'weight' in 'edge_attr'.")
            self.edge_attr = edge_attr.combine_first(self.edge_attr)

        # Compute node strengths (or degrees) and store them in node_attr
        outs, ins = node_strength_cy(
            self.edgelist, edgeweight.astype(np.float64), nnode, nedge
        )
        if not self.weighted:
            outs = outs.astype(np.int_)
            ins = ins.astype(np.int_)
        self.node_attr = DataFrame(index=range(nnode))
        if self.directed:
            self.node_attr["outs"] = outs
            self.node_attr["ins"] = ins
            if node_attr is not None:
                if "outs" in node_attr.columns:
                    print("Out-strength has overwritten 'outs' in 'node_attr'.")
                    node_attr = node_attr.drop(columns=["outs"])
                if "ins" in node_attr.columns:
                    print("In-strength has overwritten 'ins' in 'node_attr'.")
                    node_attr = node_attr.drop(columns=["ins"])
                self.node_attr = self.node_attr.combine_first(node_attr)
        else:
            self.node_attr["s"] = outs + ins
            if node_attr is not None:
                if "s" in node_attr.columns:
                    print("Strength has overwritten 's' column in 'node_attr'.")
                    node_attr = node_attr.drop(columns=["s"])
                self.node_attr = self.node_attr.combine_first(node_attr)

    def assortcoef(self) -> Union[Dict[str, Union[float, None]], Union[float, None]]:
        """
        Calculate the assortativity coefficient of the network.

        Returns:
        -------
        Union[Dict[str, Union[float, None]], Union[float, None]]
            A dictionary of floats for directed networks and a single
            float for undirected networks.

        Examples:
        --------
        >>> wdnet = WDNet(edgelist=[(1, 2), (2, 3), (3, 4)], directed=False)
        >>> wdnet.assortcoef()
        """
        src, tar = self.edgelist.T
        w = self.edge_attr["weight"].to_numpy()
        if self.directed:
            outs = self.node_attr["outs"].to_numpy()
            ins = self.node_attr["ins"].to_numpy()
            return {
                "outout": weighted_cor(
                    x=outs[src],
                    y=outs[tar],
                    w=w,
                ),
                "outin": weighted_cor(
                    x=outs[src],
                    y=ins[tar],
                    w=w,
                ),
                "inout": weighted_cor(
                    x=ins[src],
                    y=outs[tar],
                    w=w,
                ),
                "inin": weighted_cor(
                    x=ins[src],
                    y=ins[tar],
                    w=w,
                ),
            }
        else:
            s = self.node_attr["s"].to_numpy()
            return weighted_cor(
                x=s[np.concatenate([src, tar])],
                y=s[np.concatenate([tar, src])],
                w=np.concatenate([w, w]),
            )

    def centrality(self):
        # Method implementation
        print("This method is not implemented yet.")
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
            edgeweight=self.edge_attr["weight"].to_numpy(
                dtype=np.float64 if self.weighted else np.int_
            ),
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
            edgeweight=self.edge_attr["weight"].to_numpy(
                dtype=np.float64 if self.weighted else np.int_
            ),
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
            node_attr=self.node_attr.drop(
                columns=["outs", "ins"] if self.directed else ["s"]
            ),
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
            edges=[tuple(row) for row in self.edgelist],
            directed=self.directed,
        )
        for attr in self.edge_attr.columns:
            g.es[attr] = self.edge_attr[attr].tolist()
        for attr in self.node_attr.columns:
            g.vs[attr] = self.node_attr[attr].tolist()
        return g

    @classmethod
    def from_igraph(cls, g) -> "WDNet":
        """
        Convert from an igraph.Graph object.

        Parameters:
        ----------
        g (igraph.Graph): igraph.Graph object.
        """
        # Convert from an igraph.Graph object
        if not isinstance(g, Graph):
            raise ValueError("The input must be an igraph.Graph object.")
        edgelist = np.array(g.get_edgelist(), dtype=np.int_)
        directed = g.is_directed()
        if "weight" in g.es.attributes():
            edgeweight = np.array(g.es["weight"], dtype=np.float64)
        else:
            edgeweight = None
        edge_attr = DataFrame(index=range(g.ecount()))
        for attr in g.es.attributes():
            edge_attr[attr] = g.es[attr]
        node_attr = DataFrame(index=range(g.vcount()))
        for attr in g.vs.attributes():
            node_attr[attr] = g.vs[attr]
        netwk = WDNet(
            edgelist=edgelist,
            directed=directed,
            edgeweight=edgeweight,
            edge_attr=edge_attr,
            node_attr=node_attr,
        )

        is_WDNet(netwk)
        return netwk

    @classmethod
    def from_adjacency(
        cls, adj: np.ndarray, directed: bool = True, weighted: bool = True
    ) -> "WDNet":
        """
        Convert from an adjacency matrix represented as an n by n numpy array,
        where n is the number of nodes.

        Parameters:
        ----------
        adj (np.ndarray): Adjacency matrix represented as an n by n numpy array.
        directed (bool): Indicates whether the network is directed or undirected.
            Defaults to True.
        weighted (bool): Indicates whether the network is weighted or not.
            Defaults to True. If False, all non-zero elements in the adjacency matrix
            are considered to be the number of edges between the corresponding nodes.

        Returns:
        --------
        WDNet:
            A WDNet object created from the given adjacency matrix.
        """
        nnode = adj.shape[0]
        if nnode != adj.shape[1]:
            raise ValueError("The adjacency matrix must be square.")

        if not directed:
            if not np.allclose(adj, adj.T):
                raise ValueError(
                    "The adjacency matrix must be symmetric for undirected networks."
                )
            adj = np.triu(adj)

        if not weighted:
            if not np.array_equal(adj, adj.astype(int)):
                raise ValueError(
                    "The adjacency matrix must be integer-valued for unweighted networks."
                )

        if not np.all(adj >= 0):
            raise ValueError(
                "All elements in the adjacency matrix must be non-negative."
            )

        i, j = np.nonzero(adj)
        edgeweight = adj[i, j]

        if not weighted:
            edgelist = np.array(
                [(x, y) for x, y, w in zip(i, j, adj[i, j]) for _ in range(w)],
                dtype=np.int_,
            )
            edgeweight = np.ones(len(edgelist), dtype=np.int_)
        else:
            edgelist = np.column_stack((i, j))

        netwk = cls(
            edgelist=edgelist,
            directed=directed,
            edgeweight=edgeweight,
        )

        # Assuming is_WDNet is a function to validate the network
        is_WDNet(netwk)
        return netwk

    def to_edgelist(self) -> Tuple[np.ndarray, np.ndarray]:
        return self.edgelist, self.edge_attr["weight"].to_numpy(
            dtype=np.float64 if self.weighted else np.int_
        )

    def save_edgelist(self, file: str) -> None:
        edgelist_w = np.column_stack(
            (
                self.edgelist,
                self.edge_attr["weight"].to_numpy(
                    dtype=np.float64 if self.weighted else np.int_
                ),
            )
        )
        if self.weighted:
            np.savetxt(file, edgelist_w, delimiter=",", fmt=("%d", "%d", "%f"))
        else:
            np.savetxt(file, edgelist_w, delimiter=",", fmt=("%d", "%d", "%d"))

    def to_adjacency(self) -> np.ndarray:
        nnode = self.node_attr.shape[0]
        adj = np.zeros(shape=(nnode, nnode))

        for i, (src, dst) in enumerate(self.edgelist):
            adj[src, dst] = self.edge_attr["weight"][i]

        if not self.directed:
            tmp_diag = np.diag(adj).copy()
            np.fill_diagonal(adj, 0)
            adj += adj.T
            np.fill_diagonal(adj, tmp_diag)

        return adj

    def save_adjacency(self, file: str) -> None:
        np.savetxt(
            file,
            self.to_adjacency(),
            delimiter=",",
            fmt="%f" if self.weighted else "%d",
        )

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
        return (
            f"WDNet object: directed={self.directed}, "
            f"weighted={self.weighted}, "
            f"number of nodes={self.node_attr.shape[0]}, "
            f"number of edges={len(self.edgelist)}."
        )
