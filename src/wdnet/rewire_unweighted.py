# Import required libraries
import numpy as np
import pandas as pd
import cvxpy as cp
import warnings
from typing import List, Tuple, Dict, Optional, Union, Callable
from pandas import DataFrame
from math import sqrt
from .wdnet_class import WDNet


def compute_correlation(x, y, xsum, ysum, x2sum, y2sum):
    n = len(x)
    xysum = np.sum(np.multiply(x, y))
    numerator = n * xysum - xsum * ysum
    denominator = sqrt((n * x2sum - xsum**2) * (n * y2sum - ysum**2))
    return numerator / denominator


def rewire_directed(
    iteration,
    nattempts,
    tnode,
    sout,
    sin,
    tout,
    tin,
    index_s,
    index_t,
    eta,
    history,
):
    outout = np.zeros(iteration)
    outin = np.zeros(iteration)
    inout = np.zeros(iteration)
    inin = np.zeros(iteration)

    soutsum, sinsum = np.sum(sout), np.sum(sin)
    toutsum, tinsum = np.sum(tout), np.sum(tin)
    sout2sum = np.sum(sout**2)
    sin2sum = np.sum(sin**2)
    tout2sum = np.sum(tout**2)
    tin2sum = np.sum(tin**2)

    nedge = len(tnode)
    hist_row = iteration * nattempts if history else 1

    rewire_history = np.zeros((hist_row, 4), dtype=int)
    count = 0

    for n in range(iteration):
        for i in range(nattempts):
            e1 = int(np.floor(np.random.random() * nedge))
            e2 = int(np.floor(np.random.random() * nedge))

            while e1 == e2:
                e2 = int(np.floor(np.random.random() * nedge))

            if history:
                rewire_history[count, :3] = [count + 1, e1, e2]

            s1, s2 = index_s[e1], index_s[e2]
            t1, t2 = index_t[e1], index_t[e2]

            ratio = min(1, eta[s1][t2] * eta[s2][t1] / (eta[s1][t1] * eta[s2][t2]))

            u = np.random.random()
            if u <= ratio:
                index_t[e1], index_t[e2] = index_t[e2], index_t[e1]
                tnode[e1], tnode[e2] = tnode[e2], tnode[e1]
                tout[e1], tout[e2] = tout[e2], tout[e1]
                tin[e1], tin[e2] = tin[e2], tin[e1]

                if history:
                    rewire_history[count, 3] = 1

            count += 1

        outout[n] = compute_correlation(
            sout, tout, soutsum, toutsum, sout2sum, tout2sum
        )
        outin[n] = compute_correlation(sout, tin, soutsum, tinsum, sout2sum, tin2sum)
        inout[n] = compute_correlation(sin, tout, sinsum, toutsum, sin2sum, tout2sum)
        inin[n] = compute_correlation(sin, tin, sinsum, tinsum, sin2sum, tin2sum)

    return tnode, outout, outin, inout, inin, rewire_history


# Compute the standard deviation of a vector j with probability vector q
def my_sigma(j, q):
    return np.sqrt(np.sum(j**2 * q) - np.sum(j * q) ** 2)


# Check if the target assortativity coefficient is valid
def check_target_assortcoef(netwk, target_assortcoef):
    if netwk.directed:
        if target_assortcoef:
            if type(target_assortcoef) != dict:
                raise ValueError(
                    "The target assortativity coefficient must be a dictionary for directed networks."
                )
            invalid_values = [
                key for key, value in target_assortcoef.items() if not -1 <= value <= 1
            ]
            if invalid_values:
                raise ValueError(
                    "The target assortativity coefficient must be in [-1, 1]"
                )
    else:
        if target_assortcoef:
            if type(target_assortcoef) != float:
                raise ValueError(
                    "The target assortativity coefficient must be a float for undirected networks."
                )
            if not -1 <= target_assortcoef <= 1:
                raise ValueError(
                    "The target assortativity coefficient must be in [-1, 1]."
                )


# Name the rows and columns of eta for directed networks
def name_eta(eta, d1, d2, index1, index2):
    eta = pd.DataFrame(eta)
    tmp = np.array([f"{i}-{j}" for i in d1 for j in d2])
    eta.index = tmp[index1]
    eta.columns = tmp[index2]
    return eta


# Get degree distributions of the given unweighted, directed network
def get_dist_directed(netwk):
    if netwk.weighted:
        raise ValueError("The network must be unweighted.")
    if not netwk.directed:
        raise ValueError("The network must be directed.")
    outd = netwk.node_attr["outs"].values.astype(int)
    ind = netwk.node_attr["ins"].values.astype(int)

    nu = pd.crosstab(outd, ind, normalize=True)

    snodes, tnodes = np.array(netwk.edgelist).T
    tmp = pd.DataFrame(
        {
            "i": outd[snodes],
            "j": ind[snodes],
            "k": outd[tnodes],
            "l": ind[tnodes],
        }
    )
    del snodes, tnodes, outd, ind
    eta = pd.crosstab(
        index=[tmp["i"], tmp["j"]], columns=[tmp["k"], tmp["l"]], normalize=True
    )
    eta.index = eta.index.map(lambda x: f"{x[0]}-{x[1]}")
    eta.columns = eta.columns.map(lambda x: f"{x[0]}-{x[1]}")

    e = {
        "outout": pd.crosstab(index=tmp["i"], columns=tmp["k"], normalize=True),
        "outin": pd.crosstab(index=tmp["i"], columns=tmp["l"], normalize=True),
        "inout": pd.crosstab(index=tmp["j"], columns=tmp["k"], normalize=True),
        "inin": pd.crosstab(index=tmp["j"], columns=tmp["l"], normalize=True),
    }

    dout = nu.index.values
    din = nu.columns.values
    pout = nu.sum(axis=1).values
    pin = nu.sum(axis=0).values
    soutin = nu.multiply(dout, axis=0)
    soutin /= soutin.values.sum()
    toutin = nu.multiply(din, axis=1)
    toutin /= toutin.values.sum()
    qsout = soutin.sum(axis=1).values
    qsin = soutin.sum(axis=0).values
    qtout = toutin.sum(axis=1).values
    qtin = toutin.sum(axis=0).values
    soutin = np.array(soutin).ravel(order="C")
    toutin = np.array(toutin).ravel(order="C")
    return {
        "nu": nu,  # pandas DataFrame
        "eta": eta,  # pandas DataFrame
        "e": e,  # a dictionary of pandas DataFrames
        "dout": dout,  # numpy array; one-dimensional
        "din": din,
        "pout": pout,  # numpy array; one-dimensional
        "pin": pin,
        "qsout": qsout,  # numpy array; one-dimensional
        "qsin": qsin,
        "qtout": qtout,
        "qtin": qtin,
        "soutin": soutin,  # numpy array; one-dimensional
        "toutin": toutin,
    }


# Get degree distributions of the given unweighted, undirected network
def get_dist_undirected(netwk):
    if netwk.weighted:
        raise ValueError("The network must be unweighted.")
    if netwk.directed:
        raise ValueError("The network must be undirected.")
    d = netwk.node_attr["s"].values.astype(int)
    nu = pd.crosstab(index=d, columns=d, normalize=True)
    p = nu.sum(axis=1).values
    tmp1, tmp2 = np.array(netwk.edgelist).T
    d1 = d[np.concatenate((tmp1, tmp2))]
    d2 = d[np.concatenate((tmp2, tmp1))]
    del tmp1, tmp2, d
    d = nu.index.values
    eta = pd.crosstab(index=d1, columns=d2, normalize=True)
    q = eta.sum(axis=1).values
    return {
        "d": d,
        "nu": nu,
        "eta": eta,
        "p": p,
        "q": q,
    }


# Construct an eta for directed networks
def get_eta_directed(
    netwk,
    target_assortcoef,
    eta_obj,
    which_range,
    **kwargs,
):
    """
    Construct an eta for directed networks or compute the range of
    assortativity coefficients.

    Parameters
    ----------
    netwk (WDNet): An instance of the WDNet class.

    target_assortcoef (Optional[Dict[str, float]]): Target assortativity coefficient(s).
        This should be a dictionary with keys 'outout', 'outin', 'inout',
        and 'inin' corresponding to the four types of assortativity
        coefficients.

    eta_obj (Optional[Callable]): A objective function that takes eta
        (a numpy array) as input and returns a scalar value. It will be
        minimized when solving for an appropriate eta.
        Defaults to 0 (i.e., no objective function).

    which_range (Optional[str]): Which range of assortativity
        coefficients to compute.

    **kwargs: Additional keyword arguments to be passed to cvxpy when solving
        the optimization problem.
    """
    dist = get_dist_directed(netwk)
    m = len(dist["dout"])
    n = len(dist["din"])
    soutin = dist["soutin"]
    toutin = dist["toutin"]
    indexs = soutin != 0
    indext = toutin != 0

    eta = cp.Variable((np.sum(indexs), np.sum(indext)), nonneg=True)
    constrs = [
        cp.sum(eta, axis=1) == soutin[indexs],
        cp.sum(eta, axis=0) == toutin[indext],
    ]
    del soutin, toutin

    mat1 = np.kron(np.eye(m), np.ones((1, n)))
    mat2 = np.kron(np.ones((m, 1)), np.eye(n))
    e = {
        "outout": mat1[:, indexs] @ eta @ mat1[:, indext].T,
        "outin": mat1[:, indexs] @ eta @ mat2[indext, :],
        "inout": mat2[indexs, :].T @ eta @ mat1[:, indext].T,
        "inin": mat2[indexs, :].T @ eta @ mat2[indext, :],
    }
    del mat1, mat2, m, n

    sig = {
        "sout": my_sigma(dist["dout"], dist["qsout"]),
        "sin": my_sigma(dist["din"], dist["qsin"]),
        "tout": my_sigma(dist["dout"], dist["qtout"]),
        "tin": my_sigma(dist["din"], dist["qtin"]),
    }
    r = {
        "outout": dist["dout"].T
        @ (e["outout"] - np.outer(dist["qsout"], dist["qtout"]))
        @ dist["dout"]
        / sig["sout"]
        / sig["tout"],
        "outin": dist["dout"].T
        @ (e["outin"] - np.outer(dist["qsout"], dist["qtin"]))
        @ dist["din"]
        / sig["sout"]
        / sig["tin"],
        "inout": dist["din"].T
        @ (e["inout"] - np.outer(dist["qsin"], dist["qtout"]))
        @ dist["dout"]
        / sig["sin"]
        / sig["tout"],
        "inin": dist["din"].T
        @ (e["inin"] - np.outer(dist["qsin"], dist["qtin"]))
        @ dist["din"]
        / sig["sin"]
        / sig["tin"],
    }
    if target_assortcoef is not None:
        for each in target_assortcoef.keys():
            if each in r.keys():
                if target_assortcoef[each] is not None:
                    constrs.append(r[each] == target_assortcoef[each])
    if which_range is None:
        problem = cp.Problem(cp.Minimize(eta_obj(eta)), constrs)
        problem.solve(**kwargs)
        if problem.status != "optimal":
            # show a warning about solver status
            warnings.warn(f"Solver status: {problem.status}")
        eta = name_eta(eta.value, dist["dout"], dist["din"], indexs, indext)
        return eta, problem.status, problem.solver_stats
    else:
        problem1 = cp.Problem(cp.Minimize(r[which_range]), constrs)
        problem1.solve(**kwargs)
        if problem1.status != "optimal":
            # show a warning about solver status
            warnings.warn(f"Lower bound solver status: {problem1.status}")
        problem2 = cp.Problem(cp.Maximize(r[which_range]), constrs)
        problem2.solve(**kwargs)
        if problem2.status != "optimal":
            # show a warning about solver status
            warnings.warn(f"Upper bound solver status: {problem2.status}")
        return (
            [problem1.value, problem2.value],
            [problem1.status, problem1.solver_stats],
            [problem2.status, problem2.solver_stats],
        )


# Construct an eta for undirected networks
def get_eta_undirected(
    netwk,
    target_assortcoef,
    eta_obj,
    **kwargs,
):
    """
    Construct an eta for undirected networks.

    Parameters
    ----------
    netwk (WDNet): An instance of the WDNet class.

    target_assortcoef (Optional[float]): Target assortativity coefficient.
        This should be a single float value.

    eta_obj (Optional[Callable]): A objective function that takes eta
        (a numpy array) as input and returns a scalar value. It will be
        minimized when solving for an appropriate eta.
        Defaults to 0 (i.e., no objective function).

    which_range (Optional[str]): Which range of assortativity
        coefficients to compute.

    **kwargs: Additional keyword arguments to be passed to cvxpy when solving
        the optimization problem.
    """
    dist = get_dist_undirected(netwk)
    m = len(dist["d"])
    eta = cp.Variable((m, m), nonneg=True)
    constrs = [cp.sum(eta, axis=1) == dist["q"], eta == eta.T]
    del m
    sig2 = np.sum(dist["d"] ** 2 * dist["q"]) - np.sum(dist["d"] * dist["q"]) ** 2
    r = dist["d"].T @ (eta - np.outer(dist["q"], dist["q"])) @ dist["d"] / sig2
    if target_assortcoef == 0:
        eta = np.outer(dist["q"], dist["q"])
        eta = pd.DataFrame(eta)
        eta.index = dist["d"]
        eta.columns = dist["d"]
        return eta, None, None
    if target_assortcoef:
        constrs.append(r == target_assortcoef)
        problem = cp.Problem(cp.Minimize(eta_obj(eta)), constrs)
        problem.solve(**kwargs)
        if problem.status != "optimal":
            # show a warning about solver status
            warnings.warn(f"Solver status: {problem.status}")
        eta = pd.DataFrame(eta.value)
        eta.index = dist["d"]
        eta.columns = dist["d"]
        return eta, problem.status, problem.solver_stats
    else:
        problem1 = cp.Problem(cp.Minimize(r), constrs)
        problem1.solve(**kwargs)
        if problem1.status != "optimal":
            # show a warning about solver status
            warnings.warn(f"Lower bound solver status: {problem1.status}")
        problem2 = cp.Problem(cp.Maximize(r), constrs)
        problem2.solve(**kwargs)
        if problem2.status != "optimal":
            # show a warning about solver status
            warnings.warn(f"Upper bound solver status: {problem2.status}")
        return (
            [problem1.value, problem2.value],
            [problem1.status, problem1.solver_stats],
            [problem2.status, problem2.solver_stats],
        )


# Degree preserving rewiring for directed networks
def dprewire_directed(
    netwk,
    iteration,
    nattempts,
    history,
    eta,
):
    """
    Rewire an unweighted directed network towards target assortativity
    coefficient(s) while preserving node degrees.

    Parameters
    ----------
    netwk (WDNet): An instance of the WDNet class.

    iteration (int): Number of iterations to run the rewiring algorithm.
        Each iteration consists of nattempts rewiring attempts.

    nattempts (int): Number of rewiring attempts per iteration.

    history (bool): If True, the rewiring history is returned.

    eta (DataFrame): A pandas DataFrame representing the proportion of
        edges emerging from nodes with out-degree i, in-degree j, and
        pointing to nodes with out-degree k, in-degree l.

    Returns
    -------
    The rewired network (WDNet), the assortativity coefficients after
    each iteration, and the rewiring history (if history is True).

    history (Optional[DataFrame]): A pandas DataFrame representing the
        rewiring history.
    """
    snodes, tnodes = np.array(netwk.edgelist).T
    outd = netwk.node_attr["outs"].values.astype(int)
    ind = netwk.node_attr["ins"].values.astype(int)
    sout = outd[snodes]
    sin = ind[snodes]
    tout = outd[tnodes]
    tin = ind[tnodes]
    dfs = pd.DataFrame({"type": eta.index.values, "index": np.arange(eta.shape[0])})
    dfs.set_index("type", inplace=True)
    dft = pd.DataFrame({"type": eta.columns.values, "index": np.arange(eta.shape[1])})
    dft.set_index("type", inplace=True)

    sindex = [dfs.loc[f"{i}-{j}", "index"] for i, j in zip(sout, sin)]
    tindex = [dft.loc[f"{i}-{j}", "index"] for i, j in zip(tout, tin)]
    del dfs, dft
    tnode, outout, outin, inout, inin, rewire_history = rewire_directed(
        iteration=iteration,
        nattempts=nattempts,
        tnode=tnodes,
        sout=sout,
        sin=sin,
        tout=tout,
        tin=tin,
        index_s=sindex,
        index_t=tindex,
        eta=eta.values,
        history=history,
    )
    assortcoef = pd.DataFrame(
        {
            "outout": outout,
            "outin": outin,
            "inout": inout,
            "inin": inin,
        }
    )
    # insert one row at the beginning of assortcoef and reset the index
    initial_assort = pd.DataFrame([netwk.assortcoef()], columns=assortcoef.columns)
    assortcoef = pd.concat([initial_assort, assortcoef]).reset_index(drop=True)
    snodes = snodes.astype(int)
    tnodes = tnode.astype(int)
    netwk = WDNet(edgelist=np.array([snodes, tnodes]).T, directed=True)
    if not history:
        rewire_history = None
    return netwk, assortcoef, rewire_history


# Degree preserving rewiring for undirected networks
def dprewire_undirected(
    netwk,
    iteration,
    nattempts,
    history,
    eta,
):
    return None


# Main function for degree preserving rewiring
def dprewire(
    netwk: WDNet,
    target_assortcoef: Union[Dict[str, float], float],
    iteration: int = 200,
    nattempts: Optional[int] = None,
    history: bool = False,
    eta_obj: Optional[Callable] = lambda eta: 0,
    cvxpy_solver: str = "ECOS",
    eta: Optional[DataFrame] = None,
    **kwargs,
):
    """
    Rewire an unweighted directed or undirected network towards target
    assortativity coefficient(s) while preserving node degrees.

    Parameters
    ----------
    netwk (WDNet): An instance of the WDNet class.
    target_assortcoef (Optional[Union[Dict[str, float], float]]): Target assortativity coefficient(s).
        For directed networks, this should be a dictionary with keys
        'outout', 'outin', 'inout', and 'inin' corresponding to the
        four types of assortativity coefficients. For undirected
        networks, this should be a single float value.
    iteration (int): Number of iterations to run the rewiring algorithm.
        Each iteration consists of nattempts rewiring attempts.
    nattempts (Optional[int]): Number of rewiring attempts per iteration. If None,
        nattempts is set to the number of edges in the network.
    history (bool): If True, the rewiring history is returned.
    eta_obj (Optional[Callable]): A objective function that takes eta
        (a numpy array) as input and returns a scalar value. It will be
        minimized when solving for an appropriate eta.
        Defaults to 0 (i.e., no objective function).
    cvxpy_solver (str): The solver to use for solving the optimization
        problem. Defaults to 'ECOS'.
    eta (Optional[DataFrame]): A pandas DataFrame or None. By default,
        eta is computed with the given network, the target
        assortativity coefficient(s), and the given eta_obj function.
        If eta is provided, the target assortativity coefficient(s)
        and eta_obj function are ignored. For directed networks, the
        element eta.loc['i-j', 'k-l'] representing the proportion of
        edges emerging from nodes with out-degree i, in-degree j, and
        pointing to nodes with out-degree k, in-degree l. For
        undirected networks, the element eta.loc['i', 'j']
        representing the proportion of edges between nodes with degree
        i and j.
    **kwargs: Additional keyword arguments to be passed to cvxpy when solving
        the optimization problem.
        See https://www.cvxpy.org/tutorial/advanced/index.html#setting-solver-options
        for details.
    """
    if netwk.weighted:
        raise ValueError("The network must be unweighted.")
    check_target_assortcoef(netwk, target_assortcoef)
    if nattempts is None:
        nattempts = netwk.edge_attr.shape[0]
    if netwk.directed:
        if eta is None:
            eta, _, _ = get_eta_directed(
                netwk=netwk,
                target_assortcoef=target_assortcoef,
                eta_obj=eta_obj,
                which_range=None,
                solver=cvxpy_solver,
                **kwargs,
            )  # type: ignore
        return dprewire_directed(
            netwk=netwk,
            iteration=iteration,
            nattempts=nattempts,
            history=history,
            eta=eta,
        )
    else:
        if eta is None:
            eta, _, _ = get_eta_undirected(
                netwk=netwk,
                target_assortcoef=target_assortcoef,
                eta_obj=eta_obj,
                solver=cvxpy_solver,
                **kwargs,
            )  # type: ignore
        return dprewire_undirected(
            netwk=netwk,
            iteration=iteration,
            nattempts=nattempts,
            history=history,
            eta=eta,
        )


# Function to compute the range of assortativity coefficients
def dprewire_range(
    netwk: WDNet,
    which_range: Optional[str] = None,
    target_assortcoef: Optional[Union[Dict[str, float], float]] = None,
    cvxpy_solver: str = "ECOS",
    **kwargs,
):
    """
    Compute the range of assortativity coefficient for an unweighted
    directed or undirected network.

    Parameters
    ----------
    netwk (WDNet): An instance of the WDNet class.

    which_range (Optional[str]): Which range of assortativity
        coefficients to compute. For directed networks, this should be
        one of 'outout', 'outin', 'inout', and 'inin'. For undirected
        networks, this should be None.

    target_assortcoef (Optional[Union[Dict[str, float], float]]): Target
        assortativity coefficient(s). For directed networks, this should be a
        dictionary with keys 'outout', 'outin', 'inout', and 'inin'
        corresponding to the four types of assortativity
        coefficients. For undirected networks, this should be a None.

    cvxpy_solver (str): The solver to use for solving the optimization
        problem. Defaults to 'ECOS'.

    **kwargs: Additional keyword arguments to be passed to cvxpy when solving
        the optimization problem.
        See https://www.cvxpy.org/tutorial/advanced/index.html#setting-solver-options
        for details.
    """
    if which_range is None and netwk.directed:
        raise ValueError("which_range must be specified for directed networks.")
    if netwk.weighted:
        raise ValueError("The network must be unweighted.")
    check_target_assortcoef(netwk, target_assortcoef)
    if netwk.directed:
        ret, _, _ = get_eta_directed(
            netwk=netwk,
            target_assortcoef=target_assortcoef,
            eta_obj=lambda eta: 0,
            which_range=which_range,
            solver=cvxpy_solver,
            **kwargs,
        )
        return ret
    else:
        if which_range is not None:
            warnings.warn('"which_range" is ignored for undirected networks.')
        ret, _, _ = get_eta_undirected(
            netwk=netwk,
            target_assortcoef=target_assortcoef,
            eta_obj=lambda eta: 0,
            solver=cvxpy_solver,
            **kwargs,
        )
        return ret
