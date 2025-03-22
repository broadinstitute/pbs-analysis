import numpy as np
import pandas as pd
from scipy.stats import gamma, uniform
from scipy.optimize import minimize
from scipy.interpolate import interp1d


def get_distribution_parameters_with_optim(
    working_df,
    theta=0.5,
    lambda_range=np.arange(0.5, 1.1, 0.1),
    param_search_length=3,
    fix_weight=True,
    weight_value=500,
    max_dens=20,
):
    """
    Retrieve distribution parameters by fitting a gamma distribution to the data.

    Parameters
    ----------
    working_df : pandas.DataFrame
        A DataFrame containing genomic bin data with a column named `counts`.
    theta : float, optional
        The quantile threshold for determining the noise level, by default 0.5.
    lambda_range : numpy.ndarray, optional
        The range of lambda values to test, by default np.arange(0.5, 1.1, 0.1).
    param_search_length : int, optional
        The number of parameters to generate per parameter, by default 3.
    fix_weight : bool, optional
        Whether to fix the weight value, by default True.
    weight_value : float, optional
        The fixed weight value, by default 500.
    max_dens : float, optional
        The maximum density value for the x-axis, by default 20.

    Returns
    -------
    dict
        A dictionary containing the optimized parameters `beta`, `k`, `lambda`, and `RH_area`.
    """
    # Ensure the dataframe has a "counts" column
    counts_col = working_df.columns.str.contains("count", case=False)
    if counts_col.sum() != 1:
        raise ValueError('working_df must have exactly one column with a name similar to "counts".')
    working_df = working_df.rename(columns={working_df.columns[counts_col][0]: "counts"})
    working_df = working_df[working_df["counts"] > 0]

    # Get initial distribution parameters
    params_init = get_initial_distribution_parameters(
        working_df=working_df,
        lambda_range=lambda_range,
        param_search_length=param_search_length,
        theta=theta,
    )

    beta_init = params_init["beta"]
    k_init = params_init["k"]
    lambda_init = params_init["lambda"]

    if not fix_weight:
        weight_value = 0.99 * working_df["counts"].max()

    print("########## Fitting distribution #################")

    # Use optimization to calculate parameters
    def objective_function(par):
        params_df = {"beta": par[0], "k": par[1], "lambda": par[2]}
        return get_cvm_distance(
            params_df=params_df,
            weight_value=weight_value,
            working_df=working_df,
            theta=theta,
        )

    bounds = [(1e-5, None), (1e-5, None), (1e-5, 0.99)]
    result = minimize(
        objective_function,
        x0=[beta_init, k_init, lambda_init],
        method="L-BFGS-B",
        bounds=bounds,
    )

    # Calculate lambda based on the ratio of points below and above theta
    beta_opt, k_opt, lambda_opt = result.x
    p_values = gamma.sf(working_df["counts"], a=k_opt, scale=1 / beta_opt)
    ratio = np.mean(p_values > (1 - theta))
    lambda_final = ratio / theta

    params_df = {"beta": beta_opt, "k": k_opt, "lambda": lambda_final}

    # Get fit quality metric
    fit_quality = find_areas(
        working_df=working_df,
        params_df=params_df,
        xlim=max_dens,
        theta=theta,
    )
    params_df["RH_area"] = fit_quality["RH"]
    return params_df


def get_initial_distribution_parameters(
    working_df,
    theta=0.5,
    lambda_range=np.arange(0.4, 0.85, 0.05),
    fix_weight=True,
    weight_value=500,
    use_log=False,
    param_search_length=100,
    max_range_multiple=300,
    plot_data=False,
):
    """
    Generate initial values of lambda, beta, and k based on the minimum CvM statistic.

    Parameters
    ----------
    working_df : pandas.DataFrame
        A DataFrame containing genomic bin data with a column named `counts`.
    theta : float, optional
        The quantile threshold for determining the noise level, by default 0.5.
    lambda_range : numpy.ndarray, optional
        The range of lambda values to test, by default np.arange(0.4, 0.85, 0.05).
    fix_weight : bool, optional
        Whether to fix the weight value, by default True.
    weight_value : float, optional
        The fixed weight value, by default 500.
    use_log : bool, optional
        Whether to use the logarithm of counts, by default False.
    param_search_length : int, optional
        The number of parameters to generate per parameter, by default 100.
    max_range_multiple : int, optional
        The maximum range multiplier for beta and k, by default 300.
    plot_data : bool, optional
        Whether to plot the data, by default False.

    Returns
    -------
    pandas.Series
        A Series containing the initial parameters `beta`, `k`, and `lambda`.
    """
    working_df = working_df[working_df["counts"] > 0]

    if use_log:
        working_df["counts"] = np.log(working_df["counts"])
        bin_width = 0.1
        max_dens = 10
    else:
        bin_width = 10
        max_dens = 500

    mean_init = working_df["counts"].mean()
    var_init = working_df["counts"].var()
    beta_init = mean_init / var_init
    k_init = mean_init * beta_init

    params_all = get_parameter_combinations(
        beta_range=np.linspace(beta_init, max_range_multiple * beta_init, param_search_length),
        k_range=np.linspace(k_init, max_range_multiple * k_init, param_search_length),
        lambda_range=lambda_range,
    )

    print(f"########## Testing {len(params_all)} parameters for initial distribution values #################")

    if not fix_weight:
        weight_value = 0.99 * working_df["counts"].max()

    # Calculate CvM statistic for each parameter combination
    params_all["omega2"] = params_all.apply(
        lambda row: get_cvm_distance(
            params_df=row,
            working_df=working_df,
            weight_value=weight_value,
            theta=theta,
        ),
        axis=1,
    )

    # Retrieve parameters with the minimum CvM statistic
    params_init = params_all.loc[params_all["omega2"].idxmin()]
    return params_init


def get_parameter_combinations(beta_range, k_range, lambda_range):
    """
    Generate all combinations of beta, k, and lambda.

    Parameters
    ----------
    beta_range : numpy.ndarray
        The range of beta values.
    k_range : numpy.ndarray
        The range of k values.
    lambda_range : numpy.ndarray
        The range of lambda values.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing all combinations of beta, k, and lambda.
    """
    beta, k, lambda_ = np.meshgrid(beta_range, k_range, lambda_range, indexing="ij")
    params_all = pd.DataFrame(
        {
            "beta": beta.ravel(),
            "k": k.ravel(),
            "lambda": lambda_.ravel(),
        }
    )
    return params_all


def get_cvm_distance(working_df, params_df, weight_value, theta=0.5):
    """
    Calculate the Cramer-von Mises (CvM) statistic for goodness-of-fit.

    Parameters
    ----------
    working_df : pandas.DataFrame
        A DataFrame containing genomic bin data with a column named `counts`.
    params_df : dict
        A dictionary containing the parameters `beta`, `k`, and `lambda`.
    weight_value : float
        The fixed weight value.
    theta : float, optional
        The quantile threshold for determining the noise level, by default 0.5.

    Returns
    -------
    float
        The CvM statistic.
    """
    counts = working_df["counts"]
    pgamma_mix = pgamma_null(counts, params_df, weight_value)
    n = len(counts)
    U = np.sort(pgamma_mix)
    k = np.arange(1, n + 1)
    omega2 = (1 / (12 * n)) + np.sum((U - (2 * k - 1) / (2 * n)) ** 2)
    return omega2


def pgamma_null(q, params_df, weight_value):
    """
    Mixture of gamma and uniform distributions.

    Parameters
    ----------
    q : numpy.ndarray
        The input data.
    params_df : dict
        A dictionary containing the parameters `beta`, `k`, and `lambda`.
    weight_value : float
        The fixed weight value.

    Returns
    -------
    numpy.ndarray
        The mixture of gamma and uniform distributions.
    """
    pgamma_mix = (
        params_df["lambda"] * gamma.cdf(q, a=params_df["k"], scale=1 / params_df["beta"])
        + (1 - params_df["lambda"]) * uniform.cdf(q, loc=weight_value, scale=1.01 * weight_value - weight_value)
    )
    return pgamma_mix


def find_areas(working_df, params_df, xlim=20, theta=0.5):
    """
    Find areas under the curve between empirical and fitted distributions.

    Parameters
    ----------
    working_df : pandas.DataFrame
        A DataFrame containing genomic bin data with a column named `counts`.
    params_df : dict
        A dictionary containing the parameters `beta`, `k`, and `lambda`.
    xlim : float, optional
        The maximum value for the x-axis, by default 20.
    theta : float, optional
        The quantile threshold for splitting the AUC, by default 0.5.

    Returns
    -------
    dict
        A dictionary containing the left-hand (`LH`) and right-hand (`RH`) areas under the curve.
    """
    density_pts = np.histogram(working_df["counts"], bins=1024, range=(0, xlim), density=True)
    x_values = np.linspace(0, xlim, len(density_pts[0]))
    y_empirical = density_pts[0]
    y_fitted = gamma.pdf(x_values, a=params_df["k"], scale=1 / params_df["beta"]) * params_df["lambda"]

    difference_df = pd.DataFrame({"xValues": x_values, "yEmpirical": y_empirical, "yFitted": y_fitted})
    lh_area = np.sum(
        np.abs(
            difference_df.loc[difference_df["xValues"] < np.quantile(working_df["counts"], theta), "yFitted"]
            - difference_df.loc[difference_df["xValues"] < np.quantile(working_df["counts"], theta), "yEmpirical"]
        )
    )
    rh_area = np.sum(
        np.abs(
            difference_df.loc[difference_df["xValues"] > np.quantile(working_df["counts"], theta), "yFitted"]
            - difference_df.loc[difference_df["xValues"] > np.quantile(working_df["counts"], theta), "yEmpirical"]
        )
    )
    return {"LH": lh_area, "RH": rh_area}