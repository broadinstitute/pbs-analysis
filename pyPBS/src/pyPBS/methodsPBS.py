import numpy as np
import pandas as pd
from scipy.stats import gamma
from scipy.interpolate import interp1d

def get_probability_being_signal(bin_df, params_df, theta, max_dens=20):
    """
    Calculate the probability of being a signal for each bin.

    Parameters
    ----------
    bin_df : pandas.DataFrame
        A DataFrame containing genomic bin data with columns `chr`, `start`, `end`, and `counts`.
    params_df : dict
        A dictionary containing the parameters `beta` and `k` for the gamma distribution.
    theta : float
        The quantile threshold for determining the noise level.
    max_dens : float, optional
        The maximum density value for the x-axis, by default 20.

    Returns
    -------
    dict
        A dictionary containing:
        - `bin_df` (pandas.DataFrame): The input DataFrame with an additional `pbs` column.
        - `empirical_density_df` (pandas.DataFrame): The empirical density data.
        - `fitted_density_df` (pandas.DataFrame): The fitted gamma density data.
        - `beta` (float): The beta parameter of the gamma distribution.
        - `k` (float): The k (shape) parameter of the gamma distribution.
        - `lambda` (float): The lambda parameter for scaling the fitted density.
        - `RH_area` (float or None): The right-hand area under the curve, if available in `params_df`.
    """
    # Filter bins with counts > 0
    working_df = bin_df[bin_df['counts'] > 0]

    # Calculate threshold for counts
    threshold = np.quantile(working_df.loc[working_df['chr'].str.startswith("chr"), 'counts'], 0.9999)

    # Estimate empirical density
    counts_below_threshold = working_df.loc[working_df['counts'] < threshold, 'counts']
    APP_density_x = np.linspace(min(counts_below_threshold), max(counts_below_threshold), 10000)
    APP_density_y = np.histogram(counts_below_threshold, bins=10000, density=True)[0]
    empirical_density_df = pd.DataFrame({'x': APP_density_x, 'y': APP_density_y})

    # Extract beta, k from params_df
    beta = params_df['beta']
    k = params_df['k']

    # Calculate p-values using the gamma distribution
    p_values = gamma.sf(working_df['counts'], a=k, scale=1/beta)  # sf = 1 - cdf
    ratio = np.mean(p_values > (1 - theta))
    lambda_ = ratio / theta

    print("########## Assigning PBS score to each bin #################")

    # Calculate fitted density for background
    fitted_density_y = gamma.pdf(APP_density_x, a=k, scale=1/beta)
    fitted_density_df = pd.DataFrame({'x': APP_density_x, 'y': fitted_density_y})

    # Calculate PBS scores
    empirical_density_df['pbs'] = (empirical_density_df['y'] - lambda_ * fitted_density_df['y']) / empirical_density_df['y']

    # Set PBS values to 0 where unstable
    first_zero_idx = np.max(np.where((empirical_density_df['pbs'] <= 0) & np.isfinite(empirical_density_df['pbs'])))
    empirical_density_df.loc[:first_zero_idx, 'pbs'] = 0

    # Set PBS values to 1 where they are greater than first_zero_idx and not finite
    empirical_density_df.loc[(empirical_density_df['pbs'] > first_zero_idx) & ~np.isfinite(empirical_density_df['pbs']), 'pbs'] = 1

    # Add a point to represent the highest value of counts
    max_count = working_df['counts'].max()
    empirical_density_df = pd.concat([empirical_density_df, pd.DataFrame({'x': [max_count], 'y': [0], 'pbs': [1]})])

    # Approximate PBS values for bin_df
    interp_func = interp1d(empirical_density_df['x'], empirical_density_df['pbs'], bounds_error=False, fill_value=0)
    bin_df['pbs'] = interp_func(bin_df['counts'])

    return {
        "bin_df": bin_df,
        "empirical_density_df": empirical_density_df,
        "fitted_density_df": fitted_density_df,
        "beta": beta,
        "k": k,
        "lambda": lambda_,
        "RH_area": params_df.get('RH_area', None)
    }


def get_probability_being_signal_fixed_lambda(bin_df, params_df, theta, max_dens=20):
    """
    Calculate the probability of being a signal for each bin with a fixed lambda.

    Parameters
    ----------
    bin_df : pandas.DataFrame
        A DataFrame containing genomic bin data with columns `chr`, `start`, `end`, and `counts`.
    params_df : dict
        A dictionary containing the parameters `beta`, `k`, and `lambda` for the gamma distribution.
    theta : float
        The quantile threshold for determining the noise level.
    max_dens : float, optional
        The maximum density value for the x-axis, by default 20.

    Returns
    -------
    dict
        A dictionary containing:
        - `bin_df` (pandas.DataFrame): The input DataFrame with an additional `pbs` column.
        - `empirical_density_df` (pandas.DataFrame): The empirical density data.
        - `fitted_density_df` (pandas.DataFrame): The fitted gamma density data.
        - `beta` (float): The beta parameter of the gamma distribution.
        - `k` (float): The k (shape) parameter of the gamma distribution.
        - `lambda` (float): The lambda parameter for scaling the fitted density.
        - `RH_area` (float or None): The right-hand area under the curve, if available in `params_df`.
    """
    # Filter bins with counts > 0
    working_df = bin_df[bin_df['counts'] > 0]

    # Calculate threshold for counts
    threshold = np.quantile(working_df.loc[working_df['chr'].str.startswith("chr"), 'counts'], 0.9999)

    # Estimate empirical density
    counts_below_threshold = working_df.loc[working_df['counts'] < threshold, 'counts']
    APP_density_x = np.linspace(min(counts_below_threshold), max(counts_below_threshold), 10000)
    APP_density_y = np.histogram(counts_below_threshold, bins=10000, density=True)[0]
    empirical_density_df = pd.DataFrame({'x': APP_density_x, 'y': APP_density_y})

    # Extract beta, k, and lambda from params_df
    beta = params_df['beta']
    k = params_df['k']
    lambda_ = params_df['lambda']

    print("########## Assigning PBS score to each bin #################")

    # Calculate fitted density for background
    fitted_density_y = gamma.pdf(APP_density_x, a=k, scale=1/beta)
    fitted_density_df = pd.DataFrame({'x': APP_density_x, 'y': fitted_density_y})

    # Calculate PBS scores
    empirical_density_df['pbs'] = (empirical_density_df['y'] - lambda_ * fitted_density_df['y']) / empirical_density_df['y']

    # Set PBS values to 0 where unstable
    first_zero_idx = np.max(np.where((empirical_density_df['pbs'] <= 0) & np.isfinite(empirical_density_df['pbs'])))
    empirical_density_df.loc[:first_zero_idx, 'pbs'] = 0

    # Set PBS values to 1 where they are greater than first_zero_idx and not finite
    empirical_density_df.loc[(empirical_density_df['pbs'] > first_zero_idx) & ~np.isfinite(empirical_density_df['pbs']), 'pbs'] = 1

    # Add a point to represent the highest value of counts
    max_count = working_df['counts'].max()
    empirical_density_df = pd.concat([empirical_density_df, pd.DataFrame({'x': [max_count], 'y': [0], 'pbs': [1]})])

    # Approximate PBS values for bin_df
    interp_func = interp1d(empirical_density_df['x'], empirical_density_df['pbs'], bounds_error=False, fill_value=0)
    bin_df['pbs'] = interp_func(bin_df['counts'])

    return {
        "bin_df": bin_df,
        "empirical_density_df": empirical_density_df,
        "fitted_density_df": fitted_density_df,
        "beta": beta,
        "k": k,
        "lambda": lambda_,
        "RH_area": params_df.get('RH_area', None)
    }