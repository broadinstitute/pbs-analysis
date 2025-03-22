import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gamma

def plot_multiplot_histograms(pbs_obj, bin_width=0.05, max_dens=20):
    """
    Plot multiple histograms with fitted gamma distributions.

    Parameters
    ----------
    pbs_obj : dict
        A dictionary containing PBS data, including `bin_df`, `beta`, `k`, and `lambda` for both `pbs_obj` and `corr_pbs_obj`.
    bin_width : float, optional
        The width of the bins in the histogram, by default 0.05.
    max_dens : float, optional
        The maximum density value for the x-axis, by default 20.

    Returns
    -------
    None
        Displays the plots.
    """
    # Simulated data for the first plot
    x = np.linspace(0, max_dens, 500)
    y1 = gamma.pdf(x, a=pbs_obj['pbs_obj']['k'], scale=1/pbs_obj['pbs_obj']['beta'])
    sim_df1 = pd.DataFrame({'x': x, 'y': y1})

    # First plot
    bin_df = pbs_obj['pbs_obj']['bin_df']
    plt.figure(figsize=(10, 8))
    plt.subplot(2, 1, 1)
    sns.histplot(bin_df['counts'], binwidth=bin_width, stat="density", alpha=0.5, color="blue")
    plt.plot(sim_df1['x'], sim_df1['y'] * pbs_obj['pbs_obj']['lambda'], color="red")
    plt.xlim(0, max_dens)
    plt.title(f"Bin Counts Histogram - All\nbeta, k, lambda: {round(pbs_obj['pbs_obj']['beta'], 2)}, "
              f"{round(pbs_obj['pbs_obj']['k'], 2)}, {round(pbs_obj['pbs_obj']['lambda'], 2)}")
    plt.xlabel("Counts")
    plt.ylabel("Density")

    # Simulated data for the second plot
    y2 = gamma.pdf(x, a=pbs_obj['corr_pbs_obj']['k'], scale=1/pbs_obj['corr_pbs_obj']['beta'])
    sim_df2 = pd.DataFrame({'x': x, 'y': y2})

    # Second plot
    fit_regions_df = pbs_obj['corr_pbs_obj']['fit_regions_df']
    plt.subplot(2, 1, 2)
    sns.histplot(fit_regions_df['counts'], binwidth=bin_width, stat="density", alpha=0.5, color="green")
    plt.plot(sim_df2['x'], sim_df2['y'] * pbs_obj['corr_pbs_obj']['fit_regions_lambda'], color="orange")
    plt.xlim(0, max_dens)
    plt.title(f"Bin Counts Histogram - Regions to Fit\nbeta, k, lambda: {round(pbs_obj['corr_pbs_obj']['beta'], 2)}, "
              f"{round(pbs_obj['corr_pbs_obj']['k'], 2)}, {round(pbs_obj['corr_pbs_obj']['fit_regions_lambda'], 2)}")
    plt.xlabel("Counts")
    plt.ylabel("Density")

    plt.tight_layout()
    plt.show()


def plot_cvm_heatmap(cvm_stat_df):
    """
    Plot a heatmap for CVM (Cramer-von Mises) statistics.

    Parameters
    ----------
    cvm_stat_df : pandas.DataFrame
        A DataFrame containing CVM statistics with columns `k`, `beta`, and `omega2`.

    Returns
    -------
    None
        Displays the heatmap.
    """
    pivot_table = cvm_stat_df.pivot("k", "beta", "omega2")
    plt.figure(figsize=(8, 6))
    sns.heatmap(pivot_table, cmap="viridis", annot=False)
    plt.title("CVM Heatmap")
    plt.xlabel("Beta")
    plt.ylabel("K")
    plt.show()


def plot_fitted_distribution(working_df, empirical_density_df, fitted_density_df, max_dens, title_str, beta, k, lambda_):
    """
    Plot a histogram with empirical and fitted gamma distributions.

    Parameters
    ----------
    working_df : pandas.DataFrame
        A DataFrame containing the counts data.
    empirical_density_df : pandas.DataFrame
        A DataFrame containing the empirical density data with columns `x` and `y`.
    fitted_density_df : pandas.DataFrame
        A DataFrame containing the fitted density data with columns `x` and `y`.
    max_dens : float
        The maximum density value for the x-axis.
    title_str : str
        The title of the plot.
    beta : float
        The beta parameter of the gamma distribution.
    k : float
        The k (shape) parameter of the gamma distribution.
    lambda_ : float
        The lambda parameter for scaling the fitted density.

    Returns
    -------
    None
        Displays the plot.
    """
    plt.figure(figsize=(10, 6))
    sns.histplot(working_df['counts'], binwidth=0.05, color="grey", stat="density", label="Histogram")
    plt.plot(empirical_density_df['x'], empirical_density_df['y'], color="lightgreen", label="Empirical", linewidth=2)
    plt.plot(fitted_density_df['x'], lambda_ * fitted_density_df['y'], color="salmon", label="Fitted", linewidth=2)
    plt.xlim(0, max_dens)
    plt.title(f"{title_str}\nbeta, k, lambda: {round(beta, 2)}, {round(k, 2)}, {round(lambda_, 2)}")
    plt.xlabel("Counts")
    plt.ylabel("Density")
    plt.legend()
    plt.show()


def plot_auc(working_df, params_df, xlim=20, theta=0.5):
    """
    Plot areas under the curve (AUC) between empirical and fitted distributions.

    Parameters
    ----------
    working_df : pandas.DataFrame
        A DataFrame containing the counts data.
    params_df : dict
        A dictionary containing the parameters `beta`, `k`, and `lambda`.
    xlim : float, optional
        The maximum value for the x-axis, by default 20.
    theta : float, optional
        The quantile threshold for splitting the AUC, by default 0.5.

    Returns
    -------
    None
        Displays the AUC plot.
    """
    # Empirical density
    density_pts = np.histogram(working_df['counts'], bins=1024, range=(0, xlim), density=True)
    x_values = np.linspace(0, xlim, len(density_pts[0]))
    y_empirical = density_pts[0]

    # Fitted density
    y_fitted = gamma.pdf(x_values, a=params_df['k'], scale=1/params_df['beta']) * params_df['lambda']

    # Calculate differences
    difference_df = pd.DataFrame({'xValues': x_values, 'yEmpirical': y_empirical, 'yFitted': y_fitted})
    lh_area = np.sum(np.abs(difference_df[difference_df['xValues'] < np.quantile(working_df['counts'], theta)]['yFitted'] -
                            difference_df[difference_df['xValues'] < np.quantile(working_df['counts'], theta)]['yEmpirical']))
    rh_area = np.sum(np.abs(difference_df[difference_df['xValues'] > np.quantile(working_df['counts'], theta)]['yFitted'] -
                            difference_df[difference_df['xValues'] > np.quantile(working_df['counts'], theta)]['yEmpirical']))

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(difference_df['xValues'], difference_df['yEmpirical'], label="Empirical", color="blue")
    plt.plot(difference_df['xValues'], difference_df['yFitted'], label="Fitted", color="red")
    plt.axvline(x=np.quantile(working_df['counts'], theta), linestyle="dotted", color="black", label="Theta")
    plt.title(f"AUC Plot\nLH Area: {lh_area:.2f}, RH Area: {rh_area:.2f}")
    plt.xlabel("Counts")
    plt.ylabel("Density")
    plt.legend()
    plt.show()