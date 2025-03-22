import pandas as pd
from .plotPBS import plot_fitted_distribution, plot_multiplot_histograms


def get_pbs(pbs_obj):
    """
    Get the bin_df from the PBS object.

    Parameters:
        pbs_obj (dict): The PBS object.

    Returns:
        pd.DataFrame: The bin_df from the PBS object.
    """
    return pbs_obj["pbs_obj"]["bin_df"]


def get_corrected_pbs(pbs_obj):
    """
    Get the bin_df from the corrected PBS object.

    Parameters:
        pbs_obj (dict): The PBS object.

    Returns:
        pd.DataFrame: The bin_df from the corrected PBS object.
    """
    return pbs_obj["corr_pbs_obj"]["bin_df"]


def get_pbs_plot(pbs_obj, max_dens=20, title_str="PBS"):
    """
    Plot the PBS data using the fitted distribution.

    Parameters:
        pbs_obj (dict): The PBS object.
        max_dens (float): Maximum density for the x-axis.
        title_str (str): Title for the plot.

    Returns:
        None
    """
    plot_obj = pbs_obj.get("pbs_obj", pbs_obj)
    plot_fitted_distribution(
        working_df=plot_obj["bin_df"],
        empirical_density_df=plot_obj["empirical_density_df"],
        fitted_density_df=plot_obj["fitted_density_df"],
        max_dens=max_dens,
        title_str=title_str,
        beta=plot_obj["beta"],
        k=plot_obj["k"],
        lambda_=plot_obj["lambda"],
    )


def get_corrected_pbs_plot(pbs_obj, max_dens=20, title_str="PBS"):
    """
    Plot the corrected PBS data using the fitted distribution.

    Parameters:
        pbs_obj (dict): The PBS object.
        max_dens (float): Maximum density for the x-axis.
        title_str (str): Title for the plot.

    Returns:
        None
    """
    plot_obj = pbs_obj.get("corr_pbs_obj", pbs_obj)
    plot_fitted_distribution(
        working_df=plot_obj["bin_df"],
        empirical_density_df=plot_obj["empirical_density_df"],
        fitted_density_df=plot_obj["fitted_density_df"],
        max_dens=max_dens,
        title_str=title_str,
        beta=plot_obj["beta"],
        k=plot_obj["k"],
        lambda_=plot_obj["lambda"],
    )


def get_compartment_fits(pbs_obj):
    """
    Plot multiple histograms for the PBS object.

    Parameters:
        pbs_obj (dict): The PBS object.

    Returns:
        None
    """
    plot_multiplot_histograms(pbs_obj)


def write_pbs(pbs_obj, output_filename):
    """
    Write the PBS bin_df to a file.

    Parameters:
        pbs_obj (dict): The PBS object.
        output_filename (str): The output file path.

    Returns:
        None
    """
    bin_df = pbs_obj["pbs_obj"]["bin_df"]
    bin_df.to_csv(output_filename, sep="\t", index=False, header=True)


def write_corrected_pbs(pbs_obj, output_filename):
    """
    Write the corrected PBS bin_df to a file.

    Parameters:
        pbs_obj (dict): The PBS object.
        output_filename (str): The output file path.

    Returns:
        None
    """
    bin_df = pbs_obj["corr_pbs_obj"]["bin_df"]
    bin_df.to_csv(output_filename, sep="\t", index=False, header=True)