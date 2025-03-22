import pandas as pd
import argparse
from .methodsPBS import get_probability_being_signal
from .fitDistribution import get_distribution_parameters_with_optim
from .plotPBS import plot_multiplot_histograms, plot_fitted_distribution
from .utilsPBS import (
    get_pbs,
    get_corrected_pbs,
    get_pbs_plot,
    get_corrected_pbs_plot,
    get_compartment_fits,
    write_pbs,
    write_corrected_pbs,
)

def PBS(bin_df, theta=0.5, plot_pbs=True, fit_regions_df=None):
    """
    Main function to calculate PBS (Probability of Being Signal) for genomic bins.

    Parameters:
        bin_df (str or pd.DataFrame): Path to a file or a DataFrame containing genomic bin data.
        theta (float): Estimate of noise (default 0.5).
        plot_pbs (bool): Whether to plot PBS results (default True).
        fit_regions_df (str or pd.DataFrame): Path to a file or a DataFrame containing regions to fit.

    Returns:
        dict: A dictionary containing the PBS object and corrected PBS object.
    """
    # Load bin_df if it's a file path
    if isinstance(bin_df, str):
        bin_df = pd.read_csv(
            bin_df, sep="\t", header=None, names=["chr", "start", "end", "counts"]
        )
    else:
        bin_df.columns = ["chr", "start", "end", "counts"]

    # Get initial distribution parameters and calculate PBS
    params_df = get_distribution_parameters_with_optim(working_df=bin_df, theta=theta)
    pbs_obj = get_probability_being_signal(bin_df=bin_df, params_df=params_df, theta=theta)

    corr_pbs_obj = None

    # If fit_regions_df is provided, calculate corrected PBS
    if fit_regions_df is not None:
        if isinstance(fit_regions_df, str):
            fit_regions_df = pd.read_csv(
                fit_regions_df, sep="\t", header=None, names=["chr", "start", "end", "counts"]
            )
        else:
            fit_regions_df.columns = ["chr", "start", "end", "counts"]

        fit_regions_params_df = get_distribution_parameters_with_optim(
            working_df=fit_regions_df, theta=theta
        )
        corr_pbs_obj = get_probability_being_signal(
            bin_df=bin_df, params_df=fit_regions_params_df, theta=theta
        )
        corr_pbs_obj["fit_regions_df"] = fit_regions_df
        corr_pbs_obj["fit_regions_lambda"] = fit_regions_params_df["lambda"]

        if plot_pbs:
            plot_multiplot_histograms(corr_pbs_obj)
    else:
        if plot_pbs:
            plot_multiplot_histograms(pbs_obj)

    return {"pbs_obj": pbs_obj, "corr_pbs_obj": corr_pbs_obj}


# Example usage
if __name__ == "__main__":
    # Example input data
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run PBS analysis on genomic bins.")
    parser.add_argument(
        "--bin_df_path", type=str, required=True, help="Path to the bin_df file."
    )
    parser.add_argument(
        "--fit_regions_df_path",
        type=str,
        required=False,
        help="Path to the fit_regions_df file (optional).",
    )
    parser.add_argument(
        "--theta", type=float, default=0.5, help="Estimate of noise (default 0.5)."
    )
    parser.add_argument(
        "--plot_pbs",
        action="store_true",
        help="Whether to plot PBS results (default False).",
    )
    args = parser.parse_args()

    # Assign arguments to variables
    bin_df_path = args.bin_df_path
    fit_regions_df_path = args.fit_regions_df_path
    theta = args.theta
    plot_pbs = args.plot_pbs

    # Run PBS function
    results = PBS(bin_df=bin_df_path, theta=0.5, plot_pbs=True, fit_regions_df=fit_regions_df_path)

    # Access results
    pbs_obj = results["pbs_obj"]
    corr_pbs_obj = results["corr_pbs_obj"]