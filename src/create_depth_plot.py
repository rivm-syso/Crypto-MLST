import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def plot_coverage(input_file, output_pdf):
    """Generate a coverage plot from the input file and save it to a PDF."""
    colors = plt.cm.viridis(np.linspace(0, 1, 8))  # Create 8 different colors
    gene_lengths = [553, 503, 507, 526, 476, 587, 451, 558]

    # Read the coverage file
    try:
        coverage = pd.read_csv(input_file, sep="\t")
        coverage.columns = ["gene", "position", "depth"]
    except (FileNotFoundError, pd.errors.EmptyDataError):
        print(f"Error: The file {input_file} is empty or does not exist.")
        return

    # Create a PDF file and add the plot
    with PdfPages(output_pdf) as pdf:
        plt.figure(figsize=(10, 6))

        # Group by 'gene' and plot each group
        for (class_value, group), color, gene_length in zip(coverage.groupby('gene'), colors, gene_lengths):
            plt.plot(group['position'], group['depth'], label=f'{class_value}', color=color)
            plt.axvline(x=gene_length, color=color, linestyle='--', linewidth=1)  # Dashed line for gene length

        # Add chart metadata
        plt.xlabel('Position')
        plt.ylabel('Depth')
        plt.title(f'Coverage Plot for {input_file}')
        plt.legend()
        pdf.savefig()  # Save the figure to the PDF file
        plt.close()

    print(f"Coverage plot saved to {output_pdf}")

if __name__ == "__main__":
    # Set up argparse to handle command-line arguments
    parser = argparse.ArgumentParser(description="Generate a coverage plot from a coverage file.")
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input coverage file (TSV format with columns: gene, position, depth)."
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Path to the output PDF file where the plot will be saved."
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Pass the inputs to the function
    plot_coverage(args.input_file, args.output_file)
