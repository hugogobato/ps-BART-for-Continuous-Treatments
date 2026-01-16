import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data from Excel files
cover1 = pd.read_excel("BCF.xlsx")
cover2 = pd.read_excel("ps_BART.xlsx")

# Function to create and save a histogram
def create_histogram(data, bin_width, var_name, center=0.95):
    min_value = center - (center - data.min())  # Empirical min value
    max_value = 1.0
    bins = 20
    
    plt.figure(figsize=(8, 6))
    plt.hist(data, bins=bins, weights=[1/len(data)]*len(data), edgecolor='black')
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
    plt.ylim(0, 1)  # Set y-axis to range from 0 to 1
    plt.grid(axis='y')  # Only show horizontal grid lines
    plt.xlim(min_value, max_value)  # Ensure x-axis ranges from empirical min to max

    # Colors for the vertical lines
    colors = ['red', 'green','blue', 'purple']
    
    # Add vertical lines
    for i in range(1, 5):
        plt.axvline(x=center + 0.0125 * i, color=colors[i-1], linestyle='dotted', label=f'0.95 +/- {i}*0.0125')
        plt.axvline(x=center - 0.0125 * i, color=colors[i-1], linestyle='dotted')#, label=f'0.95 - {i}*0.0125')
    
    # Remove duplicate labels in the legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())

    #plt.title(f"Histogram of {var_name}")
    plt.xlabel(r'$\bar{C}$')
    plt.ylabel("Percentage")
    plt.savefig(f"{var_name}_histogram.png")
    plt.close()


# Create histograms for each variable
create_histogram(cover2["empirical_ATE_cover"], 0.05, "ps-BART_ATE")
create_histogram(cover1["empirical_ATE_cover"], 0.05, "BCF_ATE")
create_histogram(cover2["theoretical_CATE_cover"], 0.05, "ps-BART_CATE")
create_histogram(cover1["theoretical_CATE_cover"], 0.05, "BCF_CATE")

# Function to calculate the percentage of values within a given range
def calculate_percentage_within_range(data, center, increment, i):
    lower_bound = center - increment * i
    upper_bound = center + increment * i
    count_within_range = ((data >= lower_bound) & (data <= upper_bound)).sum()
    percentage_within_range = (count_within_range / len(data)) * 100
    return percentage_within_range

# Parameters
center = 0.95
increment = 0.0125

percentages_ps_bart_ATE = [calculate_percentage_within_range(cover2["empirical_ATE_cover"], center, increment, i) for i in range(1, 5)]

percentages_bcf_ATE = [calculate_percentage_within_range(cover1["empirical_ATE_cover"], center, increment, i) for i in range(1, 5)]

percentages_ps_bart_CATE = [calculate_percentage_within_range(cover2["theoretical_CATE_cover"], center, increment, i) for i in range(1, 5)]

percentages_bcf_CATE= [calculate_percentage_within_range(cover1["theoretical_CATE_cover"], center, increment, i) for i in range(1, 5)]

# Print results
print("Percentages for ps-BART ATE:", percentages_ps_bart_ATE)
print("Percentages for BCF CATE:", percentages_bcf_ATE)
print("Percentages for ps-BART CATE:", percentages_ps_bart_CATE)
print("Percentages for BCF CATE:", percentages_bcf_CATE)
