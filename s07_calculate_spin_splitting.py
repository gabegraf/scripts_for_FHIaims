#!/usr/bin/env python
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, CubicSpline

def calculate_splitting(file_path, x_range, state):
    
    # rc_081724
    # file_path takes in the address of the csv file (e.g., path/to/file/data.csv). The csv file must have three columns labelled as "State", "x-value", and "y-value"
    # x-range is the range of x-values to perform the fitting on (e.g., [0.3, 0.5]); the range should cover the lowest part of the band and should be refined iteratively to obtain a good fitting
    # the state is the index of the CBM or VBM
    
    
    ################################################################################################################
    
    plt.rcParams['axes.linewidth'] = 2
    
    # Read the csv file and convert it to a dataframe
    df=pd.read_csv(file_path)
    
    print("Data Columns:", df.columns)
    print("First few rows of data:")
    print(df.head())
    
    # Print the x_range and state for debugging
    print(f"x_range: {x_range}")
    print(f"state: {state}")
    
    # Filter the data for the specified state and x-range
    subset = df[(df['State'] == state) & (df['x-value'] >= x_range[0]) & (df['x-value'] <= x_range[1])]
    
    # Print the subset for debugging
    print("Filtered Data Subset:")
    print(subset)
    
    if subset.empty:
        print(f"No data found for state {state} with x-values in range {x_range}")
        return  # or handle the case as needed

    # Filter the data for the specified state and x-range
    subset = df[(df['State'] == state) & (df['x-value'] >= x_range[0]) & (df['x-value'] <= x_range[1])]
    
    spline = CubicSpline(subset['x-value'],subset['y-value'])
    x_fit = np.linspace(subset['x-value'].min(), subset['x-value'].max(), 200)
    y_fit = spline(x_fit)
	
    spline_derivative = spline.derivative()
    spline_critical_points = spline_derivative.roots()
    spline_second_derivative = spline.derivative(2)
    min_values = []
    for cp in spline_critical_points:
       if spline_second_derivative(cp) > 0:
          min_values.append((cp, spline(cp)))
    if len(min_values) == 1:
        vertex_x = min_values[0][0]
        vertex_y = min_values[0][1]
    else:
        print("Multiple minima detected in cubic spline")
        exit(1)

    # Δk calculation
    delta_k = vertex_x - df['x-value'].min()
    if delta_k < 0:
        delta_k = df['x-value'].max() - vertex_x        # in case the spin-splitting happens in the other direction
    
    # Plot the complete data for the state and the fitted curve
    plt.figure(figsize=(10, 6))
    plt.scatter(df[df['State'] == state]['x-value'], df[df['State'] == state]['y-value'], color='blue', label=f'State {state}')
    plt.plot(df[df['State'] == state]['x-value'], df[df['State'] == state]['y-value'], 'b-', alpha=0.3)
    plt.plot(x_fit, y_fit, 'r-', label=f'Spline Fit {state}')
    
    # Handle the other state
    other_state = df[df['State'] != state]['State'].unique()[0]
    other_data = df[df['State'] == other_state]
    plt.scatter(other_data['x-value'], other_data['y-value'], color='green', label=f'State {other_state}')
    plt.plot(other_data['x-value'], other_data['y-value'], 'g-', alpha=0.3)
    
    # Interpolation for the other state at vertex_x
    interpolator = CubicSpline(other_data['x-value'], other_data['y-value'])
    # interpolated_band = interpolator(x_fit)
    plt.plot(x_fit, interpolator(x_fit), 'k-', alpha=1, label='interpolated band')
    interpolated_y = interpolator(vertex_x)
    
    # ΔE calculation
    delta_e = np.abs(interpolated_y - vertex_y)     # absolute value for displaying the valence band splitting value appropiately
    
    # Annotate Δk and ΔE with double-headed arrows
    plt.annotate('', xy=(df['x-value'].min(), vertex_y), xytext=(vertex_x, vertex_y),
                 arrowprops=dict(arrowstyle="<->", lw=1.5, color='purple'))
    plt.annotate('', xy=(vertex_x, vertex_y), xytext=(vertex_x, interpolated_y),
                 arrowprops=dict(arrowstyle="<->", lw=1.5, color='black'))
    
    plt.xlabel('k', fontname="Arial", fontsize=16, fontweight='bold')
    plt.ylabel('Energy (eV)', fontname="Arial", fontsize=16, fontweight='bold')
    plt.xticks(size=14, fontname="Arial", fontweight='bold')
    plt.yticks(size=14, fontname="Arial", fontweight='bold')
    plt.legend()
    
    textstr = f'$\Delta k = {delta_k:.2f}$ \u212B\u207B\u00B9\n$\Delta E = {delta_e*1000:.2f}$ meV'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.gca().text(0.75, 0.15, textstr, transform=plt.gca().transAxes, fontsize=14,
                   verticalalignment='bottom', bbox=props)
    
    directory = os.path.dirname(file_path)
    plt.savefig(os.path.join(directory, f'spin_splitting_{state}.png'), dpi= 300)

    plt.show()

    
# Define variables before running
file_path = "./spin_splitting_input_for_bandmlk1002.csv"  # Path to CSV file
x_range = [0.37, 0.85]
state = 985

calculate_splitting(file_path, x_range, state)

