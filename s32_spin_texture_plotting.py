#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import click

@click.command()
@click.argument("state", type=int)
@click.option("--uniform-arrows", is_flag=True, help="If set, scale all arrows to have the same length.")
@click.option("--scale-factor", default=50.0, type=float, help="Scaling divisor for arrow length.")
@click.option("--no-border-arrows", is_flag=True, help="If set, removes arrows that start on the plot borders.")
def main(state, uniform_arrows, scale_factor, no_border_arrows): 
	# Load and filter data
    data = np.loadtxt("spin_texture.dat", comments='#')
    data = data[data[:, 4] == state]

    # Extract columns
    x = data[:, 1]
    y = data[:, 2]
    u = data[:, 6]
    v = data[:, 7]
    z = data[:, 8]  # Scalar field for background

    # Interpolate z-values onto a grid
    grid_x, grid_y = np.mgrid[x.min():x.max():200j, y.min():y.max():200j]
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')

    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    if no_border_arrows:
        # Define a small epsilon for numerical safety
        eps = 1e-8
        mask = (
            (x > x.min() + eps) & (x < x.max() - eps) &
            (y > y.min() + eps) & (y < y.max() - eps)
        )
        x = x[mask]
        y = y[mask]
        u = u[mask]
        v = v[mask]

    # Plot
    plt.figure(figsize=(8, 6))
    plt.contourf(grid_x, grid_y, grid_z, levels=100, cmap='coolwarm', vmin=-1, vmax=1)  # blue to red
    plt.colorbar(label=r'$\sigma_{z}$ magnitude')
    if uniform_arrows:
        # Normalize each vector to unit length
        magnitudes = np.sqrt(u**2 + v**2)
        u = u / magnitudes
        v = v / magnitudes

    plt.quiver(x, y, u / scale_factor, v / scale_factor,
               angles='xy', scale_units='xy', scale=1, color='black')	

    # Custom axis ticks and labels
    xticks = [x_min, (x_min + x_max) / 2, x_max]
    yticks = [y_min, (y_min + y_max) / 2, y_max]
    plt.xticks(xticks, ['-X', r'$\Gamma$', 'X'])
    plt.yticks(yticks, ['-Y', r'$\Gamma$', 'Y'])

    # Clean axis styling
    plt.xlabel("")
    plt.ylabel("")
#    plt.axis('equal')
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.savefig(f"{state}_spin_texture.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
