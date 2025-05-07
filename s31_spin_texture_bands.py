#!/usr/bin/env python

import click
import numpy as np

@click.command()
@click.argument('density', type=int)
@click.argument('x_a', type=float)
@click.argument('x_b', type=float)
@click.argument('x_c', type=float)
@click.argument('y_a', type=float)
@click.argument('y_b', type=float)
@click.argument('y_c', type=float)
def main(density, x_a, x_b, x_c, y_a, y_b, y_c):
    x = np.array([x_a, x_b, x_c])
    y = np.array([y_a, y_b, y_c])

    # Create linear steps between -y and y vectors
    y_steps = np.linspace(-1, 1, density)
    
    for t in y_steps:
        y_step = t * y
        first_point = -x + y_step
        second_point = x + y_step
        print(f"output band_mulliken "
        f"{first_point[0]:.8f} {first_point[1]:.8f} {first_point[2]:.8f} "
        f"{second_point[0]:.8f} {second_point[1]:.8f} {second_point[2]:.8f} "
        f"{density} X Y")

if __name__ == '__main__':
    main()
