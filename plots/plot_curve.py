import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    """

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


script_dir = os.path.dirname(__file__)
data = pd.read_csv(os.path.join(script_dir, "..", "data", "solution.dat"))

tangents = [data[f"t{i}"].to_numpy() for i in range(1, 4)]

curve = np.array([data[f"p{i}"].to_list() for i in range(1, 4)])

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.plot(curve[0], curve[1], curve[2])
ax.set_box_aspect([1, 1, 1])
set_axes_equal(ax)
plt.show()
