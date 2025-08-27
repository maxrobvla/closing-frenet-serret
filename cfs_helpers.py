import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import build.closing_frenet_serret as cfs


def loss(args, number_field_periods, order_zeros, break_point):
    loss = cfs.loss(
        number_field_periods,
        order_zeros,
        np.array(args[:break_point]),
        np.array(args[break_point:]),
    )
    return loss


def save_curve(args, add_args):
    cfs.save_curve(
        "solution",
        add_args[0],  # number_field_periods
        add_args[1],  # order_zeros
        np.array(args[: add_args[2]]),
        np.array(args[add_args[2] :]),
    )


def close_curve(args, add_args, save=False):
    result = minimize(loss, args, args=add_args, method="Nelder-Mead")
    if save:
        save_curve(result.x, add_args)
    return result


def close_curve_with_jac(args, add_args, save=False):
    # working methods L-BFGS-B, TNC,
    result = minimize(loss, args, args=add_args, method="L-BFGS-B", jac=jac_loss)
    if save:
        save_curve(result.x, add_args)
    return result


def loss_with_jac(args, number_field_periods, order_zeros, break_point):
    loss = cfs.jacobian_and_scalar_loss(
        number_field_periods,
        order_zeros,
        np.array(args[:break_point]),
        np.array(args[break_point:]),
    )
    return loss


def vector_loss_with_jac(args, number_field_periods, order_zeros, break_point):
    loss = cfs.jacobian_and_vector_loss(
        number_field_periods,
        order_zeros,
        np.array(args[:break_point]),
        np.array(args[break_point:]),
    )
    return loss


def jac_loss(args, number_field_periods, order_zeros, break_point):
    result = cfs.loss_jacobian(
        number_field_periods,
        order_zeros,
        np.array(args[:break_point]),
        np.array(args[break_point:]),
    )
    return result


def vector_loss(args, number_field_periods, order_zeros, break_point):
    loss = cfs.loss_vector(
        number_field_periods,
        order_zeros,
        np.array(args[:break_point]),
        np.array(args[break_point:]),
    )
    return np.array(loss)


def jac_vector_loss(args, number_field_periods, order_zeros, break_point):
    result = cfs.loss_jacobian_vector(
        number_field_periods,
        order_zeros,
        np.array(args[:break_point]),
        np.array(args[break_point:]),
    )
    return np.array(result)


def close_curve_with_fixed_params(args, add_args, save=False):
    # add_args like: number_field_periods: int, order_zeros: array[int], break_point: int,
    result = minimize(loss, args, args=add_args, method="Nelder-Mead")
    if save:
        save_curve(result.x, add_args)
    return result


# plotting stuff


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


def plot_curve(args, number_field_periods, order_zeros, break_point):
    curve = cfs.get_curve(
        number_field_periods,
        order_zeros,
        np.array(args[:break_point]),
        np.array(args[break_point:]),
    )

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot(curve[0], curve[1], curve[2])
    ax.set_box_aspect([1, 1, 1])
    set_axes_equal(ax)
    plt.show()
