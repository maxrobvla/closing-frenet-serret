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


def close_curve(args, add_args, save=True):
    result = minimize(loss, args, args=add_args, method="Nelder-Mead")
    if save:
        save_curve(result.x, add_args)
    return result


def close_curve_alt(args, add_args, save=True):
    result = minimize(loss, args, args=add_args, method="BFGS")
    if save:
        save_curve(result.x, add_args)
    return result


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
