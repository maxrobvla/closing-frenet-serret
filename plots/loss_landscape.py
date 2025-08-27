import os
import sys

import matplotlib.pyplot as plt
import numpy as np

path_root = os.path.dirname(__file__)
sys.path.append(path_root + "/..")

import cfs_helpers as cfsh

order_zeros = [1, 2]
break_point = 2
number_field_periods = 2
add_args = (number_field_periods, order_zeros, break_point)
# starting_params = [5.3, 3.0, 1.8, -0.65]
starting_params = [1.0, 1.0, 1.0, 1.0]


def scan_loss(idx, devs=None, vector=True):
    if devs is None:
        devs = np.linspace(-1, 1, 500)

    # starting_params = [5.3, 3.0, 1.8, -0.65]
    starting_params = [1.0, 1.0, 1.0, 1.0]

    x = np.ones_like(devs) * devs + np.ones_like(devs) * starting_params[idx]

    if vector:
        normal_dot_product = []
        binormal_dot_product = []
    else:
        losses = []

    for i in range(len(devs)):
        if i % 50 == 0:
            print(i)
        args = starting_params
        args[idx] = x[i]

        if vector:
            loss = cfsh.vector_loss(args, *add_args)
            normal_dot_product.append(loss[-2])
            binormal_dot_product.append(loss[-1])
        else:
            losses.append(cfsh.loss(args, *add_args))
    if vector:
        return (x, [normal_dot_product, binormal_dot_product])
    else:
        return (x, losses)


def get_tangent(x0, y0, deriv):
    def tangent(x):
        return deriv * (x - x0) + y0

    return tangent


fig, ax = plt.subplots(4, 2, figsize=(7, 15))
y0 = cfsh.vector_loss(starting_params, *add_args)
deriv = cfsh.jac_vector_loss(starting_params, *add_args)
for i in range(4):
    x, [normal_dot_product, binormal_dot_product] = scan_loss(i)
    xt = x[(len(x) // 2 - 25) : (len(x) // 2 + 25)]
    ax[i][0].plot(x, normal_dot_product, ".")
    ax[i][0].plot(xt, get_tangent(starting_params[i], y0[-2], deriv[i][-2])(xt))
    ax[i][0].axvline(starting_params[i])
    ax[i][1].plot(x, binormal_dot_product, ".")
    ax[i][1].plot(xt, get_tangent(starting_params[i], y0[-1], deriv[i][-1])(xt))
    ax[i][1].axvline(starting_params[i])
fig.savefig("plots/landscape.png")
plt.show()

# fig, ax = plt.subplots(1, 4, figsize=(10, 5))
# for i in range(4):
#     x, losses = scan_loss(i, vector=False)
#     ax[i].plot(x, losses, ".")
#     ax[i].axvline(starting_params[i])
# fig.savefig("plots/landscape.png")
# plt.show()
