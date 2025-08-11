import os

import numpy as np
from scipy.optimize import minimize

import build.closing_frenet_serret as cfs


def loss(args):
    loss = cfs.loss(
        number_field_periods,
        order_zeros,
        np.array(args[:break_point]),
        np.array(args[break_point:]),
    )
    return loss


def save(args):
    cfs.save_curve(
        "solution",
        number_field_periods,
        order_zeros,
        np.array(args[:break_point]),
        np.array(args[break_point:]),
    )


# order_zeros = [1, 2]

# break_point = 2
# number_field_periods = 1
# result = minimize(loss, (1.0, 1.0, 1.0, 1.0, 1.0))
# result = minimize(loss, (0.0, -1.0, -1.0, 1.0, 1.0))
# result = minimize(loss, (0.0, 1.0, 1.0, 1.0, 1.0))
# result = minimize(loss, (1.0, 0.0, 1.0, 1.0, 1.0))
# result = minimize(loss, (1.0, 1.0, 0.0, 1.0, 1.0))  # similar to no. 6
# result = minimize(loss, (1.0, 1.0, 1.0, 0.0, 1.0))
# result = minimize(loss, (1.0, 1.0, 1.0, 1.0, 0.0))

# order_zeros = [0, 0]

# break_point = 2
# number_field_periods = 1
# result = minimize(loss, (1.0, 0.0, 0.0, 0.0))

# break_point = 3
# number_field_periods = 1
# result = minimize(loss, (1.0, 0.0, 1.0, 1.0, 1.0, 1.0)) # similar to no. 5 + 6 from first block

# break_point = 3
# number_field_periods = 2
# result = minimize(
#     loss, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
# )  # interesting because of overlap
# result = minimize(loss, (1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0)) # a bit more normal
# result = minimize(loss, (1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0)) # intersects only in a single point

# break_point = 3
# number_field_periods = 3
# result = minimize(loss, (1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0))  # two triangles
# result = minimize(loss, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)) # similar to last one
# result = minimize(loss, (0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0))  # could not be optimized
# result = minimize(loss, (1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0))  # has some loops
# result = minimize(
#     loss, (1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0)
# )  # looks like three bladed clover from one side

# break_point = 4
# number_field_periods = 4
# result = minimize(loss, (1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0)) # visibly not fully closed, flower
# result = minimize(
#     loss, (1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
# )  # stacked saddles
# result = minimize(
#     loss, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
# )  # twirly butterfly

# break_point = 2
# number_field_periods = 3
# result = minimize(loss, (-100, 1.0, 1.0, 1.0))  # twirly hexagon

order_zeros = [2, 1]
break_point = 2
number_field_periods = 1
result = minimize(loss, (1.0, 1.0, 0.0, 1.0))

print(result)

save(result.x)


# circle
# order_zeros = [0, 0]
# break_point = 1
# number_field_periods = 1  # irrelevant
# circle_params = (1.0, 0.0)
# print(loss(circle_params))
# save(circle_params)
