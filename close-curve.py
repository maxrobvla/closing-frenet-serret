import time

from cfs_helpers import close_curve, close_curve_with_jac, jac_loss

# from cfs_helpers import jac_loss

order_zeros = [1, 2]
break_point = 2
number_field_periods = 2
add_args = (number_field_periods, order_zeros, break_point)

# result = close_curve((1.0, -0.1, -1.0, 0.0), add_args, save=True)
# result = close_curve((1.0, 0.1, 0.0, 10.0), add_args, save=True)
#
# result = close_curve((1.0, 1.0, 1.0, 1.0), add_args, save=False)
# print("Nelder Mead success", result.success)
# start = result.x
# print(result.x)
# start[0] += 0.01
# print(start)
t0 = time.time()
# result = close_curve_with_jac(start, add_args, save=True)
result = close_curve_with_jac((0.1, 0.1, 0.1, 0.1), add_args, save=True)
print(time.time() - t0)


# result = jac_loss((1.0, 1.0, 1.0, 1.0), *add_args)

print(result)

# circle
# order_zeros = [0, 0]
# break_point = 1
# number_field_periods = 1  # irrelevant
# circle_params = (1.0, 0.0)
# print(loss(circle_params))
# save(circle_params)
