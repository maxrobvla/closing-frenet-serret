from cfs_helpers import close_curve, close_curve_alt, jac_loss

# from cfs_helpers import jac_loss

order_zeros = [1, 2]
break_point = 2
number_field_periods = 2
add_args = (number_field_periods, order_zeros, break_point)

# result = close_curve((1.0, -0.1, -1.0, 0.0), add_args, save=True)
# result = close_curve((1.0, 0.1, 0.0, 10.0), add_args, save=True)
result = close_curve((1.0, 1.0, 1.0, 1.0), add_args, save=True)


# result = jac_loss((1.0, 1.0, 1.0, 1.0), *add_args)

print(result)

# circle
# order_zeros = [0, 0]
# break_point = 1
# number_field_periods = 1  # irrelevant
# circle_params = (1.0, 0.0)
# print(loss(circle_params))
# save(circle_params)
