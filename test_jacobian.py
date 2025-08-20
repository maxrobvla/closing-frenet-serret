import copy

import cfs_helpers as cfsh

order_zeros = [1, 2]
break_point = 2
number_field_periods = 2
add_args = (number_field_periods, order_zeros, break_point)
starting_params = [5.3, 3.0, 1.8, -0.65]


def get_fd_jacobian(epsilon):
    fd_jac = []
    for i in range(len(starting_params)):
        # central differences
        params = copy.deepcopy(starting_params)
        params[i] += epsilon / 2.0
        positive_loss = cfsh.vector_loss(params, *add_args)

        params = copy.deepcopy(starting_params)
        params[i] -= epsilon / 2.0
        negative_loss = cfsh.vector_loss(params, *add_args)
        fd_jac.append((positive_loss - negative_loss) / epsilon)
    return fd_jac


analytic_jac = cfsh.jac_vector_loss(starting_params, *add_args)

for eps in [10 ** (-i) for i in range(8, 10)]:
    print("Finite Differences", f"(epsilon={eps}): ", get_fd_jacobian(eps), "\n")
print("Exact: ", analytic_jac)
