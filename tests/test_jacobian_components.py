import copy

import matplotlib.pyplot as plt
import numpy as np

import cfs_helpers as cfsh

order_zeros = [1, 2]
break_point = 2
number_field_periods = 1
add_args = (number_field_periods, order_zeros, break_point)
starting_params = [5.3, 3.0, 1.8, -0.65]
# starting_params = [1.0, 1.0, 1.0, 1.0]


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

# for eps in [10 ** (-i) for i in range(10, 12)]:
# print("Finite Differences", f"(epsilon={eps}): ", get_fd_jacobian(eps), "\n")
print("Finite Differences", f"(epsilon={1e-8}): ", get_fd_jacobian(1e-8), "\n")
print("Exact: ", analytic_jac)

# parameter_idx = 3

# epsilons = np.logspace(2, -16, 20)
# dot_product1_jac = []
# dot_product2_jac = []
# test_jac = []
# for eps in epsilons:
#     fd_jac = get_fd_jacobian(eps)
#     dot_product1_jac.append(fd_jac[parameter_idx][-2])
#     dot_product2_jac.append(fd_jac[parameter_idx][-1])
#     test_jac.append(fd_jac[parameter_idx][0])
# dot_product1_jac = np.array(dot_product1_jac)
# dot_product2_jac = np.array(dot_product2_jac)
# fig, ax = plt.subplots()
# ax.plot(epsilons, dot_product1_jac, ".", label="normal")
# ax.axhline(analytic_jac[parameter_idx][-2], label="normal (ana)", color="black")
# ax.plot(epsilons, dot_product2_jac, ".", label="binormal")
# ax.axhline(analytic_jac[parameter_idx][-1], label="binormal (ana)", color="red")
# ax.plot(epsilons, test_jac, ".", label="position diff")
# ax.axhline(analytic_jac[parameter_idx][0], label="position diff (ana)")
# ax.legend()
# ax.set_xlabel(r"$\varepsilon$")
# ax.set_ylabel("FD gradient")
# ax.set_xscale("log")
# ax.set_yscale("symlog")
# fig.savefig("plots/dot_products_jac.png")
# plt.show()
