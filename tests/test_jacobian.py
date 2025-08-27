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
        positive_loss = cfsh.loss(params, *add_args)

        params = copy.deepcopy(starting_params)
        params[i] -= epsilon / 2.0
        negative_loss = cfsh.loss(params, *add_args)
        fd_jac.append((positive_loss - negative_loss) / epsilon)
    return fd_jac


analytic_jac = cfsh.jac_loss(starting_params, *add_args)

print("Finite Differences", f"(epsilon={1e-8}): ", get_fd_jacobian(1e-8), "\n")
print("Exact: ", analytic_jac)

# epsilons = np.logspace(2, -16, 20)
# losses_jacobian1 = []
# losses_jacobian2 = []
# losses_jacobian3 = []
# losses_jacobian4 = []
# for eps in epsilons:
#     fd_jac = get_fd_jacobian(eps)
#     losses_jacobian1.append(fd_jac[0])
#     losses_jacobian2.append(fd_jac[1])
#     losses_jacobian3.append(fd_jac[2])
#     losses_jacobian4.append(fd_jac[3])
# fig, ax = plt.subplots()
# ax.plot(epsilons, losses_jacobian1, ".", label="1")
# ax.axhline(analytic_jac[0], color="blue", label="1 (ana)")
# ax.plot(epsilons, losses_jacobian2, ".", label="2")
# ax.axhline(analytic_jac[1], color="red", label="2 (ana)")
# ax.plot(epsilons, losses_jacobian3, ".", label="3")
# ax.axhline(analytic_jac[2], color="green", label="3 (ana)")
# ax.plot(epsilons, losses_jacobian4, ".", label="4")
# ax.axhline(analytic_jac[3], color="black", label="4 (ana)")
# ax.legend()
# ax.set_xlabel(r"$\varepsilon$")
# ax.set_ylabel("FD gradient")
# ax.set_xscale("log")
# ax.set_yscale("symlog")
# fig.savefig("plots/scalar_loss_jac.png")
# plt.show()
