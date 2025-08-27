import os
import pickle
import signal
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

path_root = os.path.dirname(__file__)

sys.path.append(os.path.join(path_root, ".."))

from cfs_helpers import close_curve, close_curve_with_jac

total_number_params = 4
break_point = 2
order_zeros = [1, 2]
number_field_periods = 1
add_args = (number_field_periods, order_zeros, break_point)

gen = np.random.default_rng(np.random.PCG64())
samples = []
for i in range(100):
    random_signs = np.random.choice([-1, 1], (total_number_params))
    random_positive_values = np.exp(gen.uniform(0.0, 3.0, total_number_params)) - 1
    samples.append(random_signs * random_positive_values)


# a safety mechanism for very long optimizations
def handler(signum, frame):
    raise TimeoutError("Operation timed out")


signal.signal(signal.SIGALRM, handler)
timeout_duration = 10

method = "Nelder-Mead"
minimized_values = []
times = []
fail_counter = 0

alt_method = "BFGS"
minimized_values_alt = []
times_alt = []
fail_counter_alt = 0
counter = 0
for sample in samples:
    print(counter)
    # main method
    signal.alarm(timeout_duration)
    start = time.time()
    result = None
    t = None
    try:
        result = close_curve(sample, add_args, save=False).fun
        t = time.time() - start
    except TimeoutError:
        fail_counter += 1
        print("Optimization took too long")
    finally:
        if result is not None and t is not None:
            minimized_values.append(result)
            times.append(t)
            print("Success")
        signal.alarm(0)

    # second method
    signal.alarm(timeout_duration)
    start = time.time()
    result = None
    t = None
    try:
        result = close_curve_with_jac(sample, add_args, save=False).fun
        t = time.time() - start
    except TimeoutError:
        fail_counter_alt += 1
        print("Optimization took too long")
    finally:
        if result is not None and t is not None:
            minimized_values_alt.append(result)
            times_alt.append(t)
            print("Success")
        signal.alarm(0)
        counter += 1

# print some statistics
print(f"method 1: {method}")
print(f"runs longer than {timeout_duration}s: {fail_counter}")
print(f"run time: {np.mean(times)}({np.std(times)})")
print(f"loss: {np.mean(minimized_values)}({np.std(minimized_values)})\n")

print(f"method 2: {alt_method}")
print(f"runs longer than {timeout_duration}s: {fail_counter_alt}")
print(f"run time: {np.mean(times_alt)}({np.std(times_alt)})")
print(f"loss: {np.mean(minimized_values_alt)}({np.std(minimized_values_alt)})")

with open(path_root + "/../data/benchmark", "wb") as file:
    pickle.dump([[times, minimized_values], [times_alt, minimized_values_alt]], file)

fig, ax = plt.subplots()
ax.scatter(times, minimized_values, label=f"{method}", s=1)
ax.scatter(times_alt, minimized_values_alt, label=f"{alt_method}", s=1)
ax.set_xlabel(r"$t$ in seconds")
ax.set_ylabel(r"achieved value of loss")
ax.set_xlim(-0.1, 5)
ax.legend()
ax.set_yscale("log")
fig.savefig(os.path.join(path_root, "test_optimizer.png"))
plt.show()
