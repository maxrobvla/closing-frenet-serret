import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

script_dir = os.path.dirname(__file__)
data = pd.read_csv(os.path.join(script_dir, "..", "data", "solution.dat"))

tangents = [data[f"t{i}"].to_numpy() for i in range(1, 4)]
arc_lengths = data["l"].to_numpy()

curve = np.array([data[f"p{i}"].to_list() for i in range(1, 4)])

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.plot(curve[0], curve[1], curve[2])
plt.show()
