import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


data = np.genfromtxt("temperature_profile.csv", delimiter=',', names=True)
x = data["arc_length"]*1e3
T = data["T"]

fig = plt.figure(figsize=(6.4, 4))
ax = fig.add_subplot(111)
plt.plot(x, T, color="black")

h = 1250 - 380
rect1 = Rectangle((0, 370), 6, h, color="tab:grey", alpha=0.5)
ax.add_patch(rect1)

rect2 = Rectangle((6, 370), 1, h, color=(228/255, 146/255, 64/255), alpha=0.5)
ax.add_patch(rect2)

rect3 = Rectangle((7, 370), 1.55, h, color=(180/255, 95/255, 6/255), alpha=0.5)
ax.add_patch(rect3)

plt.annotate("W", (5, 1100), weight="bold", fontsize=15)
plt.annotate("Cu", (6.1, 1100), weight="bold", fontsize=15)
plt.annotate("CuCrZr", (7.1, 1100), weight="bold", fontsize=15)

plt.ylim(370, 1200)
plt.xlim(0, 8.5)

plt.xlabel("Depth (mm)", fontsize=15)
plt.ylabel("T (K)", fontsize=15)
plt.tight_layout()
plt.savefig("temperature_profile.pdf")
plt.show()
