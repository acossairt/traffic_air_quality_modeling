import sys
print(sys.executable)
import numpy as np
import matplotlib.pyplot as plt

# Parameters
f = 1
h = 0.5
l = 2

x = np.linspace(0, 2, 500)

# Curve 1: Generalized speed-density (purple)
y1 = f * np.exp(-np.log(2) * (x / h) ** l)

# Curve 2: Drake's equation (green)
y2 = f * np.exp(-0.5 * (x / h) ** 2)

plt.figure(figsize=(6, 6))
plt.plot(x, y1, color="purple", linewidth=2,
         #label="General sigmoid function"
         label=r"$V_{i,j,k}^f \cdot e^{-\ln(2)\left(\frac{C_{i,j,k}}{C_{i,j,k}^{half}}\right)^{\lambda}}$"
)
plt.plot(x, y2, color="green", linewidth=2,
         #label="Drake's model (1967)"
         label=r"$V_{i,j,k}^f \cdot e^{-\frac{1}{2}\left(\frac{K_{i,j,k}}{K_{j}}\right)^{\lambda}}$"
)

plt.xlabel(r"density ($C_{i,j,k}$ or $K$)", fontsize=14)
plt.ylabel(r"average speed ($V_{i,j,k}^{avg}$)", fontsize=14)
plt.xlim(0, 1.1)
plt.title("Speed-density curves", fontsize=16)
plt.legend(fontsize=14, loc="upper right")
plt.tight_layout()
plt.savefig("images/alternative_drake_new_eqs.png", dpi=150)
plt.show()
 