import matplotlib
import matplotlib.pyplot as plt

x, u = [], []
with open("./data.txt", "r") as file:
    for line in file:
        xi, ui = map(float, line.split())
        x.append(xi)
        u.append(ui)

plt.plot(x, u, label="u(x)", marker="o")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Wykres u(x)")
plt.legend()
plt.grid()
plt.show()