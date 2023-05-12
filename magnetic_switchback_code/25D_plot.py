import matplotlib.pyplot as plt
from matplotlib import cm
import sympy as sp
import numpy as np




# define constants
B_magn = 5
B0 = 1
psi_0 = 2
x1 = -1
x2 = 1
y1 = 1
y2 = -1
lx = 1
ly = 2

x, y, z = sp.symbols('x y z')


# get magnetic scalar potential
def get_psi(x, y, psi_0, x1, x2, y1, y2, lx, ly):
    r1_sqrd = ((x - x1) / lx)**2 + ((y - y1) / ly)**2
    r2_sqrd = ((x - x2) / lx)**2 + ((y - y2) / ly)**2
    return - psi_0 * (sp.exp(-r1_sqrd) - sp.exp(-r2_sqrd))

# get z-component of magnetic field
def get_Bz(B_magn, Bx, By):
    return sp.sqrt(B_magn**2 - Bx**2 - By**2)

# calculate divergence of vector field
def get_div(Bx, By, Bz):
    return sp.diff(Bx, x) + sp.diff(By, y) + sp.diff(Bz, z)

# calculate curl of vector field
def get_curl(Bx, By, Bz):
    return [sp.diff(Bz, y) - sp.diff(By, z), sp.diff(Bx, z) - sp.diff(Bz, x), sp.diff(By, x) - sp.diff(Bx, y)]



psi = get_psi(x, y, psi_0, x1, x2, y1, y2, lx, ly)

Bx = sp.diff(psi, y)
By = B0 - sp.diff(psi, x)
Bz = get_Bz(B_magn, Bx, By)


print("\n")
print("Divergence of B = \n")
print(f"{get_div(Bx, By, Bz)} \n")
print("Curl of B = \n")
print(f"{get_curl(Bx, By, Bz)} \n")


# convert to numerical via lambdify()
Bx = sp.lambdify([x, y], Bx)
By = sp.lambdify([x, y], By)
Bz = sp.lambdify([x, y], Bz)



N = 1000
nx = np.linspace(-5, 5, N)
ny = np.linspace(-5, 5, N)
nz = np.linspace(0, 1, N)

X, Y = np.meshgrid(nx, ny)
U = Bx(X, Y)
V = By(X, Y)
W = Bz(X, Y)

W_normed = W / W.max(axis=0)

fig = plt.figure(figsize = (5, 5))
plt.streamplot(X, Y, U, V, color = -W_normed, density = 2, linewidth = 1, arrowsize = 0.5, cmap=plt.cm.jet)
plt.title("Streamplot 2D")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.show()

N = 15
nx = np.linspace(-5, 5, N)
ny = np.linspace(-5, 5, N)
nz = np.linspace(0, 1, N)

X, Y, Z = np.meshgrid(nx, ny, nz)
U = Bx(X, Y)
V = By(X, Y)
W = Bz(X, Y)


# color by  angle between U and W
c = np.arctan(U, W)
# Flatten and normalize
c = (c.ravel() - c.min()) / c.ptp() 
# Repeat for each body line and two head lines
c = np.concatenate((c, np.repeat(c, 2)))

# Colormap
c = plt.cm.hsv(c)
fig = plt.figure(figsize = (5, 5))
ax = fig.add_subplot(projection = '3d')
ax.quiver(X, Y, Z, U, V, W, colors = c, length = 0.3, normalize = True) 
plt.title("Streamplot 3D")
ax.set_xlabel(r'$X$')
ax.set_ylabel(r'$Y$')
ax.set_zlabel(r'$z$')
plt.show()
