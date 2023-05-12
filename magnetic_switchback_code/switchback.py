from sympy import *
from sympy.abc import b
from sympy.vector import CoordSys3D, divergence, is_solenoidal, curl

# calculate divergence of vector field
def get_div(v):
    return divergence(v)

# calculate curl of vector field
def get_curl(v):
    return curl(v)

# calculate magnitude of vector field
def get_magn(v):
    return sqrt(v.dot(v))

# cartesian coordinates
N = CoordSys3D('N')
x, y, z = N.base_scalars()

# theta(x, y, z) and phi(x, y, z)
theta = b
phi = 2 * acot((x**2 + y**2 + z**2)**(1/2))

Bx = sin(theta) * cos(phi)
By = sin(theta) * sin(phi)
Bz = cos(theta)

B = Bx * N.i + By * N.j + Bz * N.k



print("\n")
print("divergence of B =  \n")
print(get_div(B))
print("\n")
print("In unicode: \n")
pprint(get_div(B), use_unicode=False)
print("\n")
print(f"Solenoidal? {is_solenoidal(B)}")
print("\n")
print("|B| = \n")
print(get_magn(B))
print("\n")
print("In unicode: \n")
pprint(get_magn(B), use_unicode=False)
print("\n")
print("curl of B =  \n")
print(get_curl(B))
print("\n")