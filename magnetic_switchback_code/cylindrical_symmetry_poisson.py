import sympy as sp
import numpy as np

rho, z = sp.symbols('rho z')

B = 1

def get_psi(rho, z):
    return 1/(1 + sp.exp(-rho**2 - z**2))

def get_B_phi(rho, z):
    return 0

# calculate Laplacian of vector field
def get_laplacian(psi, B_phi):
    return (1/rho) * sp.diff(rho * sp.diff(psi, rho), rho) + sp.diff(B_phi, z, 2)

def get_g(psi, B_phi):
    return - B * sp.diff(sp.ln(2*B**2 - (sp.diff(psi, rho)**2 + B_phi**2)), z)


psi = get_psi(rho, z)
B_phi = get_B_phi(rho, z)

print(f"Laplacian = {get_laplacian(psi, B_phi)}")
print(f"g = {get_g(psi, B_phi)}")
print(get_laplacian(psi, B_phi) == get_g(psi, B_phi))