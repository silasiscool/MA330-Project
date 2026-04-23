# Generate fourier coeifficients for a given function

# Modules
from sympy import Symbol, integrate, sin, cos, pi

# Setup
t = Symbol('t', real=True)
n = Symbol('n', real=True, positive=True, integer=True)

# Cofficients
def a0(f,P):
    L = P/2
    return 1/(2*L)*integrate(f,(t,-L,L))

def an(f,P):
    L = P/2
    return 1/L*integrate(f*cos(n*pi*t/L),(t, -L, L))

def bn(f,P):
    L = P/2
    return 1/L*integrate(f*sin(n*pi*t/L),(t,-L,L))
