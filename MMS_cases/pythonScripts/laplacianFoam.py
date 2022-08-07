#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 20:45:25 2021

@author: pc
"""

# Script developed by:
#      Computational Rheology Group at the Institute for Polymers and Composites, University of Minho, Portugal, 2022
# If you find any problems/mistakes or have any improvement recommendations, please contact: info@crheo.org


# import sympy
import sympy as sym 

sym.init_printing(pretty_print=False)

# Create symbolic variables
x, y, z = sym.symbols("x y z")

# Function to compute the gradient of a scalar
def grad (scalarQuantity):
    result = sym.Matrix([0,0,0])
    result[0] = sym.diff(scalarQuantity, x)
    result[1] = sym.diff(scalarQuantity, y)
    result[2] = sym.diff(scalarQuantity, z)
    return result

# Function to compute the divergence of a vector
def div (vectorQuantity):
    result = sym.diff(vectorQuantity[0], x) + sym.diff(vectorQuantity[1], y) + sym.diff(vectorQuantity[2], z)
    return result


# Manufactured solution
T = 150*(sym.cos(x*x + y*y) + 1.5)

# Diffusion coefficient [m^2/s]
DT = 1e-5

# Compute the source term
S_T = - div(DT*grad(T))

## Generate C code

# Expand powers
from sympy.codegen.rewriting import create_expand_pow_optimization
expandPowers = create_expand_pow_optimization(3)

# Generate C-code
code = sym.ccode(expandPowers (S_T) )

# Dictionary to replace with OF syntax
replaceDictionary = {
                        "sin":"Foam::sin",
                        "cos": "Foam::cos",
                        "exp": "Foam::exp"
                      }

# Loop through the generated C code and replace the entries in OF syntax
for key, value in replaceDictionary.items():
    code = code.replace(key, value)

# Print code for the user
print(code)

