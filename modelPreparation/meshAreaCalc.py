import trimesh
import math
import numpy
from neuron import h, rxd
from neuron.units import ms, mV
import pandas as pd
import time as clock
import cmath

h.load_file("stdrun.hoc")

def find_radius(height, surface_area_um):
    # Coefficients of the quadratic equation
    a = 2 * cmath.pi
    b = 2 * cmath.pi * height
    c = -surface_area_um

    # Calculate the discriminant
    discriminant = cmath.sqrt(b**2 - 4 * a * c)

    # Use the quadratic formula to find the roots
    root1 = (-b + discriminant) / (2 * a)
    root2 = (-b - discriminant) / (2 * a)

    return root1, root2

def meshItUp():
    pi = math.pi
    # Load the mesh for the axons here
    mesh = trimesh.load("datafiles/morphologyData/Axon_meshes/648518346475217761.obj")
    
    print("mesh loaded")
    print(type(mesh))
    surface_area = mesh.area

    print('Surface area:', surface_area, 'nm^2')
    surface_area_um = surface_area / (4.3*4.3 * 1000**2)
    print('Surface area:', surface_area_um, 'µm^2')

    r = 1
    height = (surface_area_um - 2*pi*(r**2)) / (2*pi*r)
    print(height)


    # Example usage
    height2 = 110.6
    roots = find_radius(height2, surface_area_um)

    print("Possible values for radius:")
    for root in roots:
        print(root.real)

meshItUp()