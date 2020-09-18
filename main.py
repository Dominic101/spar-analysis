# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 19:54:20 2020

@author: Dominic
"""

import yaml
import math
import matplotlib.pyplot as plt

with open('params.yaml') as f:
    data = yaml.load(f, Loader=yaml.FullLoader)

# compute static parameters that are
# used in multiple calculations

class Yaml_Parse:
    def __init__(self):
        self.center_weight = data['center_weight']
        self.N = data['N']
        self.span = data['span']
        self.wing_area = data['wing_area']
        self.taper_ratio = data['taper_ratio']
        self.m = float(data['m']) # TODO see why the negative makes it a string
        self.b = data['b']
        self.r_i = data['r_i']
        
d = Yaml_Parse()

class Statics:
    def __init__(self):
        self.kp = self.get_kp()
        self.root_cord = self.get_root_cord()
        self.m = self.get_m()
        self.b = self.get_b()
        
    def get_kp(self):
        """
        calculate proportionality constant
        """
        return d.center_weight*d.N/d.wing_area
    
    def get_root_cord(self):
        return d.wing_area/d.span*(2/(1+d.taper_ratio))
    
    def get_m(self):
        return 2*self.kp*self.root_cord*(d.taper_ratio-1)/d.span
    
    def get_b(self):
        return self.root_cord*self.kp
    
const = Statics()

def get_cord_y(y):
    """
    calculates cord length at a location y
    """
    pass

def distributed_load_y(y):
    """
    calculates distributed load at a 
    location y in N/m
    """
    return const.m*y+const.b

def r_0_y(y):
    """
    return r_0 at y
    """
    return d.m*y + d.b
    
def shear_root():
    """
    calculates shear at the root in N
    """
    return -d.N*d.center_weight/2

def moment_root():
    """
    calculates moment at the root in N*M
    """
    return (d.center_weight*d.N*d.span/12)*\
    (1+2*d.taper_ratio)/(1+d.taper_ratio)
    
def shear_y(y):
    """
    calculates shear at a location y in N
    """
    return const.m/2*y**2 + const.b*y + shear_root()

def get_I(y):
    return math.pi/4*(r_0_y(y)**4 - d.r_i**4)

def moment_y(y):
    """
    calculates moment at a location y in N*M
    """
    return const.m/6*y**3 + const.b/2*y**2 + shear_root()*y + moment_root()

def max_shear_stress_y(y):
    """
    calculates max shear stress at a location y in Mpa
    """
    return -shear_y(y)*(2/3)*(r_0_y(y)**3 - d.r_i**3)/(10**6*get_I(y)\
    *2*(r_0_y(y)-d.r_i))

def max_axial_stress_y(y):
    """
    calculates max axial stress at a location y in Mpa
    """
    return moment_y(y)*r_0_y(y)/(get_I(y)*10**6)

def max_axial_stress_wing():
    """
    calculates max axial stress thoughout 
    the entire wing in Mpa
    """
    # TODO: maybe use scipy for this
    delta = .01 # distance between points
    points = [i*delta for i in range(int(.5*d.span/delta))]
    return max([max_axial_stress_y(i) for i in points])

def max_shear_stress_wing():
    """
    calculates max shear stress throughout
    the entire wing in Mpa
    """
    # TODO: maybe use scipy for this
    delta = .01 # distance between points
    points = [i*delta for i in range(int(.5*d.span/delta))]
    return max([max_shear_stress_y(i) for i in points])

###### plotting functions ########

def gen_plot_data():
    delta = .01 # distance between points
    points = [i*delta for i in range(int(.5*d.span/delta))]
    
    shears = [shear_y(i) for i in points]
    moments = [moment_y(i) for i in points]
    shear_s = [max_shear_stress_y(i) for i in points] # shear stress
    axial_s = [max_axial_stress_y(i) for i in points] # axial stress
    
    plt.plot(points, shears)
    plt.title('Shear vs location on wing')
    plt.xlabel('distance from center (m)')
    plt.ylabel('Shear (N)')
    plt.grid()
    plt.show()
    
    plt.plot(points, moments)
    plt.title('moment vs location on wing')
    plt.xlabel('distance from center (m)')
    plt.ylabel('Moment (N*m)')
    plt.grid()
    plt.show()
    
    plt.plot(points, shear_s)
    plt.title('Shear stress vs location on wing')
    plt.xlabel('distance from center (m)')
    plt.ylabel('shear stress (MPa)')
    plt.grid()
    plt.show()
    
    plt.plot(points, axial_s)
    plt.title('axial stress vs location on wing')
    plt.xlabel('distance from center (m)')
    plt.ylabel('axial stress (MPa)')
    plt.grid()
    plt.show()
    
if True:
    gen_plot_data()
    #print(max_axial_stress_wing())
    
    
    