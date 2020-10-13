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
        #self.taper_ratio = data['taper_ratio']
        self.V = data['V']
        self.p = data['p']
        self.R_outer = data['R_outer']
        self.length = data['length']
        self.G = data['G']
        self.E = data['E']
        self.allowed_ax_stress = data['allowed_ax_stress']
        self.delta = data['delta']
        self.R_inner = data['R_inner']
        self.layup_density = data['layup_density']
        
d = Yaml_Parse()

def get_cord_y(y):
    """
    calculates cord length at a location y
    """
    if y < 11:
        return 1.89
    
    else:
        return 1.89 - (1.89-.95)/9 * (y-11)

class Statics:
    def __init__(self):
        self.kp = self.set_kp()
        self.q = self.get_q()
    
    def get_q(self):
        """
        return flight dynamic pressure
        """
        return .5*d.p*d.V**2
    
    def set_kp(self):
        sum_ = 0
        for i in range(0, int(.5*d.span/d.delta)):
            sum_ += get_cord_y(i*d.delta)*d.delta
        return d.N*d.center_weight/(2*sum_)
    
const = Statics()
    
def radius(y):

    if y <= 7.5:
        return 0.04*get_cord_y(7.5)
    elif y <= 15:
        return 0.04*get_cord_y(15)
    else:
        return 0.04*get_cord_y(d.span/2)

def distributed_load_y(y):
    """
    calculates distributed load at a 
    location y in N/m
    """
    return const.kp * get_cord_y(y)

    
def shear_root():
    """
    calculates shear at the root in N
    """
    return -d.N*d.center_weight/2

#def moment_root():
#    """
#    calculates moment at the root in N*M
#    """
#    return (d.center_weight*d.N*d.span/12)*\
#    (1+2*d.taper_ratio)/(1+d.taper_ratio)


###### tube sizing functions #########


def get_mass(thickness):
    """
    computes mass (kg) of enture tube spar caps given discretizated thickness
    """
    weight = 0
    for index, t in enumerate(thickness):
        t *= .01 # convert from cm to m
        weight += (1/3)*math.pi*2*radius(index*d.delta)*t*d.delta*d.layup_density * 2 
    print('weight of spar caps is: ', weight, ' kg')
    

def size_bending():
    """
    strength thickness, distribution load, shear, moment, deflection
    all as functions of spanwise distance.
    
    
    returns each in a data dictionary
    """
    print(const.kp, 'aaaaaaaaaaaa')
    
    y_sec = [y*d.delta for y in range(0,int(d.span/(2*d.delta)))] # discrete points along the wing
    
    # get load distribution
    dist_load = [distributed_load_y(y) for y in y_sec]
    
#    sum01 = 0
#    for i in dist_load:
#        sum01 += d.delta * i
#    print(sum01, '!!!!!!!!!!!!!!!!')
    
    # get shear distribution
    shear = [shear_root()]
    #shear = [0]
    for y in range(1,len(y_sec)):
        shear.append(shear[y-1] + dist_load[y]*d.delta)
    
    # get moment distribution
    #moment = [moment_root()]
    moment=[13676.51-2750]
    #moment=[13676.51]
    for y in range(1,len(y_sec)):
        moment.append(moment[y-1] + shear[y]*d.delta)
        
    # get cap thickness for strength constraint
    #computes thickness of tube (cm) at y based on bending strength sizing eq 40
    thickness = [100*moment[index]/(1.913*(radius(y))**2 * d.allowed_ax_stress*1000000)\
                 for index,y in enumerate(y_sec)] # in cm

    # print weight of spar caps
    get_mass(thickness)
        
    # moment of inertia distribution
    I = []
    for y in range(len(y_sec)):
        I.append(radius(y*d.delta)**3*thickness[y]*.01*1.913) # .001 converts from cm to m
        
    # calculates deflection distribution
    theta = [0]
    deflection=[0]
    for y in range(1,len(y_sec)):
        theta.append(theta[y-1]+d.delta*moment[y]/(d.E*1000000*I[y]))
    for y in range(1,len(y_sec)):
        deflection.append(deflection[y-1] + d.delta*theta[y])
        
    
    return {'y':y_sec, 'thickness':thickness, 'load': dist_load,\
            'shear':shear, 'moment':moment, 'deflection':deflection, 'I':I}
    

#def bending_deflect():
#    """
#    computes thickness of tube (cm) based on bending deflection sizing
#    eq 46
#    """ 
#    (1/8 * 3000 * d.length**3)/(math.pi*d.R_outer**2*d.deflection*d.E*1000000)
    
###### plotting functions ########
    
data = size_bending()

Ei = [d.E*i*1000000 for i in data['I']]
plus_5 = [i +5 for i in data['y']]
plt.plot(plus_5, Ei)
plt.title('EI vs spanwise distance')
plt.xlabel('spanwise distance (m)')
plt.ylabel('Ei')
plt.grid()
plt.show()

cord = [get_cord_y(i) for i in data['y']]
plt.plot(data['y'], cord)
plt.title('cord vs spanwise distance')
plt.xlabel('spanwise distance (m)')
plt.ylabel('cord')
plt.grid()
plt.show()

plt.plot(data['y'], data['load'])
plt.title('Lifting Load vs Spanwise Distance', fontsize = '18')
plt.xlabel('Spanwise Distance (m)', fontsize = '18')
plt.ylabel('Load (N/m)', fontsize = '18')
plt.grid()
plt.show()

plt.plot(data['y'], data['shear'])
plt.title('Shear vs spanwise distance')
plt.xlabel('spanwise distance (m)')
plt.ylabel('Shear (N)')
plt.grid()
plt.show()

plt.plot(data['y'], data['moment'])
plt.title('Moment vs Spanwise Distance', fontsize = '18')
plt.xlabel('Spanwise Distance (m)', fontsize = '18')
plt.ylabel('Moment (N*m)', fontsize = '18')
plt.grid()
plt.show()
    
#plt.plot(data['y'], data['thickness'])
#plt.title('Spar Cap Thickness vs Spanwise Distance', fontsize = '18')
#plt.xlabel('Spanwise Distance (m)', fontsize = '18')
#plt.ylabel('Thickness (cm)', fontsize = '18')
#plt.grid()
#plt.show()

plt.plot(data['y'], data['deflection'], color="red")
plt.title('Wing Deflection vs Spanwise Distance', fontsize = '18')
plt.xlabel('Spanwise Distance (m)', fontsize = '18')
plt.ylabel('Deflection (m)', fontsize = '18')
plt.grid()
plt.show()

fig,ax = plt.subplots()
ax.plot(data['y'], data['thickness'], color="red",linewidth=7.0)
plt.title('Spar Cap Thickness vs Spanwise Distance', fontsize = '18')
ax.set_xlabel('Spanwise Distance (m)', fontsize = '18')
ax.set_ylabel('Thickness (cm)', fontsize = '18', color="red")

ax2=ax.twinx()
rad = [radius(y) for y in data['y']]
ax2.plot(data['y'], rad, color="blue")
#plt.title('Spar Inner Radius', fontsize = '18')
#ax2.xlabel('Spanwise Distance (m)', fontsize = '18')
ax2.set_ylabel('Inner Radius of Spar (m)', fontsize = '18', color="blue")
plt.grid()
plt.show()
    