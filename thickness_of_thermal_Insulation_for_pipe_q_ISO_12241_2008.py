# -*- coding: utf-8 -*-
"""
Created on 20200525
@author: weililai@foxmail.com

    Method reference standard ISO 12241-2008
    Except for temperature (Celsius), all other units are in the International System of Units.
#1.Surface coefficient of heat transfer
    h_se = h_r + h_cv
        h_r = alpha_r*C_r
            alpha_r = ((theta_se+273.15)**4 - (theta_a+273.15)**4)/((theta_se+273.15) - (theta_a+273.15))
            C_r = epsilon*sigma
                sigma = 5.67*10**-8
        h_cv = 1.32*((theta_se - theta_a)/D_e)**0.25  #For vertical pipes for laminar, free convection, inside buildings
        h_cv = 1.74*(theta_se - theta_a)**(1/3)       #For vertical pipes for turbulent, free convection, inside buildings
        h_cv = 1.25*((theta_se - theta_a)/D_e)**0.25  #For horizontal pipes for laminar airflow, inside buildings
        h_cv = 1.21*((theta_se - theta_a)/D_e)**(1/3) #For horizontal pipes for turbulent airflow, inside buildings
        h_cv = 8.1*10**-3/D_e + 3.14*(v/D_e)**0.5     #For horizontal and vertical pipes that are outside buildings for laminar airflow
        h_cv = 8.9*v**0.9/D_e**0.1                    #For horizontal and vertical pipes that are outside buildings for turbulent airflow

#2.Linear density of heat flow rate between the externa surface of the insulation layer and the ambient
    q_l = (theta_se-theta_a)*pi*D_e*h_se
#3.Linear density of heat flow rate in insulation materials
    q_l = (theta_si-theta_se)/np.log(D_e/D_i)*2*pi*lamb
#4.Thermal conductivity of insulation materials
    lamb = k_0 + (theta_si/2+theta_se/2)*k_1 + (theta_si/2+theta_se/2)**2*k_2
#5.Density of heat flow ratelinear --> Linear density of heat flow rate
    q_l=pi*D_e*q

"""

from scipy.optimize import root,fsolve
import numpy as np
import math
from matplotlib import pyplot as plt

#Setting
vertical_pipes = 1 #1 means the tube is vertical,0 means the tube is horizontal
D_i = 0.3 #Pipe diameter
theta_si = 400 #Medium temperature 
theta_a = 20 #Ambient temperature
q_max = 204 #Maximum allowable density of heat flow rate
theta_se_max = 45 #Maximum allowable surface temperature
epsilon = 0.25 #The blackness of galvanized steel plate, other protective materials refer to Table 2 of ISO 12241 2008
v = 0 #Wind speed
k_0 = 0.03 #0th-order coefficient of thermal conductivity of thermal insulation material,when the temperature is 0 ℃ (273.15 K), the thermal conductivity is k_0
k_1 = 5*10**-5 #1th-order coefficient of thermal conductivity of thermal insulation material
k_2 = 2*10**-7 #2th-order coefficient of thermal conductivity of thermal insulation material

#The heat transfer equation system is established, and theta_se, D_e, q_l can be solved by the q
def q_to_theta_se_and_D_e_and_q_l(x):
    theta_se,D_e,q_l,lamb = x[0],x[1],x[2],x[3]
    sigma = 5.67*10**-8 
    C_r = epsilon*sigma
    alpha_r = ((theta_se+273.15)**4 - (theta_a+273.15)**4)/((theta_se+273.15) - (theta_a+273.15))
    h_r = alpha_r*C_r
    if v == 0 :
        if vertical_pipes :
            if D_e**3*(theta_se - theta_a) <= 10 :
                h_cv = 1.32*((theta_se - theta_a)/D_e)**0.25
            else :
                h_cv = 1.74*(theta_se - theta_a)**(1/3)
        else :
            if D_e**3*(theta_se - theta_a) <= 10 :
                h_cv = 1.25*((theta_se - theta_a)/D_e)**0.25
            else :
                h_cv = 1.21*((theta_se - theta_a)/D_e)**(1/3)
    elif v*D_e <= 8.55*10**-3:
        h_cv = 8.1*10**-3/D_e + 3.14*(v/D_e)**0.5
    else:
        h_cv = 8.9*v**0.9/D_e**0.1
    h_se = h_r + h_cv
    return np.array([(theta_se-theta_a)*math.pi*D_e*h_se-q_l,
                     (theta_si-theta_se)/np.log(D_e/D_i)*2*math.pi*lamb-q_l,
                     k_0 + (theta_si/2+theta_se/2)*k_1 + (theta_si/2+theta_se/2)**2*k_2-lamb,
                     math.pi*D_e*q-q_l
                     ])

#The heat transfer equation system is established, and q,D_e,q_l can be solved by the theta_se
def theta_se_to_q_and_D_e_and_q_l(x):
    q,D_e,q_l,lamb = x[0],x[1],x[2],x[3]
    sigma = 5.67*10**-8 
    C_r = epsilon*sigma
    alpha_r = ((theta_se+273.15)**4 - (theta_a+273.15)**4)/((theta_se+273.15) - (theta_a+273.15))
    h_r = alpha_r*C_r
    if v == 0 :
        if vertical_pipes :
            if D_e**3*(theta_se - theta_a) <= 10 :
                h_cv = 1.32*((theta_se - theta_a)/D_e)**0.25
            else :
                h_cv = 1.74*(theta_se - theta_a)**(1/3)
        else :
            if D_e**3*(theta_se - theta_a) <= 10 :
                h_cv = 1.25*((theta_se - theta_a)/D_e)**0.25
            else :
                h_cv = 1.21*((theta_se - theta_a)/D_e)**(1/3)
    elif v*D_e <= 8.55*10**-3:
        h_cv = 8.1*10**-3/D_e + 3.14*(v/D_e)**0.5
    else:
        h_cv = 8.9*v**0.9/D_e**0.1
    h_se = h_r + h_cv
    return np.array([(theta_se-theta_a)*math.pi*D_e*h_se-q_l,
                     (theta_si-theta_se)/np.log(D_e/D_i)*2*math.pi*lamb-q_l,
                     k_0 + (theta_si/2+theta_se/2)*k_1 + (theta_si/2+theta_se/2)**2*k_2-lamb,
                     math.pi*D_e*q-q_l
                     ])

#The heat transfer equation system is established, and theta_se,q,q_l can be solved by the D_e
def D_e_to_theta_se_and_q_and_q_l(x):
    theta_se,q,q_l,lamb = x[0],x[1],x[2],x[3]
    sigma = 5.67*10**-8 
    C_r = epsilon*sigma
    alpha_r = ((theta_se+273.15)**4 - (theta_a+273.15)**4)/((theta_se+273.15) - (theta_a+273.15))
    h_r = alpha_r*C_r
    if v == 0 :
        if vertical_pipes :
            if D_e**3*(theta_se - theta_a) <= 10 :
                h_cv = 1.32*((theta_se - theta_a)/D_e)**0.25
            else :
                h_cv = 1.74*(theta_se - theta_a)**(1/3)
        else :
            if D_e**3*(theta_se - theta_a) <= 10 :
                h_cv = 1.25*((theta_se - theta_a)/D_e)**0.25
            else :
                h_cv = 1.21*((theta_se - theta_a)/D_e)**(1/3)
    elif v*D_e <= 8.55*10**-3:
        h_cv = 8.1*10**-3/D_e + 3.14*(v/D_e)**0.5
    else:
        h_cv = 8.9*v**0.9/D_e**0.1
    h_se = h_r + h_cv
    return np.array([(theta_se-theta_a)*math.pi*D_e*h_se-q_l,
                     (theta_si-theta_se)/np.log(D_e/D_i)*2*math.pi*lamb-q_l,
                     k_0 + (theta_si/2+theta_se/2)*k_1 + (theta_si/2+theta_se/2)**2*k_2-lamb,
                     math.pi*D_e*q-q_l
                     ])

# return h_se
def out_h_se():
    sigma = 5.67*10**-8 
    C_r = epsilon*sigma
    alpha_r = ((theta_se+273.15)**4 - (theta_a+273.15)**4)/((theta_se+273.15) - (theta_a+273.15))
    h_r = alpha_r*C_r
    if v == 0 :
        if vertical_pipes :
            if D_e**3*(theta_se - theta_a) <= 10 :
                h_cv = 1.32*((theta_se - theta_a)/D_e)**0.25
            else :
                h_cv = 1.74*(theta_se - theta_a)**(1/3)
        else :
            if D_e**3*(theta_se - theta_a) <= 10 :
                h_cv = 1.25*((theta_se - theta_a)/D_e)**0.25
            else :
                h_cv = 1.21*((theta_se - theta_a)/D_e)**(1/3)
    elif v*D_e <= 8.55*10**-3:
        h_cv = 8.1*10**-3/D_e + 3.14*(v/D_e)**0.5
    else:
        h_cv = 8.9*v**0.9/D_e**0.1
    h_se = h_r + h_cv
    return h_se

#solve theta_se, D_e, q_l from the q_max
q = q_max  #
sol_root1 = root(q_to_theta_se_and_D_e_and_q_l,[50,D_i+0.20,math.pi*(D_i+0.20)*q,k_0]) 
sol_fsolve1 = fsolve(q_to_theta_se_and_D_e_and_q_l,[50,D_i+0.20,math.pi*(D_i+0.20)*q,k_0])
theta_se = sol_fsolve1[0]
D_e = sol_fsolve1[1]
q_l = sol_fsolve1[2]
lamb = sol_fsolve1[3]
print(" ")
print("solution :",sol_root1)
print(">>>>>>>>>>> When the q is ",q," :")
print("  theta_se   = ",theta_se)
print("         D_e = ",D_e)
print("         q_l = ",q_l)
print("        lamb = ",lamb)
print("        h_se = ",out_h_se())

#if theta_se > theta_se_max,solve q, D_e, q_l from the theta_se_max
if theta_se > theta_se_max:
    theta_se=theta_se_max
    sol_root2 = root(theta_se_to_q_and_D_e_and_q_l,[q,D_e,q_l,lamb]) 
    sol_fsolve2 = fsolve(theta_se_to_q_and_D_e_and_q_l,[q,D_e,q_l,lamb]) 
    q = sol_fsolve2[0]
    D_e = sol_fsolve2[1]
    q_l = sol_fsolve2[2]
    lamb = sol_fsolve2[3]
    print(" ")
    print(" ")
    print(" ")
    print("solution :",sol_root2)
    print(">>>>>>>>>>> When the theta_se is ",theta_se_max," :")
    print("    q = ",q)
    print("  D_e = ",D_e)
    print("  q_l = ",q_l)
    print("        lamb = ",lamb)
    print("        h_se = ",out_h_se())

#When the thickness of the felt is 1cm, the feasible actual insulation layer thickness
delta = round(D_e/2-D_i/2+0.005,2)
print(" ")
print(" ")
print(" ")
print(">>>>>>>>>>> When the thickness of the single-layer felt is 1cm , insulation thickness :")
print("  delta = ",delta)

#solve theta_se、q、q_l from the D_e
D_e = D_i + 2*delta
sol_root3 = root(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
sol_fsolve3 = fsolve(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
theta_se = sol_fsolve3[0]
q = sol_fsolve3[1]
q_l = sol_fsolve3[2]
lamb = sol_fsolve3[3]
print(" ")
print(" ")
print(" ")
print("solution :",sol_root3)
print(">>>>>>>>>>> When the delta is ",delta," :")
print("  delta   = ",delta)
print("theta_se  = ",theta_se)
print("        q = ",q)
print("      q_l = ",q_l)
print("        lamb = ",lamb)
print("        h_se = ",out_h_se())


#Display the density of heat flow rate and surface temperature of similar structural schemes in the form of a list
print(" ")
print(" ")
print(" ")
print(">>>>>>>>>>>Parameters of optimal insulation thickness(structure 3) and adjacent thickness :")
print("structure delta theta_se q q_l lamb h_se")
D_e = D_i + 2*(delta-0.02)
sol_root3 = root(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
sol_fsolve3 = fsolve(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
theta_se = sol_fsolve3[0]
q = sol_fsolve3[1]
q_l = sol_fsolve3[2]
lamb = sol_fsolve3[3]
print(1,round(delta-0.02,2),round(theta_se,2),round(q,2),round(q_l,2),round(lamb,4),round(out_h_se(),2))

D_e = D_i + 2*(delta-0.01)
sol_root3 = root(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
sol_fsolve3 = fsolve(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
theta_se = sol_fsolve3[0]
q = sol_fsolve3[1]
q_l = sol_fsolve3[2]
lamb = sol_fsolve3[3]
print(2,round(delta-0.01,2),round(theta_se,2),round(q,2),round(q_l,2),round(lamb,4),round(out_h_se(),2))

D_e = D_i + 2*delta
sol_root3 = root(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
sol_fsolve3 = fsolve(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
theta_se = sol_fsolve3[0]
q = sol_fsolve3[1]
q_l = sol_fsolve3[2]
lamb = sol_fsolve3[3]
print(3,round(delta,2),round(theta_se,2),round(q,2),round(q_l,2),round(lamb,4),round(out_h_se(),2))

D_e = D_i + 2*(delta+0.01)
sol_root3 = root(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
sol_fsolve3 = fsolve(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
theta_se = sol_fsolve3[0]
q = sol_fsolve3[1]
q_l = sol_fsolve3[2]
lamb = sol_fsolve3[3]
print(4,round(delta+0.01,2),round(theta_se,2),round(q,2),round(q_l,2),round(lamb,4),round(out_h_se(),2))

D_e = D_i + 2*(delta+0.02)
sol_root3 = root(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
sol_fsolve3 = fsolve(D_e_to_theta_se_and_q_and_q_l,[theta_se,q,q_l,lamb]) 
theta_se = sol_fsolve3[0]
q = sol_fsolve3[1]
q_l = sol_fsolve3[2]
lamb = sol_fsolve3[3]
print(5,round(delta+0.02,2),round(theta_se,2),round(q,2),round(q_l,2),round(lamb,4),round(out_h_se(),2))
