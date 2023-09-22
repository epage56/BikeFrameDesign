#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:28:51 2023

@author: elipage
"""


#%%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button

# Link lengths
a = 100   # mm
b = 250   # mm
c = 110   # mm
d = 244.244  # mm

def freudenstein(a, b, c, d, theta2start, theta2end):
    # Solution of the vector loop equation proven in paper and by hand - wilk picture of whiteboard
    
    theta2start = np.deg2rad(theta2start)
    theta2end = np.deg2rad(theta2end)
    
    theta2 = np.linspace(theta2start, theta2end, 100)
    
    K1 = d/a
    K2 = d/c
    K3 = (a**2 - b**2 + c**2 + d**2)/(2*a*c)

    A1 = np.cos(theta2) - K1 - K2*np.cos(theta2) + K3
    B1 = -2*np.sin(theta2)
    C1 = K1 - (K2+1)*np.cos(theta2) + K3
    
    theta4_1 = 2 * np.arctan2(-B1 + np.sqrt(B1**2 - 4*A1*C1), 2*A1) # this helps with quadrant specificity
    theta4_2 = 2 * np.arctan2(-B1 - np.sqrt(B1**2 - 4*A1*C1), 2*A1)

    # sol for theta3: (not used in the current calcs but good to have tbh)
    K4 = d/b 
    K5 = (c**2 - d**2 - a**2 - b**2)/(2*a*b)
    
    D1 = np.cos(theta2) - K1 + K4*np.cos(theta2) + K5
    E1 = -2*np.sin(theta2)
    F1 = K1 - (K4+1)*np.cos(theta2) + K5
    
    theta3_1 = 2 * np.arctan2(-E1 + np.sqrt(E1**2 - 4*D1*F1), 2*D1) # this helps with quadrant specificity
    theta3_2 = 2 * np.arctan2(-E1 - np.sqrt(E1**2 - 4*D1*F1), 2*D1)
    
    #outputs respond to the choice of +(1) or -(2) in the eqn.
    return theta4_1, theta4_2, theta2, theta3_1, theta3_2

def calc_joint_pos(inputangle, outputangle, a, b, c, d):

    xA = a * np.cos(inputangle) #theta2
    yA = a * np.sin(inputangle)

    xB = d + c * np.cos(outputangle) #theta4
    yB = 0 + c * np.sin(outputangle)
    
    x04 = 0
    y04 = d
    
    xO2 = 0
    yO2 = 0 
    
    return  xA, yA, xB, yB, x04, y04, xO2, yO2

def grashof(X1):
    # Create a list of tuples where each tuple contains the original value and its index
    indexed_values = list(enumerate(X1))
    
    # Sort the list of tuples based on the values (second element of each tuple)
    sorted_values = sorted(indexed_values, key=lambda x: x[1])
    
    # Extract the indexes from the sorted list
    indexes = [x[0] for x in sorted_values]

    return indexes, sorted_values 

#input is frequnecy of input crank, outputs time array, also theta2 array
def movin(frequency):
    f = frequency/60
    omega = 2 * np.pi * f
    t = np.linspace(0.6, 5, 100)
    y = a * np.sin(omega * t)
    x = a * np.cos(omega * t)
    theta2 = np.arctan2(y, x)
    print("HEREEEEEEEEEE")
    return theta2, t

def find_time(theta2, frequency):
    f = frequency / 60
    omega = 2 * np.pi * f
    a = 1  # Assuming a constant amplitude, you can change this as needed.
    
    # Rearrange the equation to find t
    t = np.arccos(np.cos(theta2) / a) / omega

    return t

def plot_output_paths(inputangletheta2, outputtheta4_1, outputtheta4_2, xB, yB, xdes, ydes):
    #plotting the results:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.subplots_adjust(left=0.1, right=0.9)

    #line_theta2, = ax1.plot(t, np.degrees(inputangletheta2) % 360, color='k', label=r'$θ_2$')
    #line_theta4_1, = ax1.plot(t, np.degrees(outputtheta4_1) % 360, color='b', label=r'$θ_{4_1}$')
    #line_theta4_2, = ax1.plot(t, np.degrees(outputtheta4_2) % 360, color='r', label=r'$θ_{4_2}$')
    ax1.set_xlabel('t [s]')
    ax1.set_ylabel('degrees')
    ax1.legend()

    
    ax2.plot(xdes, ydes, color='b', label='Desired Output Path')
    ax2.plot(xB, yB, color='k', label='Optimized Output Path')
    
    ax2.set_xlim(100, 500)
    ax2.set_ylim(300, 500)
    ax2.set_xlabel('Distance (mm)')
    ax2.set_ylabel('Distance (mm)')
    ax2.set_title('Output paths')
    ax2.legend()
    


# startinginput = 0
# endinginput = 10

# outputtheta4_1, outputtheta4_2, inputtheta2 = freudenstein(a, b, c, d, startinginput, endinginput)[:3]
# xA, yA, xB, yB = calc_joint_pos(inputtheta2, outputtheta4_1, a,b,c,d)[:4]

# plot_output_paths(inputtheta2, outputtheta4_1, outputtheta4_2, xB, yB)




#%% Other functions - nan remove and circle plotter 
def plot_semi_circle(center, radius, concave_up=True):
    # Generate an array of angles from 0 to pi if concave up or from pi to 0 if concave down
    if concave_up:
        theta = np.linspace(0, np.pi, 100)
    else:
        theta = np.linspace(np.pi, 0, 100)

    # Calculate x and y coordinates of points on the semi-circle
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)

    # Plot the semi-circle
    plt.figure(figsize=(6, 6))
    plt.plot(x, y)
    plt.axis('equal')  # Equal aspect ratio for a circular appearance
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Semi-circle Plot')
    plt.grid(True)
    plt.show()
    
    return x, y 
    

def create_concave_down_curve(vertex_x=300, width=200, plotter = True):
    a = -15 / ((width / 2) ** 2)
    x = np.linspace(vertex_x - width / 2, vertex_x + width / 2, 100)
    y = a * (x - vertex_x) ** 2 + 400

    if plotter: 
        plt.plot(x, y)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Concave Down Curve with Vertex at (450, 0)')
        plt.grid(True)
        plt.show()
        plt.xlim(100, 600)
        plt.ylim(-400, 500)
    
    return x, y

def remove_nan_inf_rows(input_array):
    # Identify rows with NaNs or infs
    nan_rows = np.isnan(input_array).any(axis=1)
    inf_rows = np.isinf(input_array).any(axis=1)

    # Combine the conditions to find rows with NaNs or infs
    bad_rows = nan_rows | inf_rows

    # Remove rows with NaNs or infs
    cleaned_array = input_array[~bad_rows]

    return cleaned_array


create_concave_down_curve(plotter = True)


#%%
import numpy as np
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.algorithms.soo.nonconvex.pso import PSO
from pymoo.termination.default import DefaultSingleObjectiveTermination
from pymoo.termination import get_termination


def objective_function(X1): 
    # Unpack parameters
    a, b, c, d , theta2start, theta2end = X1
    
    theta4_1, theta4_2 = freudenstein(a, b, c, d, theta2start, theta2end)[:2]
    #Get output path for set of parameters 
    output_path_xB, output_path_yB = calc_joint_pos(theta4_1, theta4_2, a, b, c, d) [2:4]
    
    target_curve_ax_x, target_curve_ax_y = create_concave_down_curve(plotter = False)
    
    outputpath2d = np.column_stack((output_path_xB, output_path_yB))
    idealoutputpath2d = np.column_stack((target_curve_ax_x, target_curve_ax_y))
    
    cleaned_ideal_ax_2d = remove_nan_inf_rows(idealoutputpath2d)
    cleaned_exp_ax_2d = remove_nan_inf_rows(outputpath2d)
    
    # Calculate the MSE between the linkage path and target curve
    mse = np.mean((output_path_xB - target_curve_ax_x) ** 2 + (output_path_yB - target_curve_ax_y) ** 2)
     
    try:
        dtw, path = fastdtw(cleaned_exp_ax_2d, cleaned_ideal_ax_2d, dist = euclidean)
    except:
        dtw = 100000 # i think the alg sorts this out when it begins 
        
    axel_path_score = dtw
    
    return axel_path_score

def objectivefunction2(X1):
    
    """
    define leverage ratio real and compare to leverage ratio curve ideal
    """
    
    a, b, c, d , theta2start, theta2end = X1
    
    leverage_score = a*b+c - np.sin(theta2start)
    
    return leverage_score

def travel_constraint_function(X1): 
    a, b, c, d , theta2start, theta2end = X1
    
    theta4_1, theta4_2, theta2 = freudenstein(a, b, c, d, theta2start, theta2end)[:3]
    output_path_xB, output_path_yB = calc_joint_pos(theta2, theta4_1, a, b, c, d) [2:4]
    #getting output path to do constraint calcs
    
    initial_disp = output_path_xB[1] - output_path_xB[0]
    
    #needs to be positive 
    
    continuous_path = np.isnan(output_path_yB).any()
    continuous_path = int(continuous_path)
    #needs to be True
    
    displacement_path_x = abs(output_path_xB[99] - output_path_xB[0])
    #needs to be more than 195, less than 205
    
    displacement_path_x = abs(output_path_xB[99] - output_path_xB[0])
    if output_path_xB.all() > 0 and output_path_yB.all() > 0:
    # The point is in the positive quadrant
        correct_quadrant = True
    else:
        correct_quadrant = False
        
    correct_quadrant = int(correct_quadrant)
    
    return continuous_path, correct_quadrant, displacement_path_x, initial_disp


class LinkageOptimizationProblem(ElementwiseProblem):
    def __init__(self):
        # Define the bounds for your design variables (a, b, c, d)
        xl = np.array([210, 390, -440, 200, 60, 100])
        xu = np.array([250, 420, -400, 430, 90, 170])
        super().__init__(n_var=6, n_obj=1, n_ieq_constr=5, xl=xl, xu=xu)

    def _evaluate(self, X1, out, *args, **kwargs):

        f1 = objective_function(X1)
        #f2 = objective_function(X2)
    
        g1 = + travel_constraint_function(X1)[0] - 0.5
        g2 = - travel_constraint_function(X1)[1] + 0.5
        g3 = + travel_constraint_function(X1)[2] - 202
        g4 = - travel_constraint_function(X1)[2] + 198
        g5 = + travel_constraint_function(X1)[3]
        
        out["F"] = [f1]
        out["G"] = np.column_stack([g1, g2, g3, g4, g5])
        
        #grashof condition? 
        # g6 = - X1[grashof(X1)[0][0]] - X1[grashof(X1)[0][3]] + X1[grashof(X1)[0][2]] + X1[grashof(X1)[0][1]]

problem = LinkageOptimizationProblem()

termination = get_termination("n_gen", 1000)
algorithm = PSO()

# #Create a GA algorithm instance
# algorithm = GA(
#     pop_size=500,
#     eliminate_duplicates=True,
# )

initial_guess = np.array([200, 400, -500, 310, 75, 150])

# Perform optimization
res = minimize(problem, algorithm, seed=1, verbose = True, x0=initial_guess) 

# Print the best solution found
print("Best solution found:\nX = %s\nF = %s" % (res.X, res.F))


a, b, c, d, theta2start, theta2end =res.X

outputtheta4_1, outputtheta4_2, inputtheta2 = freudenstein(a, b, c, d, theta2start, theta2end)[:3]
xA, yA, xB, yB = calc_joint_pos(inputtheta2, outputtheta4_1, a,b,c,d)[:4]
plot_output_paths(inputtheta2, outputtheta4_1, outputtheta4_2, xB, yB)

#%%


def travel_constraint_function(X1): 
    a, b, c, d , theta2start, theta2end = X1
    
    theta4_1, theta4_2, theta2 = freudenstein(a, b, c, d, theta2start, theta2end)[:3]
    output_path_xB, output_path_yB = calc_joint_pos(theta2, theta4_1, a, b, c, d) [2:4]
    #getting output path to do constraint calcs
    
    initial_disp = output_path_xB[1] - output_path_xB[0]
    
    #needs to be positive 
    
    continuous_path = np.isnan(output_path_yB).any()
    continuous_path = int(continuous_path)
    #needs to be True
    
    displacement_path = abs(output_path_xB[-1] - output_path_xB[0])
    
    #needs to be more than 195, less than 205
    
    if output_path_xB.all() > 0 and output_path_yB.all() > 0:
    # The point is in the positive quadrant
        correct_quadrant = True
    else:
        correct_quadrant = False
        
    correct_quadrant = int(correct_quadrant)
    
    return continuous_path, correct_quadrant, displacement_path, initial_disp


test = [250. ,390. ,-400. , 530.  , 77  ,121.]

a, b, c, d, theta2start, theta2end = test

xdes, ydes = create_concave_down_curve()

outputtheta4_1, outputtheta4_2, inputtheta2 = freudenstein(a, b, c, d, theta2start, theta2end)[:3]

xA, yA, xB, yB = calc_joint_pos(inputtheta2, outputtheta4_1, a,b,c,d)[:4]

print(xB)

print(calc_joint_pos(inputtheta2, outputtheta4_1, a,b,c,d)[2:4])

print(xB[-1] - xB[0])

plot_output_paths(inputtheta2, outputtheta4_1, outputtheta4_2, xB, yB, xdes, ydes)

# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
# fig.subplots_adjust(left=0.1, right=0.9)

# #line_theta2, = ax1.plot(t, np.degrees(inputangletheta2) % 360, color='k', label=r'$θ_2$')
# #line_theta4_1, = ax1.plot(t, np.degrees(outputtheta4_1) % 360, color='b', label=r'$θ_{4_1}$')
# line_theta4_2, = ax1.plot(np.degrees(inputtheta2), np.degrees(outputtheta4_2) % 360, color='r', label=r'$θ_{4_2}$')
# ax1.set_xlabel('t [s]')
# ax1.set_ylabel('degrees')
# ax1.legend()

# ax2.plot(xB, yB, color='k', label='Optimized Output Path')
# ax2.plot(xB, yB, color='k', label='Desired Output Path')

# ax2.set_xlabel('Distance (mm) x')
# ax2.set_ylabel('Distance (mm) y')
# ax2.set_title('Output paths')
# ax2.legend()

travel_constraint_function(test)

