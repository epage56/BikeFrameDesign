#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:28:51 2023

@author: elipage
"""

import numpy as np
import matplotlib.pyplot as plt

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

def plot_output_paths(inputangletheta2, outputtheta4_1, outputtheta4_2, xB, yB, xdes, ydes, levydes, levy):
    #plotting the results:
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.9)
    
    ax1.plot(xdes, ydes, color='b', label='Desired Output Path')
    ax1.plot(xB, yB, color='k', label='Optimized Output Path')
    
    ax1.set_xlim(100, 500)
    ax1.set_ylim(300, 500)
    ax1.set_xlabel('Distance (mm)')
    ax1.set_ylabel('Distance (mm)')
    ax1.set_title('Axel paths')
    ax1.legend()
    
    ax2.plot(xB, levydes, color='b', label='Desired Leverage Ratio')
    ax2.plot(xB, levy, color='k', label='Optimized Leverage Ratio')
    
    ax2.set_xlim(0, 230)
    ax2.set_ylim(0, 4)
    ax2.set_xlabel('Travel (mm)')
    ax2.set_ylabel('Leverage Rate')
    ax2.set_title('Leverage Ratio')
    ax2.legend()

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

def create_ideal_leverage_ratio(start_point=(0, 3.5), end_point=(200, 2.4), num_points=100, curviness = -0.01, plotter=True):
    x = np.linspace(start_point[0], end_point[0], num_points)
    b = start_point[1]
    a = curviness
    y = b * np.exp(a * x)
    
    if plotter:
        plt.plot(x, y)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Curved Exponential Curve')
        plt.grid(True)
        plt.show()

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

def calculate_differences(arr):
    differences = []
    for i in range(1, len(arr)):
        diff = arr[i] - arr[i - 1]
        differences.append(diff)
    return differences



#%%
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.algorithms.soo.nonconvex.pso import PSO
from pymoo.termination.default import DefaultSingleObjectiveTermination
from pymoo.termination import get_termination
from pymoo.algorithms.moo.nsga2 import NSGA2

import similaritymeasures

"""
constrain the min and max init starting points for the curve
implement more constraints, and even when doing the leverage ratio
 part constrain the init and final leverage curves.
 
this will help with the 6bar

"""

def objective_function(X1): 
    # Unpack parameters
    a, b, c, d , theta2start, theta2end = X1
    
    theta4_1, theta4_2 = freudenstein(a, b, c, d, theta2start, theta2end)[:2]
    #Get output path for set of parameters 
    output_path_xB, output_path_yB = calc_joint_pos(theta4_1, theta4_2, a, b, c, d) [2:4]
    
    target_curve_ax_x, target_curve_ax_y = create_concave_down_curve(plotter = False)
    
    outputpath2d = np.column_stack((output_path_xB, output_path_yB))
    idealoutputpath2d = np.column_stack((target_curve_ax_x, target_curve_ax_y))
     
    frechet = similaritymeasures.frechet_dist(outputpath2d, idealoutputpath2d)
 
    axel_path_score = frechet 
    
    return axel_path_score 

def objective_function2(X1):
    
    """
    define leverage ratio real and compare to leverage ratio curve ideal
    """
    a, b, c, d , theta2start, theta2end = X1
    
    theta4_1, theta4_2 = freudenstein(a, b, c, d, theta2start, theta2end)[:2]
    #Get output path for set of parameters 
    output_path_xB, output_path_yB = calc_joint_pos(theta4_1, theta4_2, a, b, c, d) [2:4]
    
    target_curve_ax_x, target_curve_ax_y = create_concave_down_curve(plotter = False)
    
    outputpath2d = np.column_stack((output_path_xB, output_path_yB))
    idealoutputpath2d = np.column_stack((target_curve_ax_x, target_curve_ax_y))

    area = similaritymeasures.area_between_two_curves(outputpath2d, idealoutputpath2d)
    
    #pcm = similaritymeasures.pcm(outputpath2d, idealoutputpath2d)
    
    return area

def objective_function3(X1):
    
    """
    define leverage ratio real and compare to leverage ratio curve ideal
    """
    a, b, c, d , theta2start, theta2end = X1
    
    theta4_1, theta4_2 = freudenstein(a, b, c, d, theta2start, theta2end)[:2]
    #Get output path for set of parameters 
    output_path_xB, output_path_yB = calc_joint_pos(theta4_1, theta4_2, a, b, c, d) [2:4]
    
    target_curve_ax_x, target_curve_ax_y = create_concave_down_curve(plotter = False)
    
    outputpath2d = np.column_stack((output_path_xB, output_path_yB))
    idealoutputpath2d = np.column_stack((target_curve_ax_x, target_curve_ax_y))

    #area = similaritymeasures.area_between_two_curves(outputpath2d, idealoutputpath2d)
    
    pcm = similaritymeasures.pcm(outputpath2d, idealoutputpath2d)
    
    return pcm

def leverage_ratio(X1):
    
    """
    define leverage ratio real and compare to leverage ratio curve ideal
    """
    a, b, c, d , theta2start, theta2end = X1
    
    theta4_1, theta4_2, theta2= freudenstein(a, b, c, d, theta2start, theta2end)[:3]
    #Get output path for set of parameters 
    output_path_xB, output_path_yB = calc_joint_pos(theta4_1, theta4_2, a, b, c, d) [2:4]
    
    #clever reversal of direction to get shock link pos, also
    #remember that the shock pivot is at a certain angle from the pivot
    
    uppershockx = np.cos(theta2 - np.pi - np.deg2rad(30))
    uppershocky = np.sin(theta2 - np.pi - np.deg2rad(30)) 
    
    lowershockx = 320 
    lowershocky = -30
    
    inst_shock_length = np.sqrt((uppershockx - lowershockx)**2 + (uppershocky-lowershocky)**2)
    
    delta_shocktravl = calculate_differences(inst_shock_length)
    delta_axeltravl = calculate_differences(output_path_yB)
    
    #print(type(delta_axeltravl))
    
    delta_axeltravl = np.array(delta_axeltravl)
    delta_shocktravl = np.array(delta_shocktravl)
    
    #print(type(delta_axeltravl))

    delta_axeltravl = np.append(delta_axeltravl, delta_axeltravl[-1])
    delta_shocktravl = np.append(delta_shocktravl, delta_shocktravl[-1])
    
    #print(type(delta_axeltravl))
    #print(delta_axeltravl)

    #print(len(delta_axeltravl))
    
    leverage_ratio = delta_shocktravl / delta_axeltravl
    
    levratex, levratey = create_ideal_leverage_ratio(plotter = False)
    
    ideal_leverage_ratio_2d = np.column_stack((levratex, levratey))
    act_leverage_ratio_2d = np.column_stack((output_path_xB, leverage_ratio))
    
    pcm = similaritymeasures.pcm(act_leverage_ratio_2d, ideal_leverage_ratio_2d) 
    #area = similaritymeasures.area_between_two_curves(act_leverage_ratio_2d, ideal_leverage_ratio_2d)
    frechet = similaritymeasures.frechet_dist(act_leverage_ratio_2d, ideal_leverage_ratio_2d)

    similarity_meas_lev = pcm + frechet
    
    return similarity_meas_lev   

def travel_constraint_function(X1): 
    a, b, c, d , theta2start, theta2end = X1
    
    theta4_1, theta4_2, theta2 = freudenstein(a, b, c, d, theta2start, theta2end)[:3]
    output_path_xB, output_path_yB = calc_joint_pos(theta2, theta4_1, a, b, c, d) [2:4]
    #getting output path to do constraint calcs
    
    initposx = output_path_xB[0]
    initposy = output_path_yB[0]
    
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
    
    
    return continuous_path, correct_quadrant, displacement_path_x, initial_disp, initposx, initposy


class LinkageOptimizationProblem(ElementwiseProblem):
    def __init__(self):
        # Define the bounds for your design variables (a, b, c, d)
        xl = np.array([0, 0, -410, 290, 0, 90])
        xu = np.array([500, 500, -390, 310, 70, 160])
        super().__init__(n_var=6, n_obj=4, n_ieq_constr=7, xl=xl, xu=xu)

    def _evaluate(self, X1, out, *args, **kwargs):

        f1 = objective_function(X1)
        f2 = objective_function2(X1)
        f3 = objective_function3(X1)
        f4 = leverage_ratio(X1)

    
        g1 = + travel_constraint_function(X1)[0] - 0.5
        g2 = - travel_constraint_function(X1)[1] + 0.5
        g3 = + travel_constraint_function(X1)[2] - 202
        g4 = - travel_constraint_function(X1)[2] + 198
        g5 = + travel_constraint_function(X1)[3]
        g6 = - travel_constraint_function(X1)[4] + 399
        g7 = + travel_constraint_function(X1)[4] - 401
        
        out["F"] = [f1, f2, f3,f4]
        out["G"] = np.column_stack([g1, g2, g3, g4, g5, g6, g7])
        
problem = LinkageOptimizationProblem()

from pymoo.algorithms.moo.sms import SMSEMOA

initial_guess = np.array([200, 400, -500, 310, 75, 150])
algorithm = SMSEMOA()
res = minimize(problem,
               algorithm,
               ('n_gen', 500),
               seed=1,
               verbose=True)


# Print the best solution found
print("Best solution found:\nX = %s\nF = %s" % (res.X, res.F))

# xdes, ydes = create_concave_down_curve()
# a, b, c, d, theta2start, theta2end =res.X
# outputtheta4_1, outputtheta4_2, inputtheta2 = freudenstein(a, b, c, d, theta2start, theta2end)[:3]
# xA, yA, xB, yB = calc_joint_pos(inputtheta2, outputtheta4_1, a,b,c,d)[:4]
# plot_output_paths(inputtheta2, outputtheta4_1, outputtheta4_2, xB, yB, xdes, ydes)


#%% plots optimized path vs desired path


test = [ 443.03227881 , 249.36197987, -402.2258077 ,  290.40725884 ,  69.81842202
    ,95.6788756 ]

test = [ 378.04520179,  246.14044847, -405.60990119 ,290.17168419 ,  65.24948587,
    97.1609317 ]


a, b, c, d, theta2start, theta2end = test
# 
xdes, ydes = create_concave_down_curve()

outputtheta4_1, outputtheta4_2, inputtheta2 = freudenstein(a, b, c, d, theta2start, theta2end)[:3]

xA, yA, xB, yB = calc_joint_pos(inputtheta2, outputtheta4_1, a,b,c,d)[:4]

def alt_leverage_ratio(X1):
    
    """
    define leverage ratio real and compare to leverage ratio curve ideal
    """
    a, b, c, d , theta2start, theta2end = X1
    
    theta4_1, theta4_2, theta2= freudenstein(a, b, c, d, theta2start, theta2end)[:3]
    #Get output path for set of parameters 
    output_path_xB, output_path_yB = calc_joint_pos(theta4_1, theta4_2, a, b, c, d) [2:4]
    
    #clever reversal of direction to get shock link pos, also
    #remember that the shock pivot is at a certain angle from the pivot
    
    uppershockx = np.cos(theta2 - np.pi - np.deg2rad(30))
    uppershocky = np.sin(theta2 - np.pi - np.deg2rad(30)) 
    
    lowershockx = 320 
    lowershocky = -30
    
    inst_shock_length = np.sqrt((uppershockx - lowershockx)**2 + (uppershocky-lowershocky)**2)
    
    delta_shocktravl = calculate_differences(inst_shock_length)
    delta_axeltravl = calculate_differences(output_path_yB)
    
    #print(type(delta_axeltravl))
    
    delta_axeltravl = np.array(delta_axeltravl)
    delta_shocktravl = np.array(delta_shocktravl)
    
    #print(type(delta_axeltravl))

    delta_axeltravl = np.append(delta_axeltravl, delta_axeltravl[-1])
    delta_shocktravl = np.append(delta_shocktravl, delta_shocktravl[-1])
    
    #print(type(delta_axeltravl))
    #print(delta_axeltravl)

    #print(len(delta_axeltravl))
    
    leverage_ratio = delta_shocktravl / delta_axeltravl
    
    levratex, levratey = create_ideal_leverage_ratio(plotter = False)
        
    ideal_leverage_ratio_2d = np.column_stack((levratex, levratey))
    act_leverage_ratio_2d = np.column_stack((output_path_xB, leverage_ratio))
    
    pcm = similaritymeasures.pcm(act_leverage_ratio_2d, ideal_leverage_ratio_2d) 
    #area = similaritymeasures.area_between_two_curves(act_leverage_ratio_2d, ideal_leverage_ratio_2d)
    frechet = similaritymeasures.frechet_dist(act_leverage_ratio_2d, ideal_leverage_ratio_2d)

    similarity_meas_lev = pcm + frechet
    
    ideal_leverage_ratio = levratey
    real_leverage_ratio = leverage_ratio

    
    return similarity_meas_lev, ideal_leverage_ratio, real_leverage_ratio


junk, idealylev, optylev = alt_leverage_ratio(test)

xB1 = xB[::-1] 
xB1= xB1 - xB1[1]

print(xB1)

plot_output_paths(inputtheta2, outputtheta4_1, outputtheta4_2, xB, yB, xdes, ydes, idealylev, optylev)

#%%
plt.plot(xB1, optylev)
plt.xlabel('xB1')
plt.ylabel('optylev')
plt.title('Plot of xB1 vs. optylev')
plt.grid(True)
plt.show()


#WE NEED TO PUT IN LEVERAGE RATIO UP TOP NEAR FREUDEINSTEIN / REFACTOR A LITTLE SO WE CAN CALL IT HERE SO IT CAN GRAPH


#%%Pareto front plot
from pymoo.algorithms.moo.sms import SMSEMOA
from pymoo.optimize import minimize
from pymoo.problems import get_problem
from pymoo.visualization.scatter import Scatter

plot = Scatter()
plot.add(problem.pareto_front(), plot_type="line", color="black", alpha=0.7)
plot.add(res.F, color="red")
plot.show()
#%%PCP visualization of scores, 1, 30


from pymoo.visualization.pcp import PCP


plot = PCP(title=("Run", {'pad': 30}),
           n_ticks=10,
           legend=(True, {'loc': "upper left"}),
           labels=["a", "b", "c", "d", r"Start $\theta_2$", r"End $\theta_2$"]
           )

plot.set_axis_style(color="grey", alpha=0.5)
plot.add(res.X, color="grey", alpha=0.3)
plot.add(res.X[1], linewidth=5, color="red")
plot.add(res.X[30], linewidth=5, color="blue")
plot.show()

#%%PCP visualization of lengths, 1,30

plot = PCP()
plot.set_axis_style(color="grey", alpha=0.5)
plot.add(res.F, color="grey", alpha=0.3)
plot.add(res.F[1], linewidth=5, color="red")
plot.add(res.F[30], linewidth=5, color="blue")
plot.show()

