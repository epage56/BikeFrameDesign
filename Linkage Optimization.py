#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:28:51 2023

@author: elipage
"""

#%% 

#Purpose of this code block is to establish ideal curves for leverage and axel path
#this is done by picking points(drawing) and fitting to the drawing via polynomials

# Import necessary libraries
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define a second-degree polynomial function
def second_degree_polynomial(x, a, b, c):
    return a * x**2 + b * x + c

# Define a fourth-degree polynomial function
def fourth_degree_polynomial(x, a, b, c, d, e):
    return a * x**4 + b * x**3 + c * x**2 + d * x + e

# Function to fit a quadratic equation and return it in scientific notation
def fit_leverage_curve(X, Y):
    # Fit the data to the second-degree polynomial function
    params, _ = curve_fit(second_degree_polynomial, X, Y)
    
    # Extract the coefficients (a, b, c) from the fitted parameters
    a, b, c = params
    
    # Create the second-degree polynomial equation in scientific notation
    equation = f'{a:.3e}x^2 + {b:.3e}x + {c:.3e}'
    
    return equation, params

def fit_ax_curve(X, Y):
    # Fit the data to the fourth-degree polynomial function
    params1, _ = curve_fit(fourth_degree_polynomial, X, Y)
    
    # Extract the coefficients (a, b, c, d, e) from the fitted parameters
    a, b, c, d, e = params1
    
    # Create the fourth-degree polynomial equation in scientific notation
    equation = f'{a:.3e}x^4 + {b:.3e}x^3 + {c:.3e}x^2 + {d:.3e}x + {e:.3e}'
    
    return equation, params1

# Sample data (replace these with your X and Y arrays)
Y_lev = np.linspace(3.3, 2.5, 9)
X_lev = np.array([0, 10, 20, 40, 60, 90, 120, 140, 200])

X_ax = np.linspace(0, 200, 14)
Y_ax = (-10, -5, 0, -0.5, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10)

# Call the function to fit the second-degree polynomial equation for leverage curve
quadratic_equation1 = fit_leverage_curve(X_lev, Y_lev)

# Call the function to fit the fourth-degree polynomial equation for axle path curve
quadratic_equation2 = fit_ax_curve(X_ax, Y_ax)

# Create a new X range for the fitted curve
X_lev_fit = np.linspace(min(X_lev), max(X_lev), 100)
X_ax_fit = np.linspace(min(X_ax), max(X_ax), 100)

# Calculate the corresponding Y values using the fitted coefficients
Y_lev_fit = second_degree_polynomial(X_lev_fit, quadratic_equation1[1][0], quadratic_equation1[1][1], quadratic_equation1[1][2])
Y_ax_fit = fourth_degree_polynomial(X_ax_fit, quadratic_equation2[1][0], quadratic_equation2[1][1], quadratic_equation2[1][2], quadratic_equation2[1][3], quadratic_equation2[1][4])

# Close all existing plots
plt.close("all")

# Create the first figure for the leverage curve
plt.figure()
plt.scatter(X_lev, Y_lev, label='Data')
plt.plot(X_lev_fit, Y_lev_fit, label='Fitted Curve', color='red')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.grid(True)
plt.title('Figure 1 - Fitted Leverage Curve')  # Add a title
plt.show()

# Create the second figure for the axle path curve
plt.figure()
plt.scatter(X_ax, Y_ax, label='Data')
plt.plot(X_ax_fit, Y_ax_fit, label='Fitted Curve', color='red')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.grid(True)
plt.title('Figure 2 - Fitted Axel Path Curve')  # Add a title
plt.show()

# Print the fitted polynomial equations in scientific notation
print(f'Fitted Leverage Equation: {quadratic_equation1[0]}')
print(f'Fitted Axel Path Equation: {quadratic_equation2[0]}')

#coeff arrays for the rest of the optimization problem

lev_fit_coeffs = quadratic_equation1[1]
ax_fit_coeffs = quadratic_equation2[1]


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

# test_result = np.array([50.12684447, 134.15106216, 488.36842375, 569.95798077])
# test_result = np.array([50.15336187, 118.72426894,  53.28384091, 120.27259326])
# test_result = np.array([50.32612371 ,101.16358628 , 52.01101179, 101.83780644])
# test_result = np.array([54.93845191 ,86.08845309 ,64.96986728 ,89.59000895])
# test_result = np.array([50.00183095,  72.20425701, 222.08995491, 199.88918179])
test_result = np.array([100, 200 ,200, 200])

a,b,c,d = test_result

# Input angular velocity
f = 16.7/60
omega = 2 * np.pi * f
t = np.linspace(0, 5, 100)
y = a * np.sin(omega * t)
x = a * np.cos(omega * t)
theta2 = np.arctan2(y, x)

fi2 = np.arctan2(y, x)

# Solution of the vector loop equation
K1 = d/a
K2 = d/c
K3 = (a**2 - b**2 + c**2 + d**2)/(2*a*c)
A = np.cos(fi2) - K1 - K2*np.cos(fi2) + K3
B = -2*np.sin(fi2)
C = K1 - (K2+1)*np.cos(fi2) + K3
fi4_1 = 2 * np.arctan2(-B + np.sqrt(B**2 - 4*A*C), 2*A)
fi4_2 = 2 * np.arctan2(-B - np.sqrt(B**2 - 4*A*C), 2*A)

# Create a figure with subplots for old and new graphs
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig.subplots_adjust(left=0.1, right=0.9)

# Initialize the old graph
theta4_1 = np.arccos((a**2 + b**2 - c**2 - d**2) / (2 * a * b))
theta4_2 = -theta4_1

line_theta2, = ax1.plot(t, np.degrees(fi2) % 360, color='k', label=r'$θ_2$')
line_theta4_1, = ax1.plot(t, np.degrees(fi4_1) % 360, color='b', label=r'$θ_{4_1}$')
line_theta4_2, = ax1.plot(t, np.degrees(fi4_2) % 360, color='r', label=r'$θ_{4_2}$')
ax1.set_xlabel('t [s]')
ax1.set_ylabel('degrees')
ax1.legend()

# Initialize the new animated graph
link2, = ax2.plot([], [], 'k-', linewidth=2, label='Link 2')
link3, = ax2.plot([], [], 'b-', linewidth=2, label='Link 3')
link4, = ax2.plot([], [], 'r-', linewidth=2, label='Link 4')

ax2.set_xlim(-400, 400)
ax2.set_ylim(-400, 400)
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_title('New Animated Graph')
ax2.legend()

# Initialize empty arrays for position history
x2_history = []
y2_history = []
x3_history = []
y3_history = []

# Animation function
def animate(i):
    theta2_i = theta2[i]
    fi4_i = fi4_1[i]

    x2 = a * np.cos(theta2_i)
    y2 = a * np.sin(theta2_i)

    # Calculate the positions of P3 based on theta2_i
    x3 = d + c * np.cos(fi4_i)
    y3 = 0 + c * np.sin(fi4_i)

    # Append current positions to history
    x2_history.append(x2)
    y2_history.append(y2)
    x3_history.append(x3)
    y3_history.append(y3)

    # Stationary position of the end of pivot 4
    x4_stationary = d
    y4_stationary = 0

    # Plot the history along with current positions
    link2.set_data([0, x2], [0, y2])
    link3.set_data([x2, x3], [y2, y3])
    link4.set_data([x3, x4_stationary], [y3, y4_stationary])

    # Plot the path history
    ax2.plot(x2_history, y2_history, 'k:', linewidth=0.5, alpha=0.5)
    ax2.plot(x3_history, y3_history, 'b:', linewidth=0.5, alpha=0.5)

# Create an animation
ani = FuncAnimation(fig, animate, frames=len(t), repeat=True, blit=False)

# Create a pause/play button
ax_pause_play = plt.axes([0.8, 0.01, 0.1, 0.05])
btn_pause_play = Button(ax_pause_play, 'Pause')

paused = False

def pause_play(event):
    global paused
    if paused:
        ani.event_source.start()
        btn_pause_play.label.set_text('Pause')
    else:
        ani.event_source.stop()
        btn_pause_play.label.set_text('Play')
    paused = not paused

btn_pause_play.on_clicked(pause_play)

plt.show()

#%%

plt.figure(figsize=(8, 4))  # Set the figure size
plt.plot(x111, y111, label='Semicircle', color='blue')  # Plot the data
plt.xlabel('x')  # Label for the x-axis
plt.ylabel('y')  # Label for the y-axis
plt.title('Semicircle Plot')  # Title of the plot
plt.grid(True)  # Display grid lines
plt.legend()  # Display legend
plt.axhline(0, color='black',linewidth=0.5)  # Add a horizontal line at y=0
plt.axvline(200, color='black',linewidth=0.5)  # Add a vertical line at x=200
plt.show()  # Show the plot


#%%
import numpy as np
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.problem import Problem

#### IDEA - FIT THE CURVE GIVEN TO A POLYNOMIAL REP THE output CURVE,
 ### WHICH WILL HAVE EVENLY DIST POINTS. bec u define the x range. TRY THIS NEXT TIME

plt.clf()

testarr = [100, 250, 110, 244.244]

# Define the target curve (replace with your own data)

x111 = np.linspace(-40, -20, 100)

y111 = -np.sqrt((25**2 - (x111 - 100)**2))

y111 = (1/2*x111+10)**2 - 100

plt.figure(figsize=(8, 4))  # Set the figure size
plt.plot(x111, y111, label='output curve to match', color='blue')  # Plot the data
plt.xlabel('x')  # Label for the x-axis
plt.ylabel('y')  # Label for the y-axis
plt.title('Semicircle Plot')  # Title of the plot
plt.grid(True)  # Display grid lines
plt.legend()  # Display legend
plt.axhline(0, color='black',linewidth=0.5)  # Add a horizontal line at y=0
plt.axvline(125, color='black',linewidth=0.5)  # Add a vertical line at x=200
plt.show()  # Show the plot


target_curve_ax_x = X_ax_fit
target_curve_ax_y = Y_ax_fit

target_curve_ax_x = x111
target_curve_ax_x = y111


# Define the objective function to minimize (e.g., MSE between linkage path and target curve) 
# passing x1 array (should be array of a,b,c,d) or linkage lengths
def objective_function(X1): 
    # Extract design variables (e.g., link lengths)
    a, b, c, d = X1

    # Simulate the four-bar linkage and obtain its path
    # You need to calculate the linkage path (x, y) based on the given link lengths and equations
    
    # Calculate the positions of P3 based on the four-bar linkage equations
    # Modify these equations according to your four-bar linkage model

    # Input angular velocity
    f = 16.7/60
    omega = 2 * np.pi * f
    t = np.linspace(0, np.pi, 100)
    y = a * np.sin(omega * t)
    x = a * np.cos(omega * t)
    theta2 = np.arctan2(y, x)

    fi2 = np.arctan2(y, x)

    # Solution of the vector loop equation
    K1 = d/a
    K2 = d/c
    K3 = (a**2 - b**2 + c**2 + d**2)/(2*a*c)
    A = np.cos(fi2) - K1 - K2*np.cos(fi2) + K3
    B = -2*np.sin(fi2)
    C = K1 - (K2+1)*np.cos(fi2) + K3
    fi4_1 = 2 * np.arctan2(-B + np.sqrt(B**2 - 4*A*C), 2*A)
    fi4_2 = 2 * np.arctan2(-B - np.sqrt(B**2 - 4*A*C), 2*A)
    
    theta4_1 = np.arccos((a**2 + b**2 - c**2 - d**2) / (2 * a * b))
    theta4_2 = -theta4_1
    
    # Calculate the x and y coordinates of the path points
    linkage_path_x = d + c * np.cos(fi4_1)
    linkage_path_y = c * np.sin(fi4_1)
    
    
    # Calculate the MSE between the linkage path and target curve
    mse = np.mean((linkage_path_x - target_curve_ax_x) ** 2 + (linkage_path_y - target_curve_ax_y) ** 2)
    
    #print("mse:", mse)
    
    return mse


import numpy as np
from pymoo.optimize import minimize

class LinkageOptimizationProblem(ElementwiseProblem):
    def __init__(self):
        # Define the bounds for your design variables (a, b, c, d)
        xl = np.array([20, 20, 20, 20])
        xu = np.array([500, 500, 500, 500])
        super().__init__(n_var=4, n_obj=1, n_ieq_constr=0, xl=xl, xu=xu)

    def _evaluate(self, X1, out, *args, **kwargs):
        # Call your objective_function with the design variables
        if len(X1) > 4:
            print(len(X1))
            print(len(target_curve_ax_x))
            print(len(target_curve_ax_y))
            
            raise ValueError("X should contain at most 4 values for the design variables.")
        
        f = objective_function(X1)
        
        # Assign the objective function value to the "F" key of the out dictionary
        out["F"] = [f]
        

problem = LinkageOptimizationProblem()
print("moo")

from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.algorithms.soo.nonconvex.pso import PSO

from pymoo.termination.default import DefaultSingleObjectiveTermination
from pymoo.termination import get_termination

#termination = get_termination("n_gen", 500)

# Create a GA algorithm instance
algorithm = GA(
    pop_size=5000,
    eliminate_duplicates=True,
)


initial_guess = np.array([100, 250, 110, 244.244])
# Perform optimization
res = minimize(problem, algorithm, seed=1, verbose = True, x0=initial_guess) 
#add termination to customize the termination = termination gen amnt


# Print the best solution found
print("Best solution found:\nX = %s\nF = %s" % (res.X, res.F))


test_result = res.X

a,b,c,d = test_result

# Input angular velocity
f = 16.7/60
omega = 2 * np.pi * f
t = np.linspace(0, np.pi, 50)
y = a * np.sin(omega * t)
x = a * np.cos(omega * t)
theta2 = np.arctan2(y, x)

fi2 = np.arctan2(y, x)

# Solution of the vector loop equation
K1 = d/a
K2 = d/c
K3 = (a**2 - b**2 + c**2 + d**2)/(2*a*c)
A = np.cos(fi2) - K1 - K2*np.cos(fi2) + K3
B = -2*np.sin(fi2)
C = K1 - (K2+1)*np.cos(fi2) + K3
fi4_1 = 2 * np.arctan2(-B + np.sqrt(B**2 - 4*A*C), 2*A)
fi4_2 = 2 * np.arctan2(-B - np.sqrt(B**2 - 4*A*C), 2*A)

# Create a figure with subplots for old and new graphs
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig.subplots_adjust(left=0.1, right=0.9)

# Initialize the old graph
theta4_1 = np.arccos((a**2 + b**2 - c**2 - d**2) / (2 * a * b))
theta4_2 = -theta4_1

line_theta2, = ax1.plot(t, np.degrees(fi2) % 360, color='k', label=r'$θ_2$')
line_theta4_1, = ax1.plot(t, np.degrees(fi4_1) % 360, color='b', label=r'$θ_{4_1}$')
line_theta4_2, = ax1.plot(t, np.degrees(fi4_2) % 360, color='r', label=r'$θ_{4_2}$')
ax1.set_xlabel('t [s]')
ax1.set_ylabel('degrees')
ax1.legend()

# Initialize the new animated graph
link2, = ax2.plot([], [], 'k-', linewidth=2, label='Link 2')
link3, = ax2.plot([], [], 'b-', linewidth=2, label='Link 3')
link4, = ax2.plot([], [], 'r-', linewidth=2, label='Link 4')

ax2.set_xlim(-400, 400)
ax2.set_ylim(-400, 400)
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_title('New Animated Graph')
ax2.legend()

# Initialize empty arrays for position history
x2_history = []
y2_history = []
x3_history = []
y3_history = []

# Animation function
def animate(i):
    theta2_i = theta2[i]
    fi4_i = fi4_1[i]

    x2 = a * np.cos(theta2_i)
    y2 = a * np.sin(theta2_i)

    # Calculate the positions of P3 based on theta2_i
    x3 = d + c * np.cos(fi4_i)
    y3 = 0 + c * np.sin(fi4_i)

    # Append current positions to history
    x2_history.append(x2)
    y2_history.append(y2)
    x3_history.append(x3)
    y3_history.append(y3)

    # Stationary position of the end of pivot 4
    x4_stationary = d
    y4_stationary = 0

    # Plot the history along with current positions
    link2.set_data([0, x2], [0, y2])
    link3.set_data([x2, x3], [y2, y3])
    link4.set_data([x3, x4_stationary], [y3, y4_stationary])

    # Plot the path history
    ax2.plot(x2_history, y2_history, 'k:', linewidth=0.5, alpha=0.5)
    ax2.plot(x3_history, y3_history, 'b:', linewidth=0.5, alpha=0.5)

# Create an animation
ani = FuncAnimation(fig, animate, frames=len(t), repeat=True, blit=False)

# Create a pause/play button
ax_pause_play = plt.axes([0.8, 0.01, 0.1, 0.05])
btn_pause_play = Button(ax_pause_play, 'Pause')

paused = False

def pause_play(event):
    global paused
    if paused:
        ani.event_source.start()
        btn_pause_play.label.set_text('Pause')
    else:
        ani.event_source.stop()
        btn_pause_play.label.set_text('Play')
    paused = not paused

btn_pause_play.on_clicked(pause_play)

plt.show()
#need to add something that makes sure the sol is possilbe 
#grashof cond?



    
#%%

"""
1. Link Lengths: The lengths of the four links 
(e.g., lengths of the crank, coupler, rocker, and ground links).

2. Joint Positions: The positions of the revolute joints where the links are connected 
(e.g., x and y coordinates of the joints).

3. Joint Angles: The angles at which the links are connected 
(e.g., angles of the crank, coupler, and rocker with respect to a reference axis).

4. Geometric Parameters: Other geometric parameters that define the shape
and orientation of the linkage.

5. Objective Variables: Variables that represent the performance measures 
you want to optimize (e.g., end-effector position, velocity, or torque).


I think the number of variables I need is large - 
    4 linkage lengths 
    1 Length from rocker pivot to shock
    2 points (x,y) in space determining stationary pivots 
        origin and x1,y1 i guess
    1 input crank angle determining initial 
    
    error from pre-defined axel path(polynomial rep by ax_fit_coeffs)
    error from pre-defined leverage curve(polynomial rep by lev_fit_coeffs)
    
    constrained by vector loop equation for planar linkage
    constrained by total y movement of rear pivot = 200 (-+1)
    constrained by grasshaulf condition I think
    constrained by initial distance from origin to axel pivot = 475 
    constrained by the init height diff (found by vector add) of the rear pivot and the origin
    constrained by an overall leverage ratio range [3.5, 2.1] 
    think about linkage size lower bounds in relation to stress and bolt sizes
    
    
    thinking of adding bb in here maybe and constraining off that? probably not
    
    
"""

import numpy as np
from pymoo.core.problem import ElementwiseProblem

class LinkageOptimizationProblem(ElementwiseProblem):
    
    def __init__(self, lev_fit_coeffs, ax_fit_coeffs):
        # Define the number of decision variables (DV)
        n_vars = 9  # 4 linkage lengths + 5 pivot points (x and y coordinates) + 1 crank angle + 1 angle of fixed pivot
        
        # Define the number of objectives (OBJ) and constraints (CON)
        n_obj = 2  # Two objectives: Error from axel path, Error from leverage curve
        n_con = 6  # Six constraints: Vector loop, y-movement, Grasshault, distance, pivot placement
        
        # Define lower and upper bounds for decision variables
        xl = np.array([0, 0, 0, 0, -2, -2, -2, -2, -360])  # Lower bounds for DVs
        xu = np.array([10, 10, 10, 10, 2, 2, 2, 2, 360])  # Upper bounds for DVs
    
        # Call the constructor of the parent class (ElementwiseProblem)
        super().__init__(n_var=n_vars, n_obj=n_obj, n_con=n_con, xl=xl, xu=xu)
        
        # Store the coefficients for later use
        self.lev_fit_coeffs = lev_fit_coeffs
        self.ax_fit_coeffs = ax_fit_coeffs

    def _evaluate(self, x, out, *args, **kwargs):
        # Extract decision variables
        linkage_lengths = x[:4]
        pivot_points = x[4:7]
        crank_angle = x[8]
        fixed_pivot_angle = x[9]
        
        # Implement your optimization objectives and constraints here
        # You can use the provided coefficients (lev_fit_coeffs, ax_fit_coeffs) and the decision variables (linkage_lengths, pivot_points, crank_angle) to calculate the objectives and constraints
        
        # Calculate objective 1: Error from axel path using ax_fit_coeffs
        error_axel_path = sum((y - ax_fit_coeffs[0] * x**2 - ax_fit_coeffs[1] * x - ax_fit_coeffs[2])**2 for x, y in axel_path_points)
        
        # Calculate objective 2: Error from leverage curve using lev_fit_coeffs
        error_leverage_curve = sum((Y - lev_fit_coeffs[0] * X**2 - lev_fit_coeffs[1] * X - lev_fit_coeffs[2])**2 for X, Y in zip(X_lev, Y_lev))
        
        # Calculate constraint 1: Vector loop equation
        constraint_vector_loop = ...
        
        # Calculate constraint 2: Y-movement constraint
        constraint_travel = "200mm"
        
        # Calculate constraint 3: grashof's law
        constraint_grashof = ...
        
        # Calculate constraint 4: Distance from origin to rear pivot constraint
        constraint_chainstay = ...
        
        # Calculate constraint 5: Pivot placement constraint
        constraint_no_sag_height_diff = ...
        
        # Calculate constraint 6: link size constraint (cant be too small)
        constraint_linksize = ...
        
        # Set the objectives and constraints
        out["F"] = [error_axel_path, error_leverage_curve]
        out["G"] = [constraint_vector_loop, constraint_y_movement, constraint_grasshault, constraint_distance, constraint_pivot_placement, fixed_pivot_angle]

# Create an instance of the problem class with the appropriate coefficients
problem = LinkageOptimizationProblem(lev_fit_coeffs, ax_fit_coeffs)

#%%
from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt

# Input
#lengths of links (tube testing machine actual lengths)
a = 100   #mm
b = 250   #mm
c = 110   #mm
d = 244.244  #mm

# Solution for fi2 being a time function, f(time) = angle
f = 16.7/60    #/s
omega = 2 * np.pi * f   #rad/s

t = np.linspace(0, 5, 100)
y = a * np.sin(omega * t)
x = a * np.cos(omega * t)

fi2 = np.arctan2(y, x)

# Solution of the vector loop equation
#https://scholar.cu.edu.eg/?q=anis/files/week04-mdp206-position_analysis-draft.pdf
K1 = d/a
K2 = d/c
K3 = (a**2 - b**2 + c**2 + d**2)/(2*a*c)
A = np.cos(fi2) - K1 - K2*np.cos(fi2) + K3
B = -2*np.sin(fi2)
C = K1 - (K2+1)*np.cos(fi2) + K3
fi4_1 = 2 * np.arctan2(-B + np.sqrt(B**2 - 4*A*C), 2*A)
fi4_2 = 2 * np.arctan2(-B - np.sqrt(B**2 - 4*A*C), 2*A)

# Plot the fi2 time diagram and fi4 time diagram
plt.clf()
plt.plot(t, np.degrees(fi2) % 360, color = 'k', label=r'$θ_2$')
plt.plot(t, np.degrees(fi4_1) % 360, color = 'b', label=r'$θ_{4_1}$')
plt.plot(t, np.degrees(fi4_2) % 360, color = 'r', label=r'$θ_{4_2}$')
plt.xlabel('t [s]')
plt.ylabel('degrees')
plt.legend()
plt.show()

#%%

#%%
#THIS DOESNT WORK

import math
import matplotlib.pyplot as plt
import numpy as np
import time

# Define parameters of the four-bar linkage
a = 5.0  # Length of the first link
b = 8.0  # Length of the second link
c = 7.0  # Length of the third link
d = 6.0  # Length of the fourth link

# Define an array of input link angles (in degrees)
input_link_angles_degrees = np.arange(0, 360, 10)

# Initialize lists to store output link angles
output_link_angles_degrees = []

# Create a plot to visualize the four-bar linkage and path
plt.figure(figsize=(12, 6))

# Set the figure size and axis limits to maintain a constant scale and size
plt.xlim(-5, 15)  # Adjust these values as needed
plt.ylim(-5, 5)   # Adjust these values as needed

# Iterate through the input link angles
for theta1_degrees in input_link_angles_degrees:
    # Clear the current figure
    #plt.clf()

    # Convert input angle from degrees to radians
    theta1 = math.radians(theta1_degrees)

    # Calculate the positions of the joint points
    joint_A_x = 0.0
    joint_A_y = 0.0
    joint_B_x = a * math.cos(theta1)
    joint_B_y = a * math.sin(theta1)

    # Use the vector loop equation to calculate theta2 and beta
    alpha = math.atan2(joint_B_y, joint_B_x)
    cos_beta = (a**2 + b**2 - (joint_B_x**2 + joint_B_y**2)) / (2 * a * b)
    if cos_beta > 1:
        cos_beta = 1  # Ensure cos_beta is within [-1, 1] to avoid math domain errors
    elif cos_beta < -1:
        cos_beta = -1
    beta = math.acos(cos_beta)
    theta2 = alpha + beta

    # Calculate the positions of joint C and D
    joint_C_x = joint_B_x + d * math.cos(theta2)
    joint_C_y = joint_B_y + d * math.sin(theta2)
    joint_D_x = c
    joint_D_y = 0.0

    # Store the output link angle in degrees
    output_link_angles_degrees.append(math.degrees(theta2))

    # Plot the links
    plt.plot([joint_A_x, joint_B_x], [joint_A_y, joint_B_y], 'r-', label='Link AB')
    plt.plot([joint_B_x, joint_C_x], [joint_B_y, joint_C_y], 'g-', label='Link BC')
    plt.plot([joint_C_x, joint_D_x], [joint_C_y, joint_D_y], 'b-', label='Link CD')
    plt.plot([joint_D_x, joint_A_x], [joint_D_y, joint_A_y], 'm-', label='Link DA')

    # Plot the joints
    plt.plot(joint_A_x, joint_A_y, 'ro', label='Joint A')
    plt.plot(joint_B_x, joint_B_y, 'go', label='Joint B')
    plt.plot(joint_C_x, joint_C_y, 'bo', label='Joint C')
    plt.plot(joint_D_x, joint_D_y, 'mo', label='Joint D')

    # Display the plot and introduce a 500 ms delay
    plt.title('Four-Bar Linkage with Varying Input Link Angles')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.legend()
    plt.grid(True)
    #plt.axis('equal')  # Ensure equal scaling on both axes
    plt.draw()
    plt.pause(0.5)  # Delay for 500 ms

# Print the output link angles in degrees
print("Input Link Angles (degrees):", input_link_angles_degrees)
print("Output Link Angles (degrees):", output_link_angles_degrees)

# Keep the plot window open until manually closed
plt.show()





#%%

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling

algorithm = NSGA2(
    pop_size=40,
    n_offsprings=10,
    sampling=FloatRandomSampling(),
    crossover=SBX(prob=0.9, eta=15),
    mutation=PM(eta=20),
    eliminate_duplicates=True
)

from pymoo.termination import get_termination

termination = get_termination("n_gen", 40)

from pymoo.optimize import minimize

res = minimize(problem,
               algorithm,
               termination,
               seed=1,
               save_history=True,
               verbose=True)

X = res.X
F = res.F

#%%
import matplotlib.pyplot as plt
xl, xu = problem.bounds()
plt.figure(figsize=(7, 5))
plt.scatter(X[:, 0], X[:, 1], s=30, facecolors='none', edgecolors='r')
plt.xlim(xl[0], xu[0])
plt.ylim(xl[1], xu[1])
plt.title("Design Space")
plt.show()
#%%
plt.figure(figsize=(7, 5))
plt.scatter(F[:, 0], F[:, 1], s=30, facecolors='none', edgecolors='blue')
plt.title("Objective Space")
plt.show()
