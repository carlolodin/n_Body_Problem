# n_Body_Problem
A MATLAB n-body problem simulation created by Carlo Lodin &amp; Giacomo Fundarò.

Use n_body.m to launch the main simulation.
You can change the number of bodies (n), the gravitational constant (G), the total simulation time (T), the time step (dt), and the sun mass (sun_mass)

n_body_chaotic is used to showcase the chaotic nature of the problem by running two parallel simulation with initial conditions differing by the bodies starting positions by a factor of epsilon (by default equal to 1e-4).

n_body_comparison.m is used to showcase the performances of the four integration schemes implemented.

The other files are required functions used to implement the numerical methods and to update the accelerations at each step.
