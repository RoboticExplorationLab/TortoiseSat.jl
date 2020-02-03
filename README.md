# TortoiseSat
Under-actuated control of satellite attitude using magnetorquers and trajectory optimization

## Summary
This code base was designed for research into utilizing Differential Dynamic Programming, and specifically the Trajectory Optimization package being worked on by the Robotic Exploxartion Lab at Stanford University, to magnetorquer-only control of satellites. Satellite dynamics with magnetorquers are naturally underactuated, and they can easily be formulated in the format given by the conference paper on this topic: (rexlab.stanford.edu/papers/magnetorquer_only.pdf). In this format, the dynamics look very similar to a classically underactuated problem in robotics, such as limbs and other actuators. Therefore, the solvers that are used for motion planning for these purposes were reapplied to spacecraft dynamics. 

## Results
The application of the Augmented Lagrangian/ iterative LQR control scheme developed through the TrajectoryOpitmization.jl package (found here: github.com/RoboticExplorationlab/TrajectoryOptimization.jl) was immediately useful in planning spacecraft slews. Monte-Carlo testing of the control loop found that it was applicable in all orbits with varying utility, as orbits that present more magnetic diversity are naturally more useful. 

## Getting Started 
Currently, the codebase is not terribly user friendly, but if you have your own satellite parameters and orbital parameters that you would like to test, please navigate to the TortoiseSat.jl file in the src folder. This file has the most documentation and it should be relatively straightforward to replace the values that you need to work with in order to appropriately simulate your own system for basic slews using this method. 