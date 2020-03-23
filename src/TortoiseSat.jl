#this code is for trajectory optimization of an underactuated satellite

#install packages
using Pkg
using LinearAlgebra
using Plots
using TrajectoryOptimization
using DifferentialEquations
using SatelliteToolbox
using ForwardDiff
using SparseArrays
using Interpolations
# using AttitudeController
# using IGRF
# using ApproxFun


#declare functions 
include("src/kep_ECI.jl")
include("src/NED_cart.jl")
include("src/DerivFunction.jl")
include("src/OrbitPlotter.jl")
include("src/simulator.jl")
include("src/gain_simulator.jl")
include("src/quaternion_toolbox.jl")
include("src/magnetic_toolbox.jl")
include("src/input_parameters.jl")
include("src/eigen_axis_slew.jl")
include("src/attitude_controller.jl")

R_E = 6178 #km

# Keplerian elements as a starting point
alt = 400 #km
Kep=zeros(1,6)
#enter orbital parameters
Kep[1] = 0 #eccentricity
Kep[2]= alt+R_E #semi-major axis (km)
Kep[3] = 96 #inclination (degrees)
Kep[4] = 0 #assign ascending nodes (degrees)
Kep[5] = 0 #argument of periapsis (degrees)
Kep[6] = 90 #initial true anomoly (degrees)

MJD = 58155.0 #modified julian day

# generate satellite parameters
# all satellite parameters are stored inside a struct to make it easier to reference them
# for different satellites 
p = input_parameters("1P",Kep,MJD)

# First, we need to generate the magnetic field for half an orbit to scope out how long 
# we want the simulation to go for. For harder orbits, we'll need longer simulations to 
# achieve the same control authority. For really nice orbits like Polar orbits, the 
# simulation will naturally be smaller since the satellite has more magnetic diversity to 
# work with. This is all captured in the magnetic gramian calculation

## magnetic field generation for half orbit
t0 = 0.0 #seconds
tf = 60*90 #seconds
cutoff = 50 #condition number cutoff
N = 5000

#extrapolate over the magnetic field
mag_field = igrf_data(p.alt,2019)
B_ECI_initial,a,b = magnetic_simulation(p,t0,tf,N,mag_field)

#plot magnetic field
plot(B_ECI_initial)

# Now, since we have the magnetic field for half an orbit, we calculate the gramian

# Find magnetic control gramian
B_gram = magnetic_gramian(B_ECI_initial,(tf-t0)/N)

#determine at what point the condition cutoff is satisfied (heuristic)
tf_index = condition_based_time(B_gram,cutoff)

# Finally, now we know how long we want our simulation to last, we can recalculate the 
# magnetic field along that orbit in much higher detail for our solver 

#reasign time
t_final = tf_index*(tf-t0)/N

dt = .2 #we want around 5 Hz
N = convert(Int64,floor((t_final-t0)/dt))
t = t0:dt:t_final

# find magnetic field along new simulation 
B_ECI,pos,vel = magnetic_simulation(p,t0,t_final,N,mag_field)

#plot magnetic field
plot(B_ECI[:,1:3])

# # local level frame ( if needed)
# r = zeros(3,N)
# v = zeros(3,N)
# h = zeros(3,N)
# for i = 1:N
#     r[:,i] = pos[:,i]/norm(pos[:,i])
#     v[:,i] = vel[:,i]/norm(vel[:,i])
#     h[:,i] = cross(r[:,i],v[:,i])
# end

# q_TF = zeros(4,N)
# ROT = zeros(3,3,N)
# ω_ll = zeros(3,N)
# for i = 1:N
#     ROT[:,:,i] = [h[:,i] v[:,i] r[:,i]]'
#     q_TF[:,i] = rot_quaternion(ROT[:,:,i])
#     ω_ll[:,i] = qrot(q_TF[:,i],h[:,i]*p.ω_0)
# end

# Now comes the trajectory optimization. This code uses the julia TrajectoryOptimization 
# package, which abstracts a lot of the code away from the user. This makes calling it easy 
# but sometimes opaque. For our purposes, it's nice to not have to worry about it

#trajectory optimization section
#initial
omega_int = zeros(3)
# q_0 = q_TF[:,1]
θ = 90 # angle in degrees
r = [1 0 1]'/norm([1 0 1]) #normalized rotation axis
q_0 = [cosd(θ/2);r*sind(θ/2)][:]
x0=[omega_int;q_0;t0]

#final
# q_final= q_TF[:,end]
q_final = [cosd((θ-90)/2);r*sind((θ-90)/2)][:]
omega_final= zeros(3)
xf=[omega_final;q_final;1] #km or km/s

#create model
n=8
m=3

# #assign quat type
# info = Dict{Symbol,Any}()
# quat = true
# if quat
#     info[:quat] = [4:7]
# end
# #consolidate the model
# model=Model(DerivFuanction,n,m,info)

model=Model(DerivFunction, n, m)
model_d = TrajectoryOptimization.rk3(model)

# #use eigen-axis slew as first guess for state
ω_guess, q_guess = eigen_axis_slew(x0[1:7],xf[1:7],t)
plot(ω_guess)
plot(q_guess)
X = [ω_guess';q_guess';reshape(t,1,length(t))]
U0 = zeros(3,length(t))
#U0 = rand(Float64,m,length(t))/1000 #initialize random input

#Bryson's Rule LQR weights
Q = zeros(n,n)
Qf = zeros(n,n)
ω_max = maximum(abs.(X[1:3,:]))
τ_max = maximum(p.J * diff(X[1:3,:],dims = 2)/dt)
m_max = τ_max / 1.e-5 * 1.e2
α = 1.e1
Q[1:3,1:3] = Array((α/ω_max^2)*Diagonal(I,3))
Qf[1:3,1:3] = Array((α/ω_max^2)*Diagonal(I,3))*10
β = 1.e3
Q[4:7,4:7] = Array(α*β*Diagonal(I,4))
Qf[4:7,4:7] = Array(α*β*Diagonal(I,4))*10
R = Array((1/m_max^2)*Diagonal(I,m))
obj = LQRObjective(Q, R, Qf, xf, N)

#bounds
# for the control constraints, I have resized the boundaries to make units match better
# and not run into as many numerical issues, so you enter the bounds in .01's of Am2

# So, to summarize, take your input constraint (let's say 1 Am2) and multiply it by 100 to 
# get what you should put down here

bnd = BoundConstraint(n,m,u_max=1, u_min=-1) # control limits in .01's of Am2
x_bnd=10

#goal
goal = goal_constraint(xf) # terminal constraint

constraints = Constraints(N) # define constraints at each time step
for k = 1:N-1
    constraints[k] += bnd
end
constraints[N] += goal

sat = Problem(model_d, obj, constraints = constraints, x0=x0, xf=xf, N=N, dt=dt);
initial_controls!(sat,U0); # initialize problem with controls

# Set up the solver
opts_al = AugmentedLagrangianSolverOptions{Float64}()
opts_al.opts_uncon.iterations = 50  # set iLQR iterations
opts_al.iterations = 20      # set outer loop iterations
solver = TrajectoryOptimization.AugmentedLagrangianSolver(sat,opts_al)

TrajectoryOptimization.solve!(sat,solver)

# reassign inputs and states to arrays (they are currently arrays of arrays)
X = hcat([sat.X[i] for i = 1:size(sat.X,1)]...)
U = hcat([sat.U[i] for i = 1:size(sat.U,1)]...)

#plot input
plt = plot(t[1:end-2],U'./100,label = ["U1" "U2" "U3"])

plot!(title = "Slew", xlabel = "Time (s)", ylabel = "Magnetic Moment (Am^2)")

#plot omega
plt = plot(t[1:end-1],(X[1:3,:]*180/pi)')

plot!(title = "Slew", xlabel = "Time (s)", ylabel = "Rotation Rate (deg/s)")

#plot quaternion
plot(t[1:end-1],X[4:7,:]')

plot!(title = "Slew", xlabel = "Time (s)", ylabel = "Unitless")


## Now we have an optimal trajectory, we want to throw it in a closed loop and find out how 
# well it does against an environment with simulated noise. We'll impose a closed control loop 
# in the form of Time Varying LQR (TVLQR) with the gains that we got out of the 
# TrajectoryOptimization solver 

#assigns initial rotation rate conditions
x0_lqr = zeros(size(x0))
x0_lqr[1:3] = x0[1:3]

#creates noise over initial attitude conditions
q_noise = randn(3,1)*(1*pi/180)^2 #assign noise to initial quaternion
θ_noise = norm(q_noise)
r_noise = q_noise/θ_noise
x0_lqr[4:7] = qmult(x0[4:7],[cos(θ_noise/2);r_noise*sin(θ_noise/2)])

#define integration ds
integration=:rk4


#since we are no longer interpolating, just reassign
X_lqr = X
U_lqr = U
dt_lqr = dt

#assigns simulation and gains functions for TVLQR
include("simulator.jl")
f! = simulator
f_gains! = gain_simulator

#creates Q and R matrices for TVLQR
α = 1.e1
Q_lqr=zeros(6,6)
Qf_lqr=zeros(6,6)
Q_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))
Qf_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))*100
β = 1.e1
Q_lqr[4:6,4:6] = Array(β*Diagonal(I,3))
Qf_lqr[4:6,4:6] = Array(β*Diagonal(I,3))*100

R_lqr = Array((7.5e3)*Diagonal(I,m))

#executes TVLQR
X_sim, U_sim, dX, K = attitude_simulation(f!,f_gains!,
    integration,X_lqr,U_lqr,dt_lqr,x0_lqr,t[1],t[end],Q_lqr,
    R_lqr,Qf_lqr)

# These plots represent the actual closed loop responses of the satellite 
plt = plot(t[1:end-1],(X_sim[1:3,:]*180/pi)',label = ["omega1" "omega2" "omega3"])
plot!(title = "Slew", xlabel = "Time (s)", ylabel = "Rotation Rate (deg/s)")

plt = plot(t[1:end-1],X_sim[4:7,:]',label = ["q1" "q2" "q3" "q4"])
plot!(title = "Slew", xlabel = "Time (s)", ylabel = "Quaternion")

plt = plot(t[1:end-1],U_sim[1:3,:]'./100,label = ["u1" "u2" "u3"])
plot!(title = "Slew", xlabel = "Time (s)", ylabel = "Magnetic Moment (Am^2)")

# We can also inspect the deltas (differences) between the control and the open loop 
plot(dX[1:3,:]') # rotation rate
plot(dX[4:6,:]') # quaternion 

# We can also plot the different controls on top of each other to see if they were 
# substantially different
plot(U',color="red")
plot!(U_lqr',color="blue")

# # find error for comparison
# TVLQR_error = zeros(length(t))
# for i = 1:length(t)
#     euler_angles = quaternion_euler(X_sim[4:7,i])
#     TVLQR_error[i] = norm(euler_angles-quaternion_euler(q_final))
# end

# plot(t/60/60,TVLQR_error)

# # plot(B_ECI,color="red")
# # plot!(B_ECI_initial,color="blue")

# using HDF5
# using CSV
# h5write(string("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/visualize.h5"), "X_sim", X_sim)
# h5write(string("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/visualize.h5"), "B_ECI", B_ECI)

# h5read("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/visualize.h5", "X_sim")

# function L_decomp(q)
#     #this function takes in a desired quaternion computes the L decomposition
#     L = [q[1] -q[2] -q[3] -q[4];
#          q[2] q[1] -q[4] q[3];
#          q[3] q[4] q[1] -q[2];
#          q[4] -q[3] q[2] q[1]]

#          return L

# end
