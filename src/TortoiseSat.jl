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
using AttitudeController
using IGRF
# using ApproxFun


#declare functions
include("kep_ECI.jl")
include("NED_cart.jl")
include("DerivFunction.jl")
include("OrbitPlotter.jl")
include("simulator.jl")
include("gain_simulator.jl")
include("quaternion_toolbox.jl")
include("magnetic_toolbox.jl")

#earth parameters
GM = 3.986004418E14*(1/1000)^3 #earth graviational constant (km^3 s^-2)""
R_E = 6371.0 #earth radius (km)

#satellite parameters (1U)
mass=.75 #kg
J=[00.00125 0 0;
    0 0.00125 0;
    0 0 0.00125] #kgm2
J_inv=inv(J)
BC = mass/2.2/(J[1,1]*J[2,2])

#enter orbital parameters
alt=400 #altitude in km
A=zeros(1,6)
#enter orbital parameters
A[1] = 0 #circular orbits
A[2]= alt+R_E #semi-major axis
# A[3] = rand(1)[1]*90 #inclination in degrees
A[3] = 57.6 #inclination in degrees
A[4] = 0 #assign ascending nodes
A[5] = 0 #argument of periapsis
A[6] = 90 #initial true anomoly
T=2*pi*sqrt(A[2].^3/GM)
ω_0 = sqrt(GM/A[2]^3) #rad/s

#Modified Julia Date
MJD_0 = 58155.0 #modified julian day

#magnetic field generation for half orbit
t0 = 0.0 #seconds
tf = 60*90 #seconds
cutoff = 100 #condition number cutoff
N = 5000

#extrapolate over the magnetic field
alt = 400 # km
mag_field = igrf_data(alt,2019)

B_ECI_initial,a,b = magnetic_simulation(A,t0,tf,N,mag_field,GM,MJD_0)

#plot magnetic field
plot(B_ECI_initial)

# Find magnetic control garmian
B_gram = magnetic_gramian(B_ECI_initial,(tf-t0)/N)

#determine at what point the condition cutoff is satisfied
tf_index = condition_based_time(B_gram,cutoff)

#reasign time
t_final = tf_index*(tf-t0)/N

dt = .2 #we want around 5 Hz
N = convert(Int64,floor((t_final-t0)/dt))
t = t0:dt:t_final

# find magnetic field along new time
B_ECI,pos,vel = magnetic_simulation(A,t0,t_final,N, mag_field,GM,MJD_0)

#plot magnetic field
plot(B_ECI[:,1:3])

# local level frame
r = zeros(3,N)
v = zeros(3,N)
h = zeros(3,N)
for i = 1:N
    r[:,i] = pos[:,i]/norm(pos[:,i])
    v[:,i] = vel[:,i]/norm(vel[:,i])
    h[:,i] = cross(r[:,i],v[:,i])
end

q_TF = zeros(4,N)
ROT = zeros(3,3,N)
ω_ll = zeros(3,N)
for i = 1:N
    ROT[:,:,i] = [h[:,i] v[:,i] r[:,i]]'
    q_TF[:,i] = rot_quaternion(ROT[:,:,i])
    ω_ll[:,i] = qrot(q_TF[:,i],h[:,i]*ω_0)
end

#trajectory optimization section
#initial
omega_int = zeros(3)
# q_0 = q_TF[:,1]
θ = 90 # angle in degrees
r = [0 0 1]'/norm([0 0 1]) #normalized rotation axis
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

#consolidate the model
model=Model(DerivFunction,n,m,quaternion_error,quaternion_expansion)

#use eigen-axis slew as first guess for state
ω_guess, q_guess = eigen_axis_slew(x0[1:7],xf[1:7],t)
plot(ω_guess')
plot(q_guess')
X = [ω_guess;q_guess;reshape(t,1,length(t))]

#Bryson's Rule LQR weights
Q = zeros(n,n)
Qf = zeros(n,n)
ω_max = maximum(abs.(X[1:3,:]))
τ_max = maximum(J * diff(X[1:3,:],dims = 2)/dt)
m_max = τ_max / 1.e-5 * 1.e2
α = 1.e1
Q[1:3,1:3] = Array((α/ω_max^2)*Diagonal(I,3))
Qf[1:3,1:3] = Array((α/ω_max^2)*Diagonal(I,3))*10
β = 1.e3
Q[4:7,4:7] = Array(α*β*Diagonal(I,4))
Qf[4:7,4:7] = Array(α*β*Diagonal(I,4))*10
R = Array((1/m_max^2)*Diagonal(I,m))

#bounds
u_bnd=[10;10;10]
x_bnd=10

obj = LQRObjective(Q, R, Qf, t_final, x0, xf)
obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
solver = TrajectoryOptimization.Solver(model,obj_con,dt = dt)
solver.opts.verbose=true
solver.opts.iterations_outerloop=5
solver.opts.iterations_innerloop=30
solver.opts.dJ_counter_limit = 1
solver.opts.sat_att = true
# X = rand(Float64,n,length(t))/1000 #initialize random input and state
U = rand(Float64,m,length(t))/1000 #initialize random input and state

@time results, stats = TrajectoryOptimization.solve(solver,U)
X = TrajectoryOptimization.to_array(results.X) #assigns state and input
U = TrajectoryOptimization.to_array(results.U)


#plot input
plt = plot(t[1:end-1],U[1:3,:]'./100,label = ["U1" "U2" "U3"])

plot!(title = "Swarm Slew", xlabel = "Time (s)", ylabel = "Magnetic Moment(Am2)")

#plot omega
plot(t,X[1:3,:]')

#plot quaternion
plot(t,X[4:7,:]')


#plot time
plot(X[8,:])

#now test!

#assigns initial conditions
x0_lqr = .38*pi/180
x0_lqr[1:3] = x0[1:3]

#creates noise over initial conditions
q_noise = randn(3,1)*(1*pi/180)^2 #assign noise to initial quaternion
θ_noise = norm(q_noise)
r_noise = q_noise/θ_noise
x0_lqr[4:7] = qmult(x0[4:7],[cos(θ_noise/2);r_noise*sin(θ_noise/2)])

# #assigns new time vector
# interp_factor = 5
# N_lqr = N*interp_factor
# dt_lqr = dt/interp_factor
# t_lqr = [t0:dt_lqr:tf...]
#
#define integration ds
integration=:rk4
#
# #interpolate to find new state and input
# X_lqr, U_lqr = interpolate_trajectory( integration, X, U, t_lqr)


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

#plotting
plt = plot(t,X_sim[1:3,:]',label = ["omega1" "omega2" "omega3"])
plot!(title = "Swarm Slew", xlabel = "Time (s)", ylabel = "Rotation Rate (rad/s)")

plt = plot(t,X_sim[4:7,:]',label = ["q1" "q2" "q3" "q4"])
plot!(title = "Swarm Slew", xlabel = "Time (s)", ylabel = "Quaternion")

plt = plot(t,U_sim[1:3,:]'./100,label = ["u1" "u2" "u3"])
plot!(title = "Swarm Slew", xlabel = "Time (s)", ylabel = "Magnetic Moment (Am2)")

plot(dX[1:3,:]')
plot(dX[4:6,:]')


plot(U_sim',color="red")
plot!(U_lqr',color="blue")

# find error for comparison
TVLQR_error = zeros(length(t))
for i = 1:length(t)
    euler_angles = quaternion_euler(X_sim[4:7,i])
    TVLQR_error[i] = norm(euler_angles-quaternion_euler(q_final))
end

plot(t/60/60,TVLQR_error)

# plot(B_ECI,color="red")
# plot!(B_ECI_initial,color="blue")

using HDF5
using CSV
h5write(string("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/visualize.h5"), "X_sim", X_sim)
h5write(string("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/visualize.h5"), "B_ECI", B_ECI)

h5read("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/visualize.h5", "X_sim")

function L_decomp(q)
    #this function takes in a desired quaternion computes the L decomposition
    L = [q[1] -q[2] -q[3] -q[4];
         q[2] q[1] -q[4] q[3];
         q[3] q[4] q[1] -q[2];
         q[4] -q[3] q[2] q[1]]

         return L

end
