#this code is for trajectory optimization of an underactuated satellite

#install packages
using Pkg;
using LinearAlgebra;
using Plots;
using TrajectoryOptimization;
using DifferentialEquations;
using SatelliteToolbox
using ForwardDiff
using SparseArrays
using Interpolations
using AttitudeController
using IGRF


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

# #satellite parameters (3U)
# mass=3.5 #kg
# J=[0.005256 0 0;
#     0 0.04939 0;
#     0 0 0.04939] #kgm2
# J_inv=inv(J)

#enter orbital parameters
alt=400 #altitude in km
A=zeros(1,6)
#enter orbital parameters
A[1] = 0 #circular orbits
A[2]= alt+R_E #semi-major axis
A[3] = rand(1)[1]*180 #inclination in degrees
A[4] = rand(1)[1]*360 #assign ascending nodes
A[5] = 0 #argument of periapsis
A[6] = rand(1)[1]*360 #initial true anomoly
T=2*pi*sqrt(A[2].^3/GM)

#Modified Julia Date
MJD_0 = 58155.0 #modified julian day

#magnetic field generation for half orbit
t0 = 0.0 #seconds
tf = 60*40 #seconds
cutoff = 20 #condition number cutoff
N = 1000
B_N_sim = magnetic_simulation(A,t0,tf,N)

#plot magnetic field
plot(B_N_sim[:,1:3])

# Find magnetic control garmian
B_gram = magnetic_gramian(B_N_sim)

#determine at what point gramian eigenvalues meet criteria
tf_index = condition_based_time(B_gram,cutoff)
tf = tf_index*(tf-t0)/N

#reasign time
# find magnetic field along new time
N = 1000
B_N_sim = magnetic_simulation(A,t0,tf,N)

#plot magnetic field
plot(B_N_sim[:,1:3])

#trajectory optimization section
#initial
omega_int=[0;0;0]
q_0=[0.;0;1;0]
x0=[omega_int;q_0;t0]

#final
q_final=[sqrt(2)/2;sqrt(2)/2;0;0]
omega_final=[0.;0;0]
xf=[omega_final;q_final;1] #km or km/s

#create model
n=8;
m=3;

#consolidate the model
model=Model(DerivFunction,n,m,quaternion_error,quaternion_expansion)

#LQR weights
Q=zeros(n,n)
Qf=zeros(n,n)
α = 1.e1
Q[1:3,1:3] = Array((α)*Diagonal(I,3))
Qf[1:3,1:3] = Array((α)*Diagonal(I,3))*100
β = 1.e1
Q[4:7,4:7] = Array(β*Diagonal(I,4))
Qf[4:7,4:7] = Array(β*Diagonal(I,4))*100

R = Array((1.e-3/norm(B_N_sim))*Diagonal(I,m))

#bounds
u_bnd=[1.9;1.9;1.9]
x_bnd=10

#random start to trajopt
# cost = GenericCost(quat_cost,quat_cost_terminal,quat_cost_grad,quat_cost_hess,n,m)
# obj = UnconstrainedObjective(cost,tf,x0,xf)
obj = LQRObjective(Q, R, Qf, tf, x0, xf)
obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
solver = TrajectoryOptimization.Solver(model,obj_con,N=N)
solver.opts.verbose=true
solver.opts.iterations_outerloop=3
solver.opts.iterations_innerloop=50
solver.opts.dJ_counter_limit = 1
solver.opts.sat_att = true
U = rand(Float64,m,solver.N)/1000 #initialize random input and state
# X = rand(Float64,m,solver.N)/1000

#use eigen-axis slew as first guess for state
ω_guess, q_guess = eigen_axis_slew(x0[1:7],xf[1:7],t[1:N])
X = [ω_guess[:,1:N];q_guess[:,1:N];reshape(t[1:N],1,N)]

results, stats = TrajectoryOptimization.solve(solver,U)
X = TrajectoryOptimization.to_array(results.X) #assigns state and input
U = TrajectoryOptimization.to_array(results.U)


#plot input
plt = plot(t[1:N-1],U[1:3,:]',label = ["U1" "U2" "U3"])
plot!(title = "1U CubeSat Slew", xlabel = "Time (s)", ylabel = "Magnetic Moment(Am2)")

#plot magnetic gramian for reference
plot(B_gram_cond',yaxis=:log)

#plot omega
plot(t[1:N],X[1:3,:]')

#plot quaternion
plot(t[1:N],X[4:7,:]')

#plot time
plot(X[8,:])

#now test!

#assigns initial conditions
x0_lqr = zeros(n)
x0_lqr[1:3] = x0[1:3]

#creates noise over initial conditions
q_noise = randn(3,1)*(.0001*pi/180)^2 #assign noise to initial quaternion
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
α = 1.e-1
β = 1.e3
Q_lqr = zeros(6,6)
Q_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))
Q_lqr[4:6,4:6] = Array((α*β)*Diagonal(I,3))
Qf_lqr = zeros(6,6)
Qf_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))*1.e4
Qf_lqr[4:6,4:6] = Array((α*β)*Diagonal(I,3))*1.e8
R_lqr = Array((2.e-2)*Diagonal(I,3))

#executes TVLQR
X_sim, U_sim, dX, K = attitude_simulation(f!,f_gains!,integration,X_lqr,U_lqr,dt_lqr,x0_lqr,t0,tf,Q_lqr,R_lqr,Qf_lqr)

#plotting
plt = plot(t[1:N],X_sim[1:3,:]',label = ["omega1" "omega2" "omega3"])
plot!(title = "3U CubeSat Slew", xlabel = "Time (s)", ylabel = "Rotation Rate (rad/s)")

plt = plot(t[1:N],X_sim[4:7,:]',label = ["q1" "q2" "q3" "q4"])
plot!(title = "3U CubeSat Slew", xlabel = "Time (s)", ylabel = "Quaternion")

plt = plot(t[1:N],U_sim[1:3,:]'./100,label = ["u1" "u2" "u3"])
plot!(title = "3U CubeSat Slew", xlabel = "Time (s)", ylabel = "Magnetic Moment (Am2)")

plot(dX[1:3,:]')
plot(dX[4:6,:]')


plot(U_sim',color="red")
plot!(U_lqr',color="blue")

# plot(B_N,color="red")
# plot!(B_N_sim,color="blue")


function L_decomp(q)
    #this function takes in a desired quaternion computes the L decomposition
    L = [q[1] -q[2] -q[3] -q[4];
         q[2] q[1] -q[4] q[3];
         q[3] q[4] q[1] -q[2];
         q[4] -q[3] q[2] q[1]]

         return L

end
