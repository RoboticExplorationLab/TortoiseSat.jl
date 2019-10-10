#this code is for simulating the dynamics and
#control of the "Magnetic Torquer Attitude Control
# via Asymptotic Periodic Linear Quadratic Regulation"
# paper

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
using ControlSystems


#Include funtions
include("/Users/agatherer/.julia/dev/TortoiseSat/src/kep_ECI.jl")
include("/Users/agatherer/.julia/dev/TortoiseSat/src/NED_cart.jl")
include("/Users/agatherer/.julia/dev/TortoiseSat/src/DerivFunction.jl")
include("/Users/agatherer/.julia/dev/TortoiseSat/src/OrbitPlotter.jl")
include("/Users/agatherer/.julia/dev/TortoiseSat/src/simulator.jl")
include("/Users/agatherer/.julia/dev/TortoiseSat/src/gain_simulator.jl")
include("/Users/agatherer/.julia/dev/TortoiseSat/src/magnetic_toolbox.jl")
include("/Users/agatherer/.julia/dev/TortoiseSat/src/quaternion_toolbox.jl")
include("psiaki_dynamics.jl")
include("/Users/agatherer/.julia/dev/TortoiseSat/src/attitude_dynamics.jl")

#earth parameters
GM = 3.986004418E14*(1/1000)^3 #earth graviational constant (km^3 s^-2)""
R_E = 6371.0 #earth radius (km)

#satellite parameters (1U)
mass=.75 #kg
J=Diagonal(.00125*ones(3)) #kgm2
J_inv=inv(J)
BC = mass/2.2/(J[1,1]*J[2,2])

# # satellite parameters (big boy)
# J=Diagonal([8.7,10,6.5]) #kgm2
# mass = 100 # kg


#enter orbital parameters
alt=400 #altitude in km
kep=zeros(1,6)
#enter orbital parameters
kep[1] = 0 #circular orbits
kep[2]= alt+R_E #semi-major axis
kep[3] = 88 #inclination in degrees
kep[4] = 0 #assign ascending nodes
kep[5] = 0 #argument of periapsis
kep[6] = 0 #initial true anomoly
T=2*pi*sqrt(kep[2].^3/GM)

A = kep

#calculate the orbital mean motion
ω_0 = sqrt(GM/kep[2]^3) #rad/s


#Modified Julia Date
MJD_0 = 58155.0 #modified julian day

#magnetic field generation for full orbit
t0 = 0.0 #seconds
tf = 20T #seconds
N = 100000
t = t0:(tf-t0)/N:tf

#extrapolate over the magnetic field
mag_field = igrf_data(alt,2019)

#find the ECI magnetic field and the position and velocity
B_ECI,pos,vel = magnetic_simulation(kep,t0,tf,N,mag_field,GM,MJD_0)
B_ECI = B_ECI'
plot(B_ECI[:,1:N]')

#initialize the 3 unit vectors for the trajectory frame
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
    ω_ll[:,i] = ROT[:,:,i]*h[:,i]*ω_0
end

#plot magnetic field
B_ll = zeros(3,N)
for i in 1:N
    B_ll[:,i] = ROT[:,:,i] * B_ECI[:,i]
end
plot(B_ll')

# find simplified magnetic field
B_ll_simplified = zeros(3,N)
i_m = kep[3]/180*pi # inclination in radians
a = kep[2]*1000 # semi-major in m
μ_f = 7.9E15
for i in 1:N
    B_ll_simplified[:,i] = μ_f/a^3*
        [cos(ω_0*t[i])*sin(i_m);-cos(i_m);2*sin(ω_0*t[i])*sin(i_m)]
end
plot(B_ll_simplified')

#initial
omega_int=ω_ll[:,1]
# q_0 = euler_quaternion([.34;0;0])
q_0 = [sqrt(2)/2;sqrt(2)/2;0;0]

#final
q_final=[1.;0;0;0]
# q_final = q_0
omega_final=[0.;0;0]

# plot q_final over time
# q_final_permuted = zeros(size(q_TF))
# for i = 1:N
#     q_final_permuted[:,i] = qmult(q_final,q_TF[:,i])
# end
# plot(q_final_permuted')

#generate time vector
dt = (tf-t0)/N
t = t0:dt:tf
t_psiaki = t

#initialize state
x = zeros(7,length(t))
x[:,1] = [omega_int;q_0]

#extrapolate state once to get higher order magnetic field info
B_meas=qrot(q_inv(x[4:7,1]),B_ll_simplified[:,1])
ẋ = attitude_dynamics(x[:,1],zeros(3),B_meas,J)
x[:,2] = x[:,1] + dt * ẋ

# find the magnetic control effectiveness
α = 3000
U_max = .1*ones(3)
R = Diagonal(8E6*ones(3))
R_inv = inv(R)
Q = Diagonal([1.5E-4;1.5E-4;1.5E-4;1E4;1E4;1E4])
B̃ = real.(magnetic_control_effectiveness(B_ll_simplified,J,R,t))

# check to make sure controllable!
# rank(ctrb(A,B̃))

# solve the Riccatti equaiton to finde Pss
using ControlSystems
Pss = Pss_psiaki(ω_0,J,Q,B̃)

# apply control law
U = zeros(3,N)
B_matrix = zeros(6,3,N)
x_linear = zeros(6,N)
B_meas = zeros(3,N)
for i in 1:N

    # extrapole local magnetic field
    B_meas[:,i]=qrot(x[4:7,i],B_ll_simplified[:,i])

    # calculate the B matrix in the linearized dyanmics
    B_matrix[4:6,1:3,i] = -hat(B_meas[:,i])
    B_matrix[4,:,i] *= 1/J[1,1]
    B_matrix[5,:,i] *= 1/J[2,2]
    B_matrix[6,:,i] *= 1/J[3,3]

    # extrapolate the linearized state
    x_linear[1:3,i] = quaternion_euler(x[4:7,i])
    x_linear[4:6,i] = x[1:3,i] - qrot(x[4:7,i],ω_ll[:,i])

    # calculate the control
    U[:,i] = nominal_input_psiaki(R_inv,B_matrix[:,:,i],Pss,x_linear[:,i],
        α,U_max)

    # extrapolate dynamics
    ẋ = attitude_dynamics_linear(x[:,i], U[:,i],x_linear[:,i], qrot(x[4:7,i],B_ll[:,i]), J)

    x[:,i+1] = rk4_psiaki(attitude_dynamics,x[:,i],dt,U[:,i],B_meas[:,i],J)
    x[4:7,i+1]=normalize(x[4:7,i+1])

end

# find the difference
plot(x_linear[1:3,:]')
plot(x_linear[4:6,:]')
plot(U')

# TRAJOPT SECTION ------------------------
# ----------------------------------------------------

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

# # satellite parameters (Swarm)
# mass=.2 #kg
# J=[5e-4 0 0;
#     0 4e-4 0;
#     0 0 9e-4] #kgm2
# J_inv=inv(J)
# BC = mass/2.2/(J[1,1]*J[2,2])


#enter orbital parameters
alt=400 #altitude in km
A=zeros(1,6)
#enter orbital parameters
A[1] = 0 #circular orbits
A[2]= alt+R_E #semi-major axis
# A[3] = rand(1)[1]*90 #inclination in degrees
A[3] = 88 #inclination in degrees
A[4] = 0 #assign ascending nodes
A[5] = 0 #argument of periapsis
A[6] = 0 #initial true anomoly
T=2*pi*sqrt(A[2].^3/GM)
ω_0 = sqrt(GM/A[2]^3) #rad/s

#Modified Julia Date
MJD_0 = 58155.0 #modified julian day

#magnetic field generation for half orbit
t0 = 0.0 #seconds
tf = 60*90 #seconds
cutoff = 5 #condition number cutoff
N = 5000

#extrapolate over the magnetic field
alt = 400 # km
mag_field = igrf_data(alt,2019)

B_ECI_initial,pos,vel = magnetic_simulation(A,t0,tf,N,mag_field,GM,MJD_0)

# local level frame
r = zeros(3,2N)
v = zeros(3,2N)
h = zeros(3,2N)
for i = 1:2N
    r[:,i] = pos[:,i]/norm(pos[:,i])
    v[:,i] = vel[:,i]/norm(vel[:,i])
    h[:,i] = cross(r[:,i],v[:,i])
end

q_TF = zeros(4,2N)
ROT = zeros(3,3,2N)
ω_ll = zeros(3,2N)
for i = 1:2N
    ROT[:,:,i] = [h[:,i] v[:,i] r[:,i]]'
    q_TF[:,i] = rot_quaternion(ROT[:,:,i])
    ω_ll[:,i] = qrot(q_TF[:,i],h[:,i]*ω_0)
end

# create local level frame mag field
B_ll_initial = zeros(size(B_ECI_initial))
for i = 1:2N
    B_ll_initial[i,:] = (B_ECI_initial[i,:]'*ROT[:,:,i])'
end

#plot magnetic field
plot(B_ECI_initial)
plot(B_ll_initial)

# Find magnetic control garmian
B_gram = magnetic_gramian(B_ll_initial,(tf-t0)/N)

#determine at what point the condition cutoff is satisfied
tf_index = condition_based_time(B_gram,cutoff)

#reasign time
t_final = tf_index*(tf-t0)/N

dt = .2 #we want around 5 Hz
N = convert(Int64,floor((t_final-t0)/dt))
t = t0:dt:t_final

# find magnetic field along new time
B_ECI,pos,vel = magnetic_simulation(A,t0,t_final,N, mag_field,GM,MJD_0)
B_ll = zeros(size(B_ECI))
for i = 1:size(B_ECI,1)
    B_ll[i,:] = (B_ECI[i,:]'*ROT[:,:,i])'
end

#plot magnetic field
plot(B_ECI[:,1:3])
plot(B_ll[:,1:3])
# now cheat
B_ECI = B_ll

#trajectory optimization section
#initial
# omega_int = [ω_0 * h[:,1]]
omega_int = zeros(3)
# q_0 = euler_quaternion([.34;0;0])
q_0 = [sqrt(2)/2;sqrt(2)/2;0;0]
x0=[omega_int;q_0;t0]

#final
q_final= [1.;0;0;0]
# omega_final= ω_0 * h[:,1]
omega_final = zeros(3)
xf=[omega_final;q_final;1] #km or km/s

#create model
n=8
m=3

#consolidate the model
model=Model(DerivFunction,n,m,quaternion_error,quaternion_expansion)

# #use eigen-axis slew as first guess for state
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
# R = Array(1.e-10*Diagonal(I,m))

# #LQR weights
# Q=zeros(n,n)
# Qf=zeros(n,n)
# α = 1.e1
# Q[1:3,1:3] = Array((α)*Diagonal(I,3))
# Qf[1:3,1:3] = Array((α)*Diagonal(I,3))*100
# β = 1.e1
# Q[4:7,4:7] = Array(β*Diagonal(I,4))
# Qf[4:7,4:7] = Array(β*Diagonal(I,4))*100
# R = Array((7.5e-4*t_final)*Diagonal(I,m))

#bounds
u_bnd=[10;10;10]
x_bnd=10

obj = LQRObjective(Q, R, Qf, t_final, x0, xf)
obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
solver = TrajectoryOptimization.Solver(model,obj_con,dt = dt)
solver.opts.verbose=true
solver.opts.iterations_outerloop=5
solver.opts.iterations_innerloop=50
solver.opts.dJ_counter_limit = 1
solver.opts.sat_att = true
X = rand(Float64,n,length(t))/1000 #initialize random input and state
U = rand(Float64,m,length(t))/1000 #initialize random input and state

@time results, stats = TrajectoryOptimization.solve(solver,U)
X = TrajectoryOptimization.to_array(results.X) #assigns state and input
U = TrajectoryOptimization.to_array(results.U)

t_TO = t

plt = plot(t_TO,X[1:3,:]',label = ["omega1" "omega2" "omega3"])
plot!(title = "Swarm Slew", xlabel = "Time (s)", ylabel = "Rotation Rate (rad/s)")

plt = plot(t_TO,X[4:7,:]',label = ["q1" "q2" "q3" "q4"])
plot!(title = "Swarm Slew", xlabel = "Time (s)", ylabel = "Quaternion")

plt = plot(t_TO[1:end-1],U[1:3,:]'./100,label = ["u1" "u2" "u3"])
plot!(title = "Swarm Slew", xlabel = "Time (s)", ylabel = "Magnetic Moment (Am2)")


# calculate error for comparisons
psiaki_error = zeros(length(t_psiaki))
for i = 1:length(t_psiaki)-1
    psiaki_error[i] = norm(x_linear[1:3,i])
end

TVLQR_error = zeros(length(t_psiaki))
for i = 1:length(t_psiaki)-1
    if i < length(t_TO)
        TVLQR_error[i] = 2*acos(min(1,qmult(q_inv(X[4:7,i]),q_final)[1]))
    else
        TVLQR_error[i] = 2*acos(min(1,qmult(q_inv(X[4:7,end]),q_final)[1]))
    end
end

plot_length = 8000

plot(t_psiaki[1:plot_length]/60/60,psiaki_error[1:plot_length]*180/pi,label = "APLQR")
plot!(t_psiaki[1:plot_length]/60/60,TVLQR_error[1:plot_length]*180/pi,label = "Trajectory Optimization")
xlabel!("Time (Hours)")
ylabel!("Attitude Error (deg)")
savefig("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/Psiaki_comparison.png")
