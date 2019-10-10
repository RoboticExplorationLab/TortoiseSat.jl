#this code is for simulating the dynamics and
#control of the "Design and testing of magnetic controllers
# for Satellite stabilization" published by Psiaki
# in 2005 using an Israeli satellite

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
J=[00.00125 0 0;
    0 0.00125 0;
    0 0 0.00125] #kgm2
J_inv=inv(J)
BC = mass/2.2/(J[1,1]*J[2,2])

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

#Modified Julia Date
MJD_0 = 58155.0 #modified julian day

#magnetic field generation for full orbit
t0 = 0.0 #seconds
tf = 3T #seconds

#generate time vector
dt = .2
t = t0:dt:tf
N = length(t)

#extrapolate over the magnetic field
alt = 400 # km
mag_field = igrf_data(alt,2019)

#find the ECI magnetic field and the position and velocity
B_ECI,pos,vel = magnetic_simulation(kep,t0,tf,N,mag_field,GM,MJD_0)
B_ECI = B_ECI'

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
for i = 1:N
    ROT[:,:,i] = [r[:,i] v[:,i] h[:,i]]'
    q_TF[:,i] = rot_quaternion(ROT[:,:,i])
end

#plot magnetic field
plot(B_ECI')

#initial
omega_int=[0;0;0]
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

#calculate the orbital mean motion
ω_0 = sqrt(GM/kep[2]^3) #rad/s

#initialize state
euler_angles = zeros(3,length(t))
x = zeros(7,length(t))
ω_guess, q_guess = eigen_axis_slew([omega_int;q_0],[omega_final;q_final],t)
x = [ω_guess;q_guess]
# x[:,1] = [omega_int;qmult(q_inv(q_0),q_TF[:,1])]
euler_angles[:,1] = quaternion_euler(x[4:7,1])

#extrapolate state once to get higher order magnetic field info
ẋ = attitude_dynamics(x[:,1],zeros(3),B_ECI[:,1],J)
x[:,2] = x[:,1] + dt * ẋ
euler_angles[:,2] = quaternion_euler(x[4:7,2])

#assign weights
C_1 = 1E-6
C_2 = 1E-9

#magnetic moment limit
m_limit = 0.19

# create an error quaternion matrix
q̄ = zeros(4,length(t))

#simulate
for i = 2:length(t)-1

    #find the measured magnetic field
    # B_meas=qrot(q_inv(x[4:7,i]),B_ECI[:,i])
    B_meas=qrot(q_inv(x[4:7,i]),B_ECI[:,i])

    #find the difference in omega
    # ω̄ = qrot(q_inv(x[4:7,i]),(x[1:3,i]-ω_0*h[:,i]))
    ω̄ = ω_guess[:,i] - x[1:3,i]

    #find the difference in quaternion
    # q̄[:,i] = qmult(x[4:7,i],q_inv(q_TF[:,i]))
    # q̄[:,i] = qmult(q_inv(q_final_permuted[:,i]),x[4:7,i]) #finds quaternion error
    q̄[:,i] = (qmult(x[4:7,i],q_guess[:,i])) #finds quaternion error


    #control law
    m_limit = 0.19
    m̄ = psiaki_controller(C_1, C_2, J, q̄[:,i], ω̄, B_meas, m_limit)

    #extrapolate state
    ẋ = attitude_dynamics(x[:,i],m̄,B_meas,J)
    x[:,i+1] = rk4_psiaki(attitude_dynamics,x[:,i],dt,m̄,B_meas,J)
    x[4:7,i+1]=normalize(x[4:7,i+1])

end

for i = 1:N
    euler_angles[:,i] = quaternion_euler(x[4:7,i])
end

plot(x[1:3,1:N]')
plot(x[4:7,1:N]')
plot(q̄')
plot(euler_angles[:,1:7500]')

plt = plot(t[1:100],x[4:7,1:100]',label = ["q1" "q2" "q3" "q4"])
plot!(title = "PD Control", xlabel = "Time (s)", ylabel = "Quaternion")

PD_error = zeros(length(t))
for i = 1:length(t)-1
    PD_error[i] = 2*acos(min(1,qmult(q_inv(x[4:7,i]),q_final)[1]))-pi
end

plot(t[1:100],PD_error[1:100])
