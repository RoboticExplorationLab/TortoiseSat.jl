#this code is for a Monte Carlo simulation of the TortoiseSat problem
using Pkg
using LinearAlgebra
using Plots
using TrajectoryOptimization
using DifferentialEquations
using SatelliteToolbox
using ForwardDiff
using SparseArrays
using Interpolations
using Statistics
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

number_sims = 100

N = 5000 #number of knot points
n=8
m=3#assign state and control variables
t0=0.0 #seconds

# states_temp = states
# control_inputs_temp = control_inputs
# sim_states_temp = sim_states
# sim_control_inputs_temp = sim_control_inputs
# B_ECI_total_temp = B_ECI_total
# A_temp = h5read("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/100_A.h5","100_A")

#create the arrays that will hold the simulation data
states = Array{Float64,2}[]
control_inputs = Array{Float64,2}[]
sim_states = Array{Float64,2}[]
sim_control_inputs = Array{Float64,2}[]
B_N_initial = Array{Float64,2}[]
B_ECI_total = Array{Float64,2}[]
B_gram = Array{Float64,2}[]

#create keplerian element array
A=zeros(number_sims,6)

#create array for slew time
time_step = ones(1,number_sims)
t_total = Array{Float64,1}[]
t_final = zeros(number_sims)

#define limits of slew time calculation
slew_limits = zeros(1,2)
slew_limits[1] = .1 #rad/s
slew_limits[2] = .1 #error quaternion

#initialize number of knot points and MJD_0
MJD_0 = 58155.0 #modified julian day

t0 = 0.0 #seconds
tf = 60*40 #seconds
cutoff = 100 #condition number cutoff

#extrapolate over magnetic field
alt = 400 # km
mag_field = igrf_data(alt,2019)

#plot magnetic field
B_N_complete = zeros(3,1000,1000)
B_N_complete_norm = zeros(1000,1000)
lat_complete = range(-pi/2,length = 1000, stop = pi/2)
long_complete = range(-pi, length = 1000, stop = pi)
for i = 1:1000
    for j = 1:1000
        B_N_complete[1,i,j] = mag_field(i,j,1)
        B_N_complete[2,i,j] = mag_field(i,j,2)
        B_N_complete[3,i,j] = mag_field(i,j,3)
        B_N_complete_norm[i,j] = norm(B_N_complete[:,i,j])
    end
end
plt = plot(long_complete,lat_complete,B_N_complete_norm*1E9)
title!("Normalized Magnetic Field Strength (nT)")
xlabel!("Longitude (rad)")
ylabel!("Latitude (rad)")
Plots.savefig("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/mag_field_complete.png")

#trajectory optimization section
#initial
omega_int=[0;0;0]
q_0=[0.;0;1;0]
x0=[omega_int;q_0;t0]

#final
q_final=[sqrt(2)/2;sqrt(2)/2;0;0]
omega_final=[0.;0;0]
xf=[omega_final;q_final;1] #km or km/s

for i in retry
    #enter orbital parameters
    alt=400 #altitude in km
    #enter orbital parameters
    A[i,1] = 0 #circular orbits
    A[i,2]= alt+R_E #semi-major axis
    A[i,3] = rand(1)[1]*90 #inclination in degrees
    A[i,4] = rand(1)[1]*360 #assign ascending nodes
    A[i,5] = 0 #argument of periapsis
    A[i,6] = rand(1)[1]*360 #initial true anomoly
    # T=2*pi*sqrt(A[2].^3/GM)

    #Modified Julia Date
    MJD_0 = 58155.0 #modified julian day

    #extrapolate over the magnetic field
    push!(B_N_initial,magnetic_simulation(A[i,:],t0,tf,N,mag_field,GM,MJD_0)[1])

    # Find magnetic control garmian
    B_gram = magnetic_gramian(B_N_initial[i],(tf-t0)/N)

    #determine at what point gramian eigenvalues meet criteria
    t_final[i] = condition_based_time(B_gram,cutoff)*(tf-t0)/N
    # while (t_final[i]-t0)/N > .2
    #     t_final[i] *= .99
    # end
    push!(t_total,range(t0,step = .2, stop = t_final[i]))
    time_step[i] = t_total[i][end]-t_total[i][end-1]

    #reasign time
    # find magnetic field along new time
    push!(B_ECI_total,magnetic_simulation(A[i,:],t0,t_final[i], N, mag_field,GM,MJD_0)[1])
    B_ECI = B_ECI_total[i]
    global B_ECI

    #create model
    n=8
    m=3

    #consolidate the model
    model=Model(DerivFunction,n,m,quaternion_error,quaternion_expansion)

    #use eigen-axis slew as first guess for state
    ω_guess, q_guess = eigen_axis_slew(x0[1:7],xf[1:7],t_total[i])
    X = [ω_guess;q_guess;reshape(t_total[i],1,length(t_total[i]))]

    #Bryson's Rule LQR weights
    Q = zeros(n,n)
    Qf = zeros(n,n)
    ω_max = maximum(X[1:3,:])
    τ_max = maximum(J * diff(X[1:3,:],dims = 2)/time_step[i])
    m_max = τ_max / 1.e-5 * 1.e2
    α = 1.e-1
    Q[1:3,1:3] = Array((α/ω_max^2)*Diagonal(I,3))
    Qf[1:3,1:3] = Array((α/ω_max^2)*Diagonal(I,3))*10
    β = 1.e3
    Q[4:7,4:7] = Array(α*β*Diagonal(I,4))
    Qf[4:7,4:7] = Array(α*β*Diagonal(I,4))*10
    R = Array((1/m_max^2)*.1*Diagonal(I,m))

    # #LQR weights
    # Q=zeros(n,n)
    # Qf=zeros(n,n)
    # α = 1.e1
    # Q[1:3,1:3] = Array((α)*Diagonal(I,3))
    # Qf[1:3,1:3] = Array((α)*Diagonal(I,3))*100
    # β = 1.e1
    # Q[4:7,4:7] = Array(β*Diagonal(I,4))
    # Qf[4:7,4:7] = Array(β*Diagonal(I,4))*100
    #
    # R = Array((7.5e-4*t_final[i])*Diagonal(I,m))

    #bounds
    u_bnd=[19;19;19]
    x_bnd=10

    #random start to trajopt
    # cost = GenericCost(quat_cost,quat_cost_terminal,quat_cost_grad,quat_cost_hess,n,m)
    # obj = UnconstrainedObjective(cost,tf,x0,xf)
    obj = LQRObjective(Q, R, Qf, t_final[i], x0, xf)
    obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
    solver = TrajectoryOptimization.Solver(model,obj_con,dt = time_step[i])
    solver.opts.verbose= false
    solver.opts.iterations_outerloop=3
    solver.opts.iterations_innerloop=50
    solver.opts.dJ_counter_limit = 1
    solver.opts.sat_att = true
    U = rand(Float64,m,length(t_total[i]))/1000 #initialize random input and state
    # X = rand(Float64,m,solver.N)/1000

    results, stats = TrajectoryOptimization.solve(solver,U)
    X = TrajectoryOptimization.to_array(results.X) #assigns state and input
    U = TrajectoryOptimization.to_array(results.U)

    push!(states,X)
    push!(control_inputs,U)

    #TVLQR
    x0_lqr = zeros(n)
    x0_lqr[1:3] = x0[1:3]

    q_noise = randn(3,1)*(.0001*pi/180)^2 #assign noise to initial quaternion
    θ_noise = norm(q_noise)
    r_noise = q_noise/θ_noise
    x0_lqr[4:7] = qmult(x0[4:7],[cos(θ_noise/2);r_noise*sin(θ_noise/2)])
    integration=:rk4

    dt_lqr = .2

    include("simulator.jl")
    f! = simulator
    f_gains! = gain_simulator

    α = 1.e1
    Q_lqr=zeros(6,6)
    Qf_lqr=zeros(6,6)
    Q_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))
    Qf_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))*100
    β = 1.e1
    Q_lqr[4:6,4:6] = Array(β*Diagonal(I,3))
    Qf_lqr[4:6,4:6] = Array(β*Diagonal(I,3))*100

    R_lqr = Array((7.5e3)*Diagonal(I,m))

    X_sim, U_sim, dX, K = attitude_simulation(f!,f_gains!,integration,X,U,dt_lqr,x0_lqr,t0,t_final[i],Q_lqr,R_lqr,Qf_lqr)

    push!(sim_states,X_sim)
    push!(sim_control_inputs,U_sim)
    println("There have been ",i," iterations")
end

slew_time = zeros(number_sims) .+ t_final
omega_norm = []
error_quat = []
fails = zeros(number_sims)

for i = 1:number_sims
    #figure out how much time is required to turn
    omega_norm_vec = zeros(size(sim_states[i],2))
    error_quat_vec = zeros(4,size(sim_states[i],2))
    for j = 1:size(sim_states[i],2)
        omega_norm_vec[j] = norm(sim_states[i][1:3,i])
        q_final=[sqrt(2)/2;sqrt(2)/2;0;0]
        error_quat_vec[:,j] = qmult(q_inv(q_final),sim_states[i][4:7,j])
        error_angle = 2*acos(min((error_quat_vec[1,j]),1.))
        if j>10 && omega_norm_vec[j] < slew_limits[1] && error_angle < slew_limits[2] && slew_time[i] == t_final[i]
            slew_time[i] = time_step[i]*j
        end
    end
    push!(omega_norm,omega_norm_vec)
    push!(error_quat,error_quat_vec)

    if slew_time[i] == t_final[i]
        fails[i] = 1
    else fails[i] = 0
    end
end

slew_time
t_final
fails
retry = findall(x -> x == 1.,fails)
for i in retry
    println(i)
end

# construct data for the heatmap
z = zeros(40,26)
for i = 1:number_sims
    x_segment = convert(Int,ceil(t_final[i]/50))
    y_segment = convert(Int,ceil(slew_time[i]/50))
    z[x_segment,y_segment] += 1
end

xs = [string(i) for i = 0:100:1900]
ys = [string(i) for i = 0:100:1200]
z
heatmap(1:40, 1:26, z')
xlabel!("Magnetic Cutoff (s)")
ylabel!("Slew Time (s)")

scatter(t_final,slew_time,legend = false)
xlabel!("Cutoff (s)")
ylabel!("Slew Time (s)")

#plot inputs and make sure they don't match
k = 98
plot(t_total[k][1:size(control_inputs[k],2)],control_inputs[k]')
plot(t_total[k][1:size(sim_control_inputs[k],2)],sim_control_inputs[k]')

plot(t_total[k][1:size(sim_states[k],2)],sim_states[k][4:7,:]')


#plot B_N
k = 292
plot(t_total[1:N,k],B_ECI_total[1:N,:,k])

#plot omega
k = 292
plot(t_total[1:N,k],sim_states[1:3,:,k]')


#plot quaternion agains the sim quaterion to make sure TVLQR is working
k = 167
plot(t_total[1:N,k],sim_states[4:7,:,k]')

plot(sim_states[4,:,2])
plot!(sim_states[5,:,2])
plot!(sim_states[6,:,2])
plot!(sim_states[7,:,2])

plot(states[4,:,2])
plot!(states[5,:,2])
plot!(states[6,:,2])
plot!(states[7,:,2])

histogram(slew_time[101:200], bins = 50, legend = false)
xaxis!("Slew Time (s)")
yaxis!("# of tests")
savefig("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/100_tests.png")

slew_time_mean = 0.0
a = 0.0
for i = 1:number_sims
    if slew_time[i] > 0.0
        slew_time_mean += slew_time[i]
        a += 1.0
    end
    global slew_time_mean
    global a
end

slew_time_mean/a

function istf(x)
    if x == t
        true
    else
        false
    end
end
fails = findall(istf,slew_time)

using HDF5
using CSV
h5write(string("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/100_A.h5"), "A", A)
for i in 1:100
    h5write(string("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/100_states_",i,".h5"),"one_state", sim_states[i])
    h5write(string("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/100_states_",i,".h5"), "states", sim_states[i])
    h5write(string("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/100_control_",i,".h5"), "control", sim_control_inputs[i])
    h5write(string("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/100_B_N_",i,".h5"), "B_ECI", B_ECI_total[i])
    h5write(string("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/100_t_total_",i,".h5"), "t_total", t_total[i])
end
