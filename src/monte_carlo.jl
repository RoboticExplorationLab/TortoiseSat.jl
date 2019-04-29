#this code is for a Monte Carlo simulation of the TortoiseSat problem
using Pkg;
using LinearAlgebra;
using Plots;
using TrajectoryOptimization;
using DifferentialEquations;
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


number_sims = 2
N=1000 #initialize the ultimate number of knot points
n=8
m=3#assign state and control variables

#create the arrays that will hold the simulation data
states = zeros(n,N,number_sims)
control_inputs = zeros(m,N-1,number_sims)
sim_states = zeros(n,N,number_sims)
sim_control_inputs = zeros(m,N,number_sims)
B_N_total = zeros(2N,3,number_sims)

#create keplerian element array
A=zeros(number_sims,6)

#create array for slew time
time_step = ones(1,number_sims)
slew_time = ones(1,number_sims)
t_total = zeros(N,number_sims)

#initial
omega_int=[0;0;0]
q_0=[0.;0;1;0]
x0=[omega_int;q_0;0]

#final
q_final=[sqrt(2)/2;sqrt(2)/2;0;0]
omega_final=[0.;0;0]
xf=[omega_final;q_final;1] #km or km/s

#define limits of slew time calculation
slew_limits = zeros(1,2)
slew_limits[1] = .001 #rad/s
slew_limits[2] = .1 #error quaternion

#extrapolate over the magnetic field
alt = 400 # km
mag_field = igrf_data(alt,2019)

#initialize number of knot points and MJD_0
MJD_0 = 58155.0 #modified julian day
t0 = 0.0

let N = N, mag_field = mag_field
    for i = 1:number_sims

        #time initialization
        t0=0.0 #seconds
        tf=60.0*40+t0 #seconds
        dt=(tf-t0)/(N)

        #enter orbital parameters
        A[i,1] = 0 #circular orbits
        A[i,2]= alt+R_E #semi-major axis
        A[i,3] = rand(1)[1]*180 #inclination in degrees
        A[i,4] = rand(1)[1]*360 #assign ascending nodes
        A[i,5] = 0 #argument of periapsis
        A[i,6] = rand(1)[1]*360 #initial true anomoly
        T=2*pi*sqrt(A[i,2].^3/GM)

        #Modified Julia Date
        MJD_0 = 58155.0 #modified julian day

        #magnetic field generation for half orbit
        t0 = 0.0 #seconds
        tf = 60*40 #seconds
        cutoff = 20 #condition number cutoff
        N = 1000
        B_N_sim = magnetic_simulation(A[i,:],t0,tf,N)

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

        #assignm magnetic field outside
        B_N_total[:,:,i] = B_N_sim

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

        R = Array((4.e-3/norm(B_N_sim))*Diagonal(I,m))

        #bounds
        u_bnd=[19;19;19]
        x_bnd=10

        #non warm start (cold start?)
            obj = LQRObjective(Q, R, Qf, tf, x0, xf)
            obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
            solver = TrajectoryOptimization.Solver(model,obj,N=N)
            solver.opts.verbose=true
            solver.opts.iterations_outerloop=3
            solver.opts.iterations_innerloop=100
            solver.opts.cost_tolerance = 1.e-2
            solver.opts.sat_att=true
            solver.opts.dJ_counter_limit = 1
                U = rand(Float64,m,solver.N)/1000 #initialize random input and state
                X = rand(Float64,m,solver.N)/1000

            #use eigen-axis slew as first guess for state
            ω_guess,  q_guess = eigen_axis_slew(x0[1:7],xf[1:7],t[1:N])
            X = [ω_guess[:,1:N];q_guess[:,1:N];reshape(t[1:N],1,N)]

            results, stats = TrajectoryOptimization.solve(solver,U)
            X = TrajectoryOptimization.to_array(results.X) #assigns state and input
            U = TrajectoryOptimization.to_array(results.U)

            #assigns trajectories for external use
            states[:,:,i] = X
            control_inputs[:,:,i] = U

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

        # sim_states[:,:,i] = X
        # sim_control_inputs[:,:,i] = U

        #since we are no longer interpolating, reassign
        X_lqr = X
        U_lqr = U
        dt_lqr = dt

        #assigns simulation and gains functions for TVLQR
        include("simulator.jl")
        f! = simulator
        f_gains! = gain_simulator

        #creates Q and R matrices for TVLQR
        α = 1.e0
        β = 1.e1
        Q_lqr = zeros(6,6)
        Q_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))
        Q_lqr[4:6,4:6] = Array((α*β)*Diagonal(I,3))
        Qf_lqr = zeros(6,6)
        Qf_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))
        Qf_lqr[4:6,4:6] = Array((α*β)*Diagonal(I,3))
        R_lqr = Array((1.e0)/norm(B_N_sim)*Diagonal(I,3))

        #executes TVLQR
        X_sim, U_sim, dX, K = attitude_simulation(f!,f_gains!,integration,X_lqr,U_lqr,dt_lqr,x0_lqr,t0,tf,Q_lqr,R_lqr,Qf_lqr)

        #assigns simulations for external use
        sim_states[:,:,i] = X_sim
        sim_control_inputs[:,:,i] = U_sim
        print("There have been ",i," iterations")

    end
end

for i = 1:number_sims

    #figure out how much time is required to turn
    omega_norm = zeros(1,N)
    for j = 1:N
        omega_norm[j] = norm(sim_states[1:3,j,i])
        error_quat = [zeros(3,1) Array(1.0*Diagonal(I,3))]*qmult(q_inv(q_final),sim_states[4:7,j,i])
        if j>10 && norm(omega_norm[j-10:j],1) < slew_limits[1]*2
            slew_time[i] = time_step[i]*j
            break
        end
    end
    if slew_time[i] <= 1.0
        slew_time[i] = t_total[end,i]
    end


end

#plot inputs and make sure they don't match
plot(t_total[1:N-1,1],control_inputs[1:3,:,1]')

plot(t_total[1:N-1,2],control_inputs[1:3,:,2]')

k = 100
plot(t_total[1:N-1,k],control_inputs[1:3,:,k]')
plot(B_N_total[:,:,k])

#plot omega

k = 1
plot(t_total[1:N,k],states[1:3,:,k]')


#plot quaternion agains the sim quaterion to make sure TVLQR is working
k = 1
plot(t_total[1:N,k],states[4:7,:,k]')

plot(sim_states[4,:,2])
plot!(sim_states[5,:,2])
plot!(sim_states[6,:,2])
plot!(sim_states[7,:,2])

plot(states[4,:,2])
plot!(states[5,:,2])
plot!(states[6,:,2])
plot!(states[7,:,2])

histogram(slew_time[1:100], bins = 50)

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
    if x == tf
        true
    else
        false
    end
end
fails = findall(istf,slew_time)
