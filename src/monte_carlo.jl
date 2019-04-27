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


number_sims = 1000
N=10000 #initialize the ultimate number of knot points
n=8
m=3#assign state and control variables

#create the arrays that will hold the simulation data
states = zeros(n,N,number_sims)
control_inputs = zeros(m,N-1,number_sims)
sim_states = zeros(n,N,number_sims)
sim_control_inputs = zeros(m,N,number_sims)
B_N_total = zeros(2N+1,3,number_sims)

#create keplerian element array
A=zeros(number_sims,6)

#time initialization
MJD_0 = 58155.0 #modified julian day
t0=0.0 #seconds
tf=60.0*40+t0 #seconds
dt=(tf-t0)/(N)

#create array for slew time
slew_time = ones(1,number_sims)*tf

#initial
omega_int=[0;0;0]
q_0=[0.;0;1;0]
x0=[omega_int;q_0;t0]

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

for i = 1:number_sims

    #enter orbital parameters
    A[i,1] = 0 #circular orbits
    A[i,2]= alt+R_E #semi-major axis
    A[i,3] = rand(1)[1]*180 #inclination in degrees
    A[i,4] = rand(1)[1]*360 #assign ascending nodes
    A[i,5] = 0 #argument of periapsis
    A[i,6] = rand(1)[1]*360 #initial true anomoly
    T=2*pi*sqrt(A[i,2].^3/GM)

    #initial
    pos_0=kep_ECI(A[i,:],t0,GM) #km and km/s
    u0=[pos_0[1,:];pos_0[2,:]] #km or km/s

    # find magnetic field along 1 orbit
    tspan=(t0,tf*2) #go a little over for the ForwardDiff downstream
    prob=DiffEqBase.ODEProblem(OrbitPlotter,u0,tspan)
    sol=DiffEqBase.solve(prob,dt=dt,adaptive=false,RK4())
    pos=sol[1:3,:]
    vel=sol[4:6,:]

    #convert from ECI to ECEF for IGRF
    t = [t0:dt:2*tf...]
    GMST = (280.4606 .+ 360.9856473*(t/24/60/60 .+ (MJD_0.+t./(24*60*60)) .- 51544.5))./180*pi #calculates Greenwich Mean Standard Time
    ROT = zeros(3,3,length(t)) #cretes rotation matrix for every time step
    pos_ecef = zeros(3,length(t)) #populates ecef position vector
    lat = zeros(length(t),1)
    long = zeros(length(t),1) # creates latitude and longitude
    for i = 1:length(t)-1
        ROT[:,:,i] = [cos(GMST[i]) sin(GMST[i]) 0;
        -sin(GMST[i]) cos(GMST[i]) 0;
        0 0 1]
        pos_ecef[:,i] = ROT[:,:,i]*pos[:,i]
        lat[i] = pos_ecef[3,i]/norm(pos_ecef[:,i])
        long[i] = atan(pos_ecef[2,i],pos_ecef[1,i])
    end

    #find the IGRF magnetic field for the orbit!
    B_N_sim = zeros(size(pos,2)-1,3)
    for i = 1:size(pos,2)-1
        if (long[i]+π)/(2*π)*size(mag_field,1) < 1
            long[i] = 1/size(mag_field,1)*2*π-π
        end
        B_N_sim[i,1] = mag_field((lat[i]+π/2)/(π)*size(mag_field,1),(long[i]+π)/(2*π)*size(mag_field,1),1)
        B_N_sim[i,2] = mag_field((lat[i]+π/2)/(π)*size(mag_field,1),(long[i]+π)/(2*π)*size(mag_field,1),2)
        B_N_sim[i,3] = mag_field((lat[i]+π/2)/(π)*size(mag_field,1),(long[i]+π)/(2*π)*size(mag_field,1),3)
    end
    B_N_sim = convert(Array{Float64},B_N_sim)
    global B_N_sim
    B_N_total[:,:,i] = B_N_sim

    #create model
    n=8;
    m=3;

    #define expansion function
    function quaternion_expansion(cost::QuadraticCost, x::AbstractVector{Float64}, u::AbstractVector{Float64})
        #renormalize for quaternion
        qk = x
        sk = qk[4]
        vk = qk[5:7]

        qn = x
        sn = qn[4]
        vn = qn[5:7]

        Gk = [-vk'; sk*Array(1.0*Diagonal(I,3)) + hat(vk)]
        Gn = [-vn'; sn*Array(1.0*Diagonal(I,3)) + hat(vn)]

        #recompute A and B matrices through permutation matrices
        perm_Gn = zeros(7,8)
        perm_Gk = zeros(8,7)
        perm_Gn[1:3,1:3]= Array(1.0*Diagonal(I,3))
        perm_Gn[4:6,4:7] = Gn'
        perm_Gk[1:3,1:3]= Array(1.0*Diagonal(I,3))
        perm_Gk[4:7,4:6] = Gk
        return perm_Gn*cost.Q*perm_Gk, cost.R, cost.H*perm_Gk, perm_Gn*(cost.Q*x+cost.q), cost.R*u+ cost.r
    end

    function quaternion_expansion(cost::QuadraticCost, xN::AbstractVector{Float64})
        #renormalize for quaternion at boundary
        qn = xN
        sn = qn[4]
        vn = qn[5:7]

        Gn = [-vn'; sn*Array(1.0*Diagonal(I,3)) + hat(vn)]
        perm_Gn = zeros(7,8)
        perm_Gn[1:3,1:3]= Array(1.0*Diagonal(I,3))
        perm_Gn[4:6,4:7] = Gn'

        return perm_Gn*cost.Qf*perm_Gn', perm_Gn*(cost.Qf*xN+cost.qf)
    end

    function quaternion_error(X1,X2)
    #for a state where:
    #1-3: omega
    #4-7: quaternion
    #8: time
        δx = zeros(7)
        δx[1:3] = X1[1:3] - X2[1:3]

        # #calculate the error quaternion
        # δx[4:6] = [zeros(3,1) Array(1.0*Diagonal(I,3))]*qmult(q_inv(X2[4:7]),X1[4:7])

        #calculate the modified Rodrigues parameters and find the error
        q_e = qmult(q_inv(X2[4:7]),X1[4:7]) #finds quaternion error
        MRP = q_e[2:4]/(1+q_e[1]) #finds Modified Rodrigues Parameter
        δx[4:6] = MRP

        return δx
    end

    model=Model(DerivFunction,n,m,quaternion_error,quaternion_expansion);

    #LQR
    Q=zeros(n,n)
    Qf=zeros(n,n)
    α = 1.e0
    Q[1:3,1:3] = Array((α)*Diagonal(I,3))
    Qf[1:3,1:3] = Array((α)*Diagonal(I,3))
    β = 1.e0
    Q[4:7,4:7] = Array(β*Diagonal(I,4))
    Qf[4:7,4:7] = Array(β*Diagonal(I,4))*1.e0

    R = Array((5.e-3/norm(B_N_sim))*Diagonal(I,m))

    #bounds
    u_bnd=1.9
    x_bnd=1

    #non warm start (cold start?)
        obj = LQRObjective(Q, R, Qf, tf, x0, xf)
        obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
        solver = TrajectoryOptimization.Solver(model,obj,N=N)
        solver.opts.verbose=false
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

    #figure out how much time is required to turn
    omega_norm = zeros(1,N)
    for j = 1:N
        omega_norm[j] = norm(sim_states[1:3,j,i])
        error_quat = [zeros(3,1) Array(1.0*Diagonal(I,3))]*qmult(q_inv(q_final),sim_states[4:7,j,i])
        if j>10 && norm(omega_norm[j-10:j],1) < slew_limits[1]*2
            slew_time[i] = dt*j
            break
        end
    end
    print("There have been ",i," iterations")

end

#plot inputs and make sure they don't match
plot(control_inputs[1,:,1])
plot!(control_inputs[2,:,1])
plot!(control_inputs[3,:,1])

plot(control_inputs[1,:,2])
plot!(control_inputs[2,:,2])
plot!(control_inputs[3,:,2])

plot(control_inputs[1,:,3])
plot!(control_inputs[2,:,3])
plot!(control_inputs[3,:,3])

plot(control_inputs[1,:,4])
plot!(control_inputs[2,:,4])
plot!(control_inputs[3,:,4])

plot(control_inputs[1,:,5])
plot!(control_inputs[2,:,5])
plot!(control_inputs[3,:,5])

k = 100
plot(control_inputs[1:3,:,k]')

plot(B_N_total[:,:,k])

#plot omega

k = 100
plot(states[1,:,k])
plot!(states[2,:,k])
plot!(states[3,:,k])

plot(sim_states[1,:,1])
plot!(sim_states[2,:,1])
plot!(sim_states[3,:,1])

plot(sim_states[1,:,2])
plot!(sim_states[2,:,2])
plot!(sim_states[3,:,2])

plot(sim_states[1,:,3])
plot!(sim_states[2,:,3])
plot!(sim_states[3,:,3])

plot(sim_states[1,:,4])
plot!(sim_states[2,:,4])
plot!(sim_states[3,:,4])

plot(sim_states[1,:,5])
plot!(sim_states[2,:,5])
plot!(sim_states[3,:,5])



#plot quaternion agains the sim quaterion to make sure TVLQR is working
k = 1
plot(states[4,:,k])
plot!(states[5,:,k])
plot!(states[6,:,k])
plot!(states[7,:,k])

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
