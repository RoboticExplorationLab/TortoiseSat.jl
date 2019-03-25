#this code is for a Monte Carlo simulation of the TortoiseSat problem
using Pkg;
using LinearAlgebra;
using Plots;
Pkg.add("Formatting")
Pkg.add("TrajectoryOptimization")
using TrajectoryOptimization;
Pkg.add("DifferentialEquations");
using DifferentialEquations;
Pkg.add("SatelliteToolbox");
using SatelliteToolbox
Pkg.add("ForwardDiff")
using ForwardDiff
Pkg.add("SparseArrays")
using SparseArrays
Pkg.add("Interpolations")
using Interpolations
Pkg.add("Statistics")
using Statistics
Pkg.add("AttitudeController")
using AttitudeController

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

number_sims = 100
N=100 #initialize the ultimate number of knot points
n=8
m=3#assign state and control variables

#create the arrays that will hold the simulation data
states = zeros(n,N,number_sims)
control_inputs = zeros(m,N-1,number_sims)
sim_states = zeros(n,N,number_sims)
sim_control_inputs = zeros(m,N-1,number_sims)
B_N_total = zeros(N*2,3,number_sims)

#create keplerian element array
A=zeros(number_sims,6)

#time initialization
MJD_0 = 58155.0 #modified julian day
t0=0.0 #seconds
tf=60.0*5+t0 #seconds
dt=(tf-t0)/(N)

#create array for slew time
slew_time = ones(1,number_sims)*tf

#initial
omega_int=[0;0;0]
q_0=[0.;0;1;0]
x0=[omega_int;q_0;t0]

#final
q_final=[0.;1;0;0]
omega_final=[0.;0;0]
xf=[omega_final;q_final;1] #km or km/s

#define limits of slew time calculation
slew_limits = zeros(1,2)
slew_limits[1] = .01 #rad/s
slew_limits[2] = .1 #error quaternion

for i = 1:number_sims

    #enter orbital parameters
    alt=400 #altitude in km
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
    t = [t0:dt:tf+dt*100...]
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
        long[i] = atan(pos_ecef[2,i]/pos_ecef[1,i])
    end

    #find the IGRF magnetic field for the orbit!
    #magnetic field (IGRF)
    B_N_sim = zeros(size(pos,2)-1,3)
    for i = 1:size(pos,2)-1
        B_N_sim[i,:] = SatelliteToolbox.igrf12(2019,norm(pos_ecef[:,i])*1000,lat[i],long[i])/1.e9
    end
    B_N_sim = convert(Array{Float64},B_N_sim)
    global B_N_sim
    B_N_total[:,:,i] = B_N_sim

    #trajectory optimization section

    #enter into warm start
    N_vect = N/5:N/5:N #create initial number of knot points
    dt_vect=(tf-t0)./(N_vect)

    #create model
    n=8;
    m=3;
    model=Model(DerivFunction,n,m);

    #LQR
    Q=zeros(n,n)
    Qf=zeros(n,n)
    Q[1:3,1:3] = Array((1.e-1)*Diagonal(I,3))
    Qf[1:3,1:3] = Array((1.e-3)*Diagonal(I,3))
    α = 3.e3
    Q[4:7,4:7] = Array(α*Diagonal(I,4))
    Qf[4:7,4:7] = Array(α*Diagonal(I,4))*1.e3

    R = Array((100*norm(B_N_sim))*Diagonal(I,m))

    #bounds
    u_bnd=19
    x_bnd=10

    #define the quaternion cost function
    function quat_cost(x,u,xf)
        error_quat = [zeros(3,1) Array(1.0*Diagonal(I,3))]*qmult(q_inv(xf[4:7]),x[4:7])
        0.5*(x[1:3] - xf[1:3])'*Q[1:3,1:3]*(x[1:3] - xf[1:3]) + 0.5*u'*R*u + Q[4,4]*(error_quat'*error_quat) + c
    end

    #non warm start (cold start?)
        obj = LQRObjective(Q, R, Qf, tf, x0, xf)
        obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
        solver = TrajectoryOptimization.Solver(model,obj_con,N=N)
        solver.opts.verbose=true
        solver.opts.iterations_outerloop=5
        solver.opts.iterations_innerloop=50
        # solver.opts.sat_att=true
        solver.opts.dJ_counter_limit = 1
            U = rand(Float64,m,solver.N)/1000 #initialize random input and state
            X = rand(Float64,m,solver.N)/1000
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

    sim_states[:,:,i] = X
    sim_control_inputs[:,:,i] = U

    #since we are no longer interpolating, reassign
    X_lqr = X
    U_lqr = U
    dt_lqr = dt

    # #assigns simulation and gains functions for TVLQR
    # include("simulator.jl")
    # f! = simulator
    # f_gains! = gain_simulator
    #
    # #creates Q and R matrices for TVLQR
    # α = 1.e-1
    # β = 1.e3
    # Q_lqr = zeros(6,6)
    # Q_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))
    # Q_lqr[4:6,4:6] = Array((α*β)*Diagonal(I,3))
    # Qf_lqr = zeros(6,6)
    # Qf_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))*1.e4
    # Qf_lqr[4:6,4:6] = Array((α*β)*Diagonal(I,3))*1.e8
    # R_lqr = Array((2.e-2)*Diagonal(I,3))
    #
    # #executes TVLQR
    # X_sim, U_sim, dX, K = attitude_simulation(f!,f_gains!,integration,X_lqr,U_lqr,dt_lqr,x0_lqr,t0,tf,Q_lqr,R_lqr,Qf_lqr)
    #
    # #assigns simulations for external use
    # sim_states[:,:,i] = X_sim
    # sim_control_inputs[:,:,i] = U_sim
    #
    #figure out how much time is required to turn
    omega_norm = zeros(1,N)
    for j = 1:N
        omega_norm[j] = norm(sim_states[1:3,j,i])
        error_quat = [zeros(3,1) Array(1.0*Diagonal(I,3))]*qmult(q_inv(q_final),sim_states[4:7,j,i])
        if j>10 && norm(omega_norm[j-10:j],1) < slew_limits[1] && norm(error_quat) < slew_limits[2]
            slew_time[i] = dt*j
            break
        end
    end

end

# for i = 1:number_sims
#     omega_norm = zeros(1,N)
#     for j = 1:N
#         omega_norm[j] = norm(sim_states[1:3,j,i])
#         if j>10 && norm(omega_norm[j-10:j],1) < .01
#             slew_time[i] = dt*j
#             break
#         end
#     end
# end

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

#plot omega
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
plot(states[4,:,1])
plot!(states[5,:,1])
plot!(states[6,:,1])
plot!(states[7,:,1])

plot(sim_states[4,:,2])
plot!(sim_states[5,:,2])
plot!(sim_states[6,:,2])
plot!(sim_states[7,:,2])

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
