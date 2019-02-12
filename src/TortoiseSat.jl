#this code is for trajectory optimization of an underactuated satellite

using Pkg;
using LinearAlgebra;
using Plots;
Pkg.add("Formatting")
Pkg.develop(PackageSpec(url="https://github.com/GathererA/TrajectoryOptimization.jl"))
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

#satellite parameters
mass=.75 #kg
J=[00.00125 0 0;
    0 0.00125 0;
    0 0 0.00125]
J_inv=inv(J)
BC = 50 #kg/m^2

#control parameters
m_max=1E-4

#enter orbital parameters
alt=400 #altitude in km
A=zeros(1,6)
A[1] = 0 #circular orbits
A[2]= alt+R_E #semi-major axis
A[3] = 51.64 #inclination in degrees
A[4] = 0 #assign ascending nodes
A[5] = 0 #argument of periapsis
A[6] = 0 #initial true anomoly
T=2*pi*sqrt(A[2].^3/GM)

#time initialization
MJD_0 = 58155.0 #modified julian day
t0=0 #seconds
tf=60.0*10+t0 #seconds

#initial
pos_0=kep_ECI(A,t0,GM) #km and km/s
u0=[pos_0[1,:];pos_0[2,:]] #km or km/s

# find magnetic field along 1 orbit
N=5000 #initialize the ultimate number of knot points
dt=(tf-t0)/(N)
tspan=(t0,tf+dt*99) #go a little over for the ForwardDiff downstream
prob=DiffEqBase.ODEProblem(OrbitPlotter,u0,tspan)
sol=DiffEqBase.solve(prob,dt=dt,adaptive=false,RK4())
pos=sol[1:3,:]
vel=sol[4:6,:]
plot(sol.t,pos[1,:])
plot!(sol.t,pos[2,:])
plot!(sol.t,pos[3,:])

#convert from ECI to ECEF for IGRF
t = [t0:dt:tf+dt*99...]
GMST = (280.4606 .+ 360.9856473*(t/24/60/60 .+ (MJD_0.+t./(24*60*60)) .- 51544.5))./180*pi #calculates Greenwich Mean Standard Time
ROT = zeros(3,3,length(t)) #cretes rotation matrix for every time step
pos_ecef = zeros(3,length(t)) #populates ecef position vector
lat = zeros(length(t),1)
long = zeros(length(t),1) # creates latitude and longitude
for i = 1:length(t)
    ROT[:,:,i] = [cos(GMST[i]) sin(GMST[i]) 0;
    -sin(GMST[i]) cos(GMST[i]) 0;
    0 0 1]
    pos_ecef[:,i] = ROT[:,:,i]*pos[:,i]
    lat[i] = pos_ecef[3,i]/norm(pos_ecef[:,i])
    long[i] = atan(pos_ecef[2,i]/pos_ecef[1,i])
end


#find dipole model of magnetic field
B0=-8E15/(1000)^3
B_N=zeros(size(pos,2),3)
for i=1:size(pos,2)
    B_N[i,:]=[3*pos[1,i]*pos[3,i]*(B0/norm(pos[:,i])^5);3*pos[2,i]*pos[3,i]*(B0/norm(pos[:,i])^5);(3*pos[3,i]^2-norm(pos[:,i])^2)*(B0/norm(pos[:,i])^5)];
end

#find the IGRF magnetic field for the orbit!
#magnetic field (IGRF)
B_N_sim = zeros(size(pos,2),3)
for i = 1:size(pos,2)
    B_N_sim[i,:] = SatelliteToolbox.igrf12(2019,norm(pos_ecef[:,i])*1000,lat[i],long[i])/1.e9
end
B_N_sim = convert(Array{Float64},B_N_sim)
plot(sol.t,B_N_sim[:,1])
plot!(sol.t,B_N_sim[:,2])
plot!(sol.t,B_N_sim[:,3])


#trajectory optimization section
#initial
omega_int=[0;0;0]
q_0=[0.;0;1;0]
x0=[omega_int;q_0;t0]

#final
q_final=[0.;1;0;0]
omega_final=[0.;0;0]
xf=[omega_final;q_final;1] #km or km/s

#enter into warm start
N_vect = N/10:N/10:N #create initial number of knot points
dt_vect=(tf-t0)./(N_vect)

#create model
n=8;
m=3;
model=Model(DerivFunction,n,m);

#LQR
Q=zeros(n,n);
Qf=zeros(n,n);
Q[1:3,1:3] = Array((1.e-1)*Diagonal(I,3));
Qf[1:3,1:3] = Array((1.e-3)*Diagonal(I,3));
α = 3.e3;
Q[4:7,4:7] = Array(α*Diagonal(I,4));
Qf[4:7,4:7] = Array(α*Diagonal(I,4))*1.e3;

R = Array((1.e2)*Diagonal(I,m));

#bounds
u_bnd=.005
x_bnd=10


#begin on warm start

for i = 1:length(N_vect)

    obj = TrajectoryOptimization.UnconstrainedObjective(Q, R, Qf, tf, x0, xf)
    obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
    solver = TrajectoryOptimization.Solver(model,obj,N=convert(Int64,N_vect[i]+1))
    solver.opts.verbose=true
    solver.opts.live_plotting=false
    solver.opts.iterations_outerloop=40
    solver.opts.sat_att=true
    solver.opts.γ=10

    if i == 1
        U = zeros(Float64,m,solver.N) #clears state and input
        X = zeros(Float64,m,solver.N)
        U = rand(Float64,m,solver.N)/1000 #initialize random input and state
        X = rand(Float64,m,solver.N)/1000
    end
    results, stats = TrajectoryOptimization.solve(solver,U)
    X = TrajectoryOptimization.to_array(results.X) #assigns state and input
    U = TrajectoryOptimization.to_array(results.U)
    K_lqr = TrajectoryOptimization.to_array(results.K)

    integration=:rk4
    t_warm = t0:dt_vect[i+1]:tf
    X, U = interpolate_trajectory( integration, X, U, t_warm) #splines over U and X

    print("Warm start has completed pass ",i, " of ", length(N_vect))
end

#plot input
plot(U[1,:])
plot!(U[2,:])
plot!(U[3,:])

#plot omega
plot(X[1,:])
plot!(X[2,:])
plot!(X[3,:])

#plot quaternion
plot(X[4,:])
plot!(X[5,:])
plot!(X[6,:])
plot!(X[7,:])

#plot time
plot(X[8,:])

#now test!

#assigns initial conditions
x0_lqr = zeros(n)
x0_lqr[1:3] = x0[1:3]

#creates noise over initial conditions
q_noise = randn(3,1)*(10*pi/180)^2
θ_noise = norm(q_noise)
r_noise = q_noise/θ_noise
x0_lqr[4:7] = qmult(x0[4:7],[cos(θ_noise/2);r_noise*sin(θ_noise/2)])

# #assigns new time vector
# interp_factor = 5
# N_lqr = N*interp_factor
# dt_lqr = dt/interp_factor
# t_lqr = [t0:dt_lqr:tf...]
#
# #define integration ds
# integration=:rk4
#
# #interpolate to find new state and input
# X_lqr, U_lqr = interpolate_trajectory( integration, X, U, t_lqr)


#since we are no longer interpolating, just reassign
X_lqr = X
U_lqr = U
dt_lqr = dt

#assigns simulation and gains functions for TVLQR
f! = simulator
f_gains! = gain_simulator

#creates Q and R matrices for TVLQR
α = 2.e-1
β = 1.e5
Q_lqr = zeros(6,6)
Q_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))
Q_lqr[4:6,4:6] = Array((α*β)*Diagonal(I,3))
Qf_lqr = zeros(6,6)
Qf_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))*1.e10
Qf_lqr[4:6,4:6] = Array((α*β)*Diagonal(I,3))*1.e5
R_lqr = Array((1.e1)*Diagonal(I,3))

#executes TVLQR
X_sim, U_sim, dX, K = attitude_simulation(f!,f_gains!,integration,X_lqr,U_lqr,dt_lqr,x0_lqr,tf,Q_lqr,R_lqr,Qf_lqr)

#plotting
plot(X_sim[1:3,:]')
plot(X_sim[4:7,:]')
# plot(X_sim[8,:])

plot(dX[1:3,:]')
plot(dX[4:6,:]')


plot(U_sim',color="red")
plot!(U_lqr',color="blue")

# plot(B_N,color="red")
# plot!(B_N_sim,color="blue")
