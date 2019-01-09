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


#declare functions
include("kep_cart.jl")
include("cart_latlong.jl")
include("NED_cart.jl")
include("DerivFunction.jl")
include("OrbitPlotter.jl")

#earth parameters
GM = 3.986004418E14*(1/1000)^3; #earth graviational constant (km^3 s^-2)""
R_E = 6371.0; #earth radius (km)

#satellite parameters
mass=4.0; #kg
J=[1E-3 0 0;
    0 1E-3 0;
    0 0 1E-3];
J_inv=inv(J)

#control parameters
m_max=1E-4;

#enter orbital parameters
num_sats=1;
alt=600; #altitude in km
A=zeros(num_sats,6);
A[1] = 0; #circular orbits
A[2]= alt+R_E; #semi-major axis
A[3] = 45; #inclination in degrees
A[4] = 0; #assign ascending nodes
A[5] = 0; #argument of periapsis
A[6] = 0; #initial true anomoly
T=2*pi*sqrt(A[2].^3/GM);

#time initialization
t0=0.0;
tf=60.0*60*1.6; #seconds
dt=0.1; #seconds

#initial
pos_0=kep_cart(A,t0,GM); #km and km/s
u0=[pos_0[1,:];pos_0[2,:]]; #km or km/s

# find magnetic field along 1 orbit
N=100;
tspan=(t0,tf);
dt=(tf-t0)/(N-1.001)
prob=DiffEqBase.ODEProblem(OrbitPlotter,u0,tspan);
sol=DiffEqBase.solve(prob,dt=dt,Euler())
pos=sol[1:3,:]
vel=sol[4:6,:]
# plot(sol.t,pos[1,:])
# plot!(sol.t,pos[2,:])
# plot!(sol.t,pos[3,:])

B0=-8E15*(1E-3)^3;
B_N=zeros(size(pos,2),3)
for i=1:size(pos,2)
    B_N[i,:]=[3*pos[1,i]*pos[3,i]*(B0/norm(pos[:,i])^5);3*pos[2,i]*pos[3,i]*(B0/norm(pos[:,i])^5);(3*pos[3,i]^2-norm(pos[:,i])^2)*(B0/norm(pos[:,i])^5)];
end
# plot(sol.t,B_N[:,1])
# plot!(sol.t,B_N[:,2])
# plot!(sol.t,B_N[:,3])
# dt=(size(B_N,1)-1)/tf;

#trajectory optimization section
#initial
omega_int=[0.002;0;0]
q_0=[0.;0;1;0]
x0=[omega_int;q_0;t0];

#final
q_final=[0.;1;0;0]
omega_final=[0.;0;0]
xf=[omega_final;q_final;1]; #km or km/s

#create model
n=8;
m=3;
model=Model(DerivFunction,n,m);

#LQR
Q=zeros(n,n);
Qf=zeros(n,n);
Q[1:3,1:3] = Array((1.e-1)*Diagonal(I,3));
Qf[1:3,1:3] = Array((1.e-3)*Diagonal(I,3));
α = 1.e2;
Q[4:7,4:7] = Array(α*Diagonal(I,4));
Qf[4:7,4:7] = Array(α*Diagonal(I,4))*1.e2;

R = Array((1.e6)*Diagonal(I,m));

#bounds
u_bnd=.005
x_bnd=10

obj = TrajectoryOptimization.UnconstrainedObjective(Q, R, Qf, tf, x0, xf);
obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
solver = TrajectoryOptimization.Solver(model,obj_con,N=N)
solver.opts.verbose=true;
# solver.opts.use_static = false
# solver.opts.min_dt = dt
# solver.opts.max_dt = dt
solver.opts.live_plotting=false
solver.opts.iterations_outerloop=40
solver.opts.sat_att=true
solver.opts.γ=10

# U = rand(Float64,m,solver.N)/1000;
# U = zeros(m,solver.N);


results, stats = TrajectoryOptimization.solve(solver,U)
X = TrajectoryOptimization.to_array(results.X)
U = TrajectoryOptimization.to_array(results.U)

# dircol
# dt=1.
# results, stats = TrajectoryOptimization.solve_dircol(solver,X,U)

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
