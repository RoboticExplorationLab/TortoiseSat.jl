#this code is for trajectory optimization of an underactuated satellite

using Pkg;
using LinearAlgebra;
using Plots;
Pkg.add("Formatting")
Pkg.add("TrajectoryOptimization");
using TrajectoryOptimization;
Pkg.add("DifferentialEquations");
using DifferentialEquations;


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
J=[1E-4 0 0;
    0 1E-4 0;
    0 0 1E-4];

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

#find magnetic field along 1 orbit
tspan=(t0,tf);
prob=DiffEqBase.ODEProblem(OrbitPlotter,u0,tspan);
sol=DiffEqBase.solve(prob,dt=1.0,Euler())
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
plot(sol.t,B_N[:,1])
plot!(sol.t,B_N[:,2])
plot!(sol.t,B_N[:,3])
dt=(size(B_N,1)-1)/tf;

#trajectory optimization section
#initial
omega_int=[0.001;0;0]
q_N_B=[0.;0;1;0]
x0=[vel[1:3,1]/R_E;pos[1:3,1]/R_E;omega_int;q_N_B];

#final
q_N_B_final=[0.;1;0;0]
omega_final=[0.;0;0]
xf=[vel[1:3,end]/R_E;pos[1:3,end]/R_E;omega_final;q_N_B_final]; #km or km/s

#create model
n=13;
m=3;
model=Model(DerivFunction,n,m);

#LQR
Q=zeros(n,n);
Q[1:9,1:9] = Array((1.e-9)*Diagonal(I,9));
Q[10:13,10:13] = Array((1.e-5)*Diagonal(I,4));
Qf = Array((1.)*Diagonal(I,n));
R = Array((1.e-5)*Diagonal(I,m));

#bounds
u_bnd=.0001
x_bnd=100

obj = TrajectoryOptimization.UnconstrainedObjective(Q, R, Qf, tf, x0, xf);
obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
solver = TrajectoryOptimization.Solver(model,obj_con,N=100);
solver.opts.verbose=true;
# solver.opts.use_static = false
# solver.opts.min_dt = dt
# solver.opts.max_dt = dt
# solver.opts.cost_intermediate_tolerance=1.0e-3
# solver.opts.cost_tolerance=1.0e-4
# solver.opts.gradient_intermediate_tolerance=1.0e-5
# solver.opts.gradient_tolerance=1.0e-5
solver.opts.live_plotting=false
solver.opts.iterations_outerloop=5
solver.opts.sat_att=true


U = zeros(m,solver.N);


results, stats = TrajectoryOptimization.solve(solver,U)
X = TrajectoryOptimization.to_array(results.X)
U = TrajectoryOptimization.to_array(results.U)

#plot input
plot(U[1,:])
plot!(U[2,:])
plot!(U[3,:])

#plot velocity
plot(X[1,:]*R_E)
plot!(X[2,:]*R_E)
plot!(X[3,:]*R_E)

#plot position
plot(X[4,:]*R_E)
plot!(X[5,:]*R_E)
plot!(X[6,:]*R_E)

#plot omega
plot(X[7,:])
plot!(X[8,:])
plot!(X[9,:])

#plot quaternion
plot(X[10,:])
plot!(X[11,:])
plot!(X[12,:])
plot!(X[13,:])
