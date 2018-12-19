#this code is for trajectory optimization of an underactuated satellite

using Pkg;
using LinearAlgebra;
using Plots;
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

#initialize vectors
r_N_1=[0 0 -1];
r_N_2=[0 1 0];
r_B_1=[0 0 -1];
r_B_2=[1 0 0];

#q-method
num_vect=2;
weights=ones(1,2);
normal_vects=[r_N_1;r_N_2]';
body_vects=[r_B_1;r_B_2]';
B = zeros(3,3);
z=zeros(3,1);
for i in 1:num_vect
    global B
    global z
    B=B+weights[i]*body_vects[:,i]*normal_vects[:,i]';
    z=z+weights[i]*cross(body_vects[:,i],normal_vects[:,i]);
end
K=[B+B'-tr(B)*I z;
     z' tr(B)];

eig_vals=eigvals(K);
eig_vects=eigvecs(K);
(eig__val_0,location)=findmax(eig_vals);
eig_vect_0=eig_vects[:,location];
q_N_B=eig_vect_0/norm(eig_vect_0);

#initial
pos_0=kep_cart(A,t0,GM); #km and km/s
att_int=[.2; .5; 0]; #rads
omega_int=[.02; 0; 0]; #rads/s
#find initial magnetic field
B0=-8E15*(1E-3)^3;
u0=[pos_0[1,:];pos_0[2,:];att_int;omega_int;q_N_B]; #km or km/s

#find magnetic field along 1 orbit
tspan=(t0,tf);
prob=DiffEqBase.ODEProblem(OrbitPlotter,u0,tspan);
sol=DiffEqBase.solve(prob,dt=1.0,Euler())
pos=sol[1:3,:]
# plot(sol.t,pos[1,:])
# plot!(sol.t,pos[2,:])
# plot!(sol.t,pos[3,:])

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
omega_int=[.2;0;0];
q_N_B=[0;0;0;0];
x0=[omega_int;q_N_B;1];

#final
q_N_B_final=[0;0;0;0]
omega_final=[0;0;0];
xf=[omega_final;q_N_B_final;size(B_N,1)]; #km or km/s

#bounds
u_bnd=1;
x_bnd=1;

#create model
n=8;
m=3;
model=Model(DerivFunction,n,m);

#LQR
Q = Array((1.e-9)*Diagonal(I,n));
Qf = Array(1.e-9*Diagonal(I,n));
R = Array((1.e-9)*Diagonal(I,m));

#constant stage cost
c=1.0

#determine number of nodes
N=size(B_N,1);

obj = TrajectoryOptimization.UnconstrainedObjective(Q, R, Qf, c, tf, x0, xf);
obj_con = ConstrainedObjective(obj,u_min=-u_bnd, u_max=u_bnd, x_min=-x_bnd, x_max=x_bnd);
solver = TrajectoryOptimization.Solver(model,obj,N=N);
solver.opts.verbose=true;
solver.opts.use_static = false
solver.opts.min_dt = dt
solver.opts.max_dt = dt
# solver.opts.iterations=5
# solver.opts.iterations_outerloop=5
# solver.opts.cost_intermediate_tolerance=1.0e-2;
# solver.opts.cost_tolerance=1.0e-2;
# solver.opts.gradient_intermediate_tolerance=1.0e-2
# solver.opts.gradient_tolerance=1.0e-2
# solver.opts.live_plotting=true


U = zeros(m,solver.N);


results, stats = TrajectoryOptimization.solve(solver,U)
X = TrajectoryOptimization.to_array(results.X)
plot(X[1,:])
plot!(X[2,:])
plot!(X[3,:])

plot(X[8,:])
