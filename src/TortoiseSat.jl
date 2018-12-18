#this code is for trajectory optimization of an underactuated satellite

using Pkg;
using LinearAlgebra;
Pkg.add("DifferentialEquations");
using DifferentialEquations;
using Plots;
Pkg.add("TrajectoryOptimization");
Pkg.add("SatProp")
using TrajectoryOptimization


#declare functions
include("legendre.jl")
include("dlegendre.jl")
include("kep_cart.jl")
include("cart_latlong.jl")
include("igrf.jl")
include("igrf12_coefs.jl")
include("igrf12syn_coefs.jl")
include("NED_cart.jl")
include("DerivFunction.jl")

#earth parameters
GM = 3.986004418E14*(1/1000)^3; #earth graviational constant (km^3 s^-2)
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
t_int=0.0;
t_fin=60.0*60*5; #seconds

#initialize rotation vector

#earth rotation
omega_E=360/(23*60*60+56*60+4); #degree rotation per second

# #initialize spice kernels
# addpath('/Users/agatherer/Documents/Old_Documents/MATLAB/mice/src/mice/')
# addpath('/Users/agatherer/Documents/Old_Documents/MATLAB/mice/lib/' )
# cspice_furnsh('naif0012.tls');
# cspice_furnsh('pck00010.tpc');
# cspice_furnsh('de430.bsp');
# et = cspice_str2et( '2019 JAN 31 01:00' ); %initial UTC time

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
pos_0=kep_cart(A,t_int,GM); #km
att_int=[.2; 0; 0]; #rads
omega_int=[.1; 0; 0]; #rads/s
(lat_0,long_0)=cart_latlong(pos_0[1,1:3]);
B_0=NED_cart(igrf12(2018,(alt+R_E)*1000,lat_0*pi/180,long_0*pi/180)*1E-9,lat_0,long_0);
u0=[pos_0[1,:];pos_0[2,:];att_int;omega_int;q_N_B]; #km or km/s

tspan=(t_int,t_fin);
p=(J,GM,mass,m_max);
prob=ODEProblem(DerivFunction!,u0,tspan,p);
sol=solve(prob,reltol=1e-5,abstol=1e-5);


#plot satellite path
plot3d(sol[1,:],sol[2,:],sol[3,:])

#plot omega
plot(sol.t,[sol[10,:] sol[11,:] sol[12,:]])

#create lat and long vectors
lat=zeros(size(sol,2),1);
long=zeros(size(sol,2),1);
for i=1:size(sol,2)
    (lat_temp,long_temp)=cart_latlong(sol[1:3,i]);
    lat[i]=lat_temp;
    long[i]=rem(long_temp-omega_E*i*(t_fin-t_int)/length(sol.t)+180,360)+180;
end

#plot earth track
plot(long,lat)
