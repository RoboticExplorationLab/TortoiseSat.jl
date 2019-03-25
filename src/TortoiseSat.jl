#this code is for trajectory optimization of an underactuated satellite

using Pkg;
using LinearAlgebra;
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
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

# #satellite parameters (1U)
# mass=.75 #kg
# J=[00.00125 0 0;
#     0 0.00125 0;
#     0 0 0.00125] #kgm2
# J_inv=inv(J)

#satellite parameters (3U)
mass=3.5 #kg
J=[0.005256 0 0;
    0 0.04939 0;
    0 0 0.04939] #kgm2
J_inv=inv(J)


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
t0=0.0 #seconds
tf=60.0*6+t0 #seconds

#initial
pos_0=kep_ECI(A,t0,GM) #km and km/s
u0=[pos_0[1,:];pos_0[2,:]] #km or km/s

# find magnetic field along 1 orbit
N=1000 #initialize the ultimate number of knot points
dt=(tf-t0)/(N)
tspan=(t0,tf*2) #go a little over for the ForwardDiff downstream
prob=DiffEqBase.ODEProblem(OrbitPlotter,u0,tspan)
sol=DiffEqBase.solve(prob,dt=dt,adaptive=false,RK4())
pos=sol[1:3,:]
vel=sol[4:6,:]
plot(sol.t,pos[1,:])
plot!(sol.t,pos[2,:])
plot!(sol.t,pos[3,:])

#convert from ECI to ECEF for IGRF
t = t0:dt:2*tf
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
plot(sol.t[1:end-1],B_N_sim[:,1])
plot!(sol.t[1:end-1],B_N_sim[:,2])
plot!(sol.t[1:end-1],B_N_sim[:,3])


#trajectory optimization section
#initial
omega_int=[0;0;0]
q_0=[0.;0;1;0]
x0=[omega_int;q_0;t0]

#final
q_final=[0.;1;0;0]
omega_final=[0.;0;0]
xf=[omega_final;q_final;1] #km or km/s

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
Q[1:3,1:3] = Array((1.e3)*Diagonal(I,3))
Qf[1:3,1:3] = Array((1.e1)*Diagonal(I,3))
α = 1e3
Q[4:7,4:7] = Array(α*Diagonal(I,4))
Qf[4:7,4:7] = Array(α*Diagonal(I,4))*1.e3

R = Array((10*norm(B_N_sim))*Diagonal(I,m))

#bounds
u_bnd=[2;2;2]
x_bnd=10

# ## now new cost function with minimum
# function quat_cost(x::AbstractVector, u::AbstractVector)
#     #this function uses the quaternion error rather than LQR for cost
#     0.5*(x[1:3] - xf[1:3])'*Q[1:3,1:3]*(x[1:3] - xf[1:3]) + 0.5*u'*R*u +
#     Q[5,5]*min(1-xf[4:7]'*x[4:7], 1+xf[4:7]'*x[4:7])
# end
#
# function quat_cost_terminal(x::AbstractVector)
#     #this function uses the quaternion error rather than LQR for cost, terminal constraint
#     0.5*(x[1:3] - xf[1:3])'*Qf[1:3,1:3]*(x[1:3] - xf[1:3]) +
#     Q[5,5]*min(1-xf[4:7]'*x[4:7], 1+xf[4:7]'*x[4:7])
# end
#
# function quat_cost_grad(x::AbstractVector,u::AbstractVector)
#     #this function calculates the gradient of the quaternion error cost funciton
#     if 1-xf[4:7]'*x[4:7] >= 1+xf[4:7]'*x[4:7]
#         q_grad = [(x[1:3] - xf[1:3])'*Q[1:3,1:3] Q[5,5]*xf[4:7]'*L_decomp(x[4:7]) 0]
#     else
#         q_grad = [(x[1:3] - xf[1:3])'*Q[1:3,1:3] -Q[5,5]*xf[4:7]'*L_decomp(x[4:7]) 0]
#     end
#     r_grad = u'*R
#     return q_grad,r_grad
# end
#
# function quat_cost_grad(xN::AbstractVector)
#     #this function calculates the gradient of the quaternion error cost funciton
#     if 1-xf[4:7]'*xN[4:7] >= 1+xf[4:7]'*xN[4:7]
#         qf_grad = [(xN[1:3] - xf[1:3])'*Q[1:3,1:3] Q[5,5]*xf[4:7]'*L_decomp(xN[4:7]) 0]
#     else
#         qf_grad = [(xN[1:3] - xf[1:3])'*Q[1:3,1:3] -Q[5,5]*xf[4:7]'*L_decomp(xN[4:7]) 0]
#     end
#     return qf_grad
# end
#
# function quat_cost_hess(x::AbstractVector,u::AbstractVector)
#     #this function calculates the hessian of the error quaternion cost function
#     H_hess = zeros(3,7)
#     Q_hess = [Q[1:3,1:3]' zeros(3,4);
#         zeros(4,7);
#         zeros(1,7)]
#     R_hess = R'
#     return Q_hess,R_hess,H_hess
# end
#
# function quat_cost_hess(xN::AbstractVector)
#     #this function calculates the hessian of the error quaternion cost function
#     H_hess = zeros(3,7)
#     Qf_hess = [Q[1:3,1:3]' zeros(3,4);
#         zeros(4,7);
#         zeros(1,7)]
#     return Qf_hess
# end

## now new cost function with minimum
function quat_cost(x::AbstractVector, u::AbstractVector)
    #this function uses the quaternion error rather than LQR for cost
    0.5*(x[1:3] - xf[1:3])'*Q[1:3,1:3]*(x[1:3] - xf[1:3]) + 0.5*u'*R*u +
    Q[5,5]*min(1-xf[4:7]'*x[4:7], 1+xf[4:7]'*x[4:7])
end

function quat_cost_terminal(x::AbstractVector)
    #this function uses the quaternion error rather than LQR for cost, terminal constraint
    0.5*(x[1:3] - xf[1:3])'*Qf[1:3,1:3]*(x[1:3] - xf[1:3]) +
    Q[5,5]*min(1-xf[4:7]'*x[4:7], 1+xf[4:7]'*x[4:7])
end

function quat_cost_grad(x::AbstractVector,u::AbstractVector)
    #this function calculates the gradient of the quaternion error cost funciton
    if 1-xf[4:7]'*x[4:7] >= 1+xf[4:7]'*x[4:7]
        q_grad = [(x[1:3] - xf[1:3])'*Q[1:3,1:3] Array(Q[5,5]*xf[4:7]') 0]
    else
        q_grad = [(x[1:3] - xf[1:3])'*Q[1:3,1:3] Array(-Q[5,5]*xf[4:7]') 0]
    end
    r_grad = u'*R
    return q_grad[1:end],r_grad[1:end]
end

function quat_cost_grad(xN::AbstractVector)
    #this function calculates the gradient of the quaternion error cost funciton
    if 1-xf[4:7]'*xN[4:7] >= 1+xf[4:7]'*xN[4:7]
        qf_grad = [(xN[1:3] - xf[1:3])'*Q[1:3,1:3] Array(Q[5,5]*xf[4:7]') 0]
    else
        qf_grad = [(xN[1:3] - xf[1:3])'*Q[1:3,1:3] Array(-Q[5,5]*xf[4:7]') 0]
    end
    return qf_grad[1:end]
end

function quat_cost_hess(x::AbstractVector,u::AbstractVector)
    #this function calculates the hessian of the error quaternion cost function
    H_hess = zeros(3,8)
    perm = [zeros(1,3);hat(xf[5:7])+Array(Diagonal(I,3))*xf[4]]'
    Q_hess = [Q[1:3,1:3]' zeros(3,5);
        zeros(4,8);
        zeros(1,8)]
    R_hess = R'
    return Q_hess,R_hess,H_hess
end

function quat_cost_hess(xN::AbstractVector)
    #this function calculates the hessian of the error quaternion cost function
    H_hess = zeros(8,3)
    perm = [zeros(1,3);hat(xf[5:7])+Array(Diagonal(I,3))*xf[4]]'
    Qf_hess = [Q[1:3,1:3]' zeros(3,5);
        zeros(4,8);
        zeros(1,8)]
    return Qf_hess
end


#random start to trajopt
# cost = GenericCost(quat_cost,quat_cost_terminal,quat_cost_grad,quat_cost_hess,n,m)
# obj = UnconstrainedObjective(cost,tf,x0,xf)
obj = LQRObjective(Q, R, Qf, tf, x0, xf)
obj_con=TrajectoryOptimization.ConstrainedObjective(obj,x_min=-x_bnd,x_max=x_bnd,u_min=-u_bnd,u_max=u_bnd)
solver = TrajectoryOptimization.Solver(model,obj,N=N)
solver.opts.verbose=true
solver.opts.iterations_outerloop=3
solver.opts.iterations_innerloop=100
solver.opts.dJ_counter_limit = 1
solver.opts.sat_att = true
U = rand(Float64,m,solver.N)/1000 #initialize random input and state
X = rand(Float64,m,solver.N)/1000
results, stats = TrajectoryOptimization.solve(solver,U)
X = TrajectoryOptimization.to_array(results.X) #assigns state and input
U = TrajectoryOptimization.to_array(results.U)


#plot input
plt = plot(t[1:N-1],U[1:3,:]'./100,label = ["U1" "U2" "U3"])
plot!(title = "3U CubeSat Slew", xlabel = "Time (s)", ylabel = "Magnetic Moment(Am2)")

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
# #interpolate to find new state and input
# X_lqr, U_lqr = interpolate_trajectory( integration, X, U, t_lqr)


#since we are no longer interpolating, just reassign
X_lqr = X
U_lqr = U
dt_lqr = dt

#assigns simulation and gains functions for TVLQR
include("simulator.jl")
f! = simulator
f_gains! = gain_simulator

#creates Q and R matrices for TVLQR
α = 1.e-1
β = 1.e3
Q_lqr = zeros(6,6)
Q_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))
Q_lqr[4:6,4:6] = Array((α*β)*Diagonal(I,3))
Qf_lqr = zeros(6,6)
Qf_lqr[1:3,1:3] = Array((α)*Diagonal(I,3))*1.e4
Qf_lqr[4:6,4:6] = Array((α*β)*Diagonal(I,3))*1.e8
R_lqr = Array((2.e-2)*Diagonal(I,3))

#executes TVLQR
X_sim, U_sim, dX, K = attitude_simulation(f!,f_gains!,integration,X_lqr,U_lqr,dt_lqr,x0_lqr,t0,tf,Q_lqr,R_lqr,Qf_lqr)

#plotting
plt = plot(t[1:N],X_sim[1:3,:]',label = ["omega1" "omega2" "omega3"])
plot!(title = "3U CubeSat Slew", xlabel = "Time (s)", ylabel = "Rotation Rate (rad/s)")

plt = plot(t[1:N],X_sim[4:7,:]',label = ["q1" "q2" "q3" "q4"])
plot!(title = "3U CubeSat Slew", xlabel = "Time (s)", ylabel = "Quaternion")

plt = plot(t[1:N],U_sim[1:3,:]'./100,label = ["u1" "u2" "u3"])
plot!(title = "3U CubeSat Slew", xlabel = "Time (s)", ylabel = "Magnetic Moment (Am2)")

plot(dX[1:3,:]')
plot(dX[4:6,:]')


plot(U_sim',color="red")
plot!(U_lqr',color="blue")

# plot(B_N,color="red")
# plot!(B_N_sim,color="blue")


function L_decomp(q)
    #this function takes in a desired quaternion computes the L decomposition
    L = [q[1] -q[2] -q[3] -q[4];
         q[2] q[1] -q[4] q[3];
         q[3] q[4] q[1] -q[2];
         q[4] -q[3] q[2] q[1]]

         return L

end
