#define expansion function for TrajOpt

function normal_expansion(cost::QuadraticCost, x::AbstractVector, u::AbstractVector)

    return cost.Q, cost.R, cost.H, (cost.Q*x+cost.q), cost.R*u+ cost.r

end

function normal_expansion(cost::QuadraticCost, xN::AbstractVector{Float64})

    return cost.Qf, (cost.Qf*xN+cost.qf)
end


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

function normal_error(X1,X2)

    return X1 - X2

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

# now new cost function with minimum
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

function quaternion_euler(q)
    #this function converts a quaternion into its euler angles

    euler = [atan(2*(q[1]*q[2] + q[3]*q[4]),1-2*(q[2]^2+q[3]^2));
            asin(2*(q[1]*q[3] - q[4]*q[2]));
            atan(2*(q[1]*q[4] + q[2]*q[3]),1-2*(q[3]^2+q[4]^2))]

    return euler

end

function euler_quaternion(euler)
    #this function (poorly) converts a set of euler angles into a quaternion

    ϕ = euler[1]
    θ = euler[2]
    ψ = euler[3]

    quat = [cos(ϕ/2)*cos(θ/2)*cos(ψ/2) +
    sin(ϕ/2)*sin(θ/2)*sin(ψ/2);
    sin(ϕ/2)*cos(θ/2)*cos(ψ/2) -
    cos(ϕ/2)*sin(θ/2)*sin(ψ/2);
    cos(ϕ/2)*sin(θ/2)*cos(ψ/2) +
    sin(ϕ/2)*cos(θ/2)*sin(ψ/2);
    cos(ϕ/2)*cos(θ/2)*sin(ψ/2) -
    sin(ϕ/2)*sin(θ/2)*cos(ψ/2)]

    return quat

end

function quaternion_rot(q)
    # this function takes in a quaternion and spits out a rotation matrix
    w = q[1]
    x = q[2]
    y = q[3]
    z = q[4]
    ROT = [1-2y^2-2z^2 2x*y-2z*w 2x*z+2y*w;
       2x*y+2z*w 1-2x^2-2z^2 2y*z-2x*w;
       2x*z-2y*w 2y*z+2x*w 1-2x^2-2y^2]

end

function rot_quaternion(ROT)
    # this function takes in a rotation matrix and spits out a quaternion

    # w= sqrt(1 + ROT[1,1] + ROT[2,2] + ROT[3,3])/2
    # x = (ROT[3,2] - ROT[2,3])/(4*w)
    # y = (ROT[1,3] - ROT[3,1])/(4*w)
    # z = (ROT[2,1] - ROT[1,2])/(4*w)
    #
    # q = [w;x;y;z]

    r11 = ROT[1,1]
    r12 = ROT[1,2]
    r13 = ROT[1,3]
    r21 = ROT[2,1]
    r22 = ROT[2,2]
    r23 = ROT[2,3]
    r31 = ROT[3,1]
    r32 = ROT[3,2]
    r33 = ROT[3,3]


    # shepperd's method

    i = sortperm([(r11+r22+r33);r11;r22;r33])[end]

    if i == 1
        r1 = sqrt(1+r11+r22+r33)
        q = .5*[r1;(r32-r23)/(r1);(r13-r31)/r1;(r21-r12)/r1]
    elseif i == 2
        r2 = sqrt(1+r11-r22-r33)
        q = .5*[(r32-r23)/(r2);r2;(r12+r21)/r2;(r31+r13)/r2]
    elseif i == 3
        r3 = sqrt(1-r11+r22-r33)
        q = .5*[(r13-r31)/r3;(r12+r21)/r3;r3;(r23+r32)/r3]
    elseif i == 4
        r4 = sqrt(1-r11-r22+r33)
        q = .5*[(r21-r12)/r4;(r31+r13)/r4;(r32+r23)/r4;r4]
    else
        println("error in Shepperd's method")
    end

    return q

end
