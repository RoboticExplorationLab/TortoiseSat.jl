#define expansion function for TrajOpt
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
