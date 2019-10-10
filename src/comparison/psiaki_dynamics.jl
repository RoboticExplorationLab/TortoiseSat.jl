function psiaki_controller(C_1, C_2, J, q, ω̄, B_meas, m_limit)
    #this function returns the magnetic momment based on a PD controller
    #this is in the body frame
    # B̃ is the linearized B matrix
    # C_1 and C_2 are the weights
    # B_meas and B_exp are the measured magnetic field and the
    # expected magnetic field

    #required torque
    T_req = -(C_1 * ω̄ + C_2 * inv(J) * q[2:4])
    # T_req = -(C_1 * ω̄ + C_2 * inv(J) * quaternion_euler(q))

    #feasible angular momentum
    m̄ = (cross(B_meas, T_req))/(norm(B_meas)^2)

    #resize depending on limits
    # for i in 1:3
    #     if abs(m̄[i]) > m_limit && m̄[i] > 0
    #         m̄[i] = m_limit
    #     elseif abs(m̄[i]) > m_limit && m̄[i] < 0
    #         m̄[i] = -m_limit
    #     end
    # end

    return m̄
end

function psiaki_dynamics(A,B,x,m)
    #state vector x:
    # [ϕ_dot θ_dot ψ_dot ϕ θ ψ]'
    x_dot = A*x+B*m
end

function psiaki_linearization(J,ω_0,B)
    #constructs A and B matrix
    #ω_0 is the orbital mean motion
    #J is the moment of inertia matrix
    #B is the magnetic field

    σ = [(J[2,2]-J[3,3])/J[1,1],
    (J[3,3]-J[1,1])/J[2,2],
    (J[1,1]-J[2,2])/J[3,3]]

    A = zeros(6,6)
    A[1,3] = ω_0 * (1-σ[1])
    A[1,4] = -4 * ω_0^2 * σ[1]
    A[2,5] = 3 * ω_0^2 * σ[2]
    A[3,1] = -ω_0 * (1+σ[3])
    A[3,6] = ω_0^2 * σ[3]
    A[4:6,1:3] = Matrix{Float64}(I,3,3)

    B̃ = zeros(6,3)
    B̃[1,2] = B[3]/J[1,1]
    B̃[1,3] = -B[2]/J[1,1]
    B̃[2,1] = -B[3]/J[2,2]
    B̃[2,3] = B[1]/J[2,2]
    B̃[3,1] = B[2]/J[3,3]
    B̃[3,2] = -B[1]/J[3,3]

    return A, B̃
end

function rk4_psiaki(f,x,dt,u,B_B,J)

    u = u
    B_B = B_B
    ḟ_1 = f(x,u,B_B,J)
    ḟ_2 = f(x + .5*ḟ_1*dt,u,B_B,J)
    ḟ_3 = f(x + .5*ḟ_2*dt,u,B_B,J)
    ḟ_4 = f(x + ḟ_3*dt,u,B_B,J)

    x_next = x + 1/6 * (ḟ_1 + 2*ḟ_2 + 2*ḟ_3 + ḟ_4) * dt
end

function magnetic_control_effectiveness(B,J,R,t)
    # this function takes in a R weight and the
    # magnetic field and finds the control
    # effectivenss matrix
    # also uses a time vector t

    T = t[end] - t[1]
    dt = t[end] - t[end-1]

    α = zeros(6,6)
    R_inv = inv(R)
    for i in 1:size(B,2)
        B_matrix = zeros(6,3)
        B_matrix[4:6,1:3] = -hat(B[:,i])
        B_matrix[4,:] *= 1/J[1,1]
        B_matrix[5,:] *= 1/J[2,2]
        B_matrix[6,:] *= 1/J[3,3]

        α += B_matrix * R_inv * B_matrix' * dt
    end
    α = α/T

    return sqrt(α)
end


function Pss_psiaki(ω_0,J,Q,B̃)
    # this function calcualtes the steady state
    # Pss matrix using linearized dynamics
    # ω_0 is the orbital mean motion
    # J is the MOI matrix
    # Q is a weight
    # B̃ is the control effectiveness matrix

    σ = [(J[2,2]-J[3,3])/J[1,1];
        (J[3,3]-J[1,1])/J[2,2];
        (J[1,1]-J[2,2])/J[3,3]]

    A = zeros(6,6)
    A[1:3,4:6] = [1 0 0;0 1 0;0 0 1]
    A[4,1] = -4*ω_0^2*σ[1]
    A[4,6] = ω_0*(1-σ[1])
    A[5,2] = 3*ω_0^2*σ[2]
    A[6,3] = ω_0^2*σ[3]
    A[6,4] = -ω_0*(1+σ[3])

    Pss = ControlSystems.care(A, B̃, Q, Array(1.0*Diagonal(I,6)))

    return Pss

end

function nominal_input_psiaki(R_inv,B_matrix,Pss,x_linear,α,U_max)
    # this function find the nominal input given a contstraint

    u_nom = -α*R_inv*B_matrix'*Pss*x_linear

    β = maximum(abs.(u_nom)./U_max)

    u = zeros(3)
    if β <= 1
        u = u_nom
    elseif β > 1
        u = 1/β * u_nom
    end

    return u

end
