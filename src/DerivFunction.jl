function DerivFunction(dx,x,u)
#Summary of this function goes here

omega=x[1:3]
q=x[4:7]./norm(x[4:7])
t=x[8]

#calculate rotation matrix from quaternion
# q_inv=[-q[1:3]; q[4]];
#
# v_hat=[0 q[3] -q[2];
#    -q[3] 0 q[1];
#    q[2] -q[1] 0];
#
# q=I+2*(v_hat)*(q[4]*I+v_hat);
skew_omega=[0 omega[3] -omega[2];
   -omega[3] 0 omega[1];
   omega[2] -omega[1] 0];
# q_dot=q*-skew_omega;

#find quaternion derivative
# q_dot_temp=[0*Array(Diagonal(I,3)*1.)+skew_omega omega;-omega' 0];
# q_dot=q_dot_temp*q;
q_dot=0.5*qmult(q,[0; omega])

#magnetic field
# B_N=[3*pos[1]*pos[3]*(B0/norm(pos[:])^5);3*pos[2]*pos[3]*(B0/norm(pos[:])^5);(3*pos[3]^2-norm(pos[:])^2)*(B0/norm(pos[:])^5)];
temp=qrot(q,B_N[floor(Int,t/tf*N+1),:])
# temp=qrot(q,B_N[1,:])
B_B=temp[1:3];

#m_c=zeros(3);
#for i=1:3
#   temp=cross(-omega,B_B1);
#   m_c[i]=-m_max*sign(temp[i]);
#end

tau_c=cross(u[1:3]/100,B_B); #N-m
# tau_c=u[1:3]

#rotation
omega_dot=J_inv*(tau_c-cross(omega,J*omega)); #rad/s^2

#redo concatination
dx[1:8] = [omega_dot;q_dot;1/(tf-t0)];
#print(omega);


end

function qrot(q,r)
      r + 2*cross(q[2:4],cross(q[2:4],r) + q[1]*r)
end

function qmult(q1,q2)
      [q1[1]*q2[1] - q1[2:4]'*q2[2:4]; q1[1]*q2[2:4] + q2[1]*q1[2:4] + cross(q1[2:4],q2[2:4])]
end
