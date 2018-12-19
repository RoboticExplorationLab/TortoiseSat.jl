function DerivFunction(dx,x,u)
#Summary of this function goes here
omega=x[1:3];
q_N_B=x[4:7];
count=x[8]

#calculate rotation matrix from quaternion
q_N_B_inv=[-q_N_B[1:3]; q_N_B[4]];

v_hat=[0 q_N_B[3] -q_N_B[2];
   -q_N_B[3] 0 q_N_B[1];
   q_N_B[2] -q_N_B[1] 0];

Q_N_B=I+2*(v_hat)*(q_N_B[4]*I+v_hat);
skew_omega=[0 omega[3] -omega[2];
   -omega[3] 0 omega[1];
   omega[2] -omega[1] 0];
Q_N_B_dot=Q_N_B*-skew_omega;

#find quaternion derivative
q_N_B_dot=.5*q_N_B.*[omega;0];
q_N_B_dot_inv=[-q_N_B_dot[1:3];q_N_B_dot[4]];

#magnetic field
B_B1=Q_N_B*B_N[convert(Int64,N),:];

#m_c=zeros(3);
#for i=1:3
#   temp=cross(-omega,B_B1);
#   m_c[i]=-m_max*sign(temp[i]);
#end

# tau_c=cross(u[1:3],B_B1); #N-m
tau_c=u[1:3]

#rotation
M=tau_c;
omega_dot=inv(J)*(M-cross(omega,J*omega)); #rad/s^2

#redo concatination
dx[1:8] = [omega_dot;q_N_B_dot;1];
#print(omega);


end
