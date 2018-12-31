function DerivFunction(dx,x,u)
#Summary of this function goes here
vel=x[1:3];
pos=x[4:6];
omega=x[7:9];
q_N_B=x[10:13];

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
B_N=[3*pos[1]*pos[3]*(B0/norm(pos[:])^5);3*pos[2]*pos[3]*(B0/norm(pos[:])^5);(3*pos[3]^2-norm(pos[:])^2)*(B0/norm(pos[:])^5)];
B_B1=Q_N_B*B_N;

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

#force balance
#gravity
#2-body gravity
f_grav=GM*mass/(norm(pos)^2)*-pos/(norm(pos)) #km*kg/s^2

#J2 oblateness term
J2=0.0010826359;
#J2=0;
f_J2=[J2*pos[1]/norm(pos)^7*(6*pos[3]-1.5*(pos[1]^2+pos[2]^2))
     J2*pos[2]/norm(pos)^7*(6*pos[3]-1.5*(pos[1]^2+pos[2]^2))
     J2*pos[3]/norm(pos)^7*(3*pos[3]-4.5*(pos[1]^2+pos[2]^2))];


#acceleration
a=(f_grav+f_J2)/mass; #km/s^2

#redo concatination
dx[1:13] = [a;vel;omega_dot;q_N_B_dot];
#print(omega);


end
