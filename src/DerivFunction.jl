function DerivFunction!(du,u,p,t)
#Summary of this function goes here
r = u[1:3]; #km
v = u[4:6]; #km/s
att=u[7:9];
omega=u[10:12];
q_N_B=u[13:16];

J,GM,mass,m_max = p;

#force balance
#gravity
#2-body gravity
f_grav=GM*mass/(norm(r)^2)*-r/(norm(r)); #km*kg/s^2

#J2 oblateness term
J2=0.0010826359;
#J2=0;
f_J2=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))];


#acceleration
a=(f_grav+f_J2)/mass; #km/s^2

    #torques

        #rotation matrices
#         r_N_1=r/norm(r);
#         r_N_2=cross(r,v)/norm(cross(r,v));
#         pointing=[cos(att(2))*cos(att(1)) sin(att(2))*cos(att(1)) sin(att(1))];
#         r_B_1=[cos(att(2))*cos(att(1)) sin(att(2))*cos(att(1)) sin(att(1))];
#         r_B_2=cross(pointing,v)/norm(cross(pointing,v));


#magnetic field
dt=.01; #s
r2=r-v*dt;
r3=r-2*v*dt;
r4=r-3*v*dt;
(lat1,long1)=cart_latlong(r);
(lat2,long2)=cart_latlong(r2);
(lat3,long3)=cart_latlong(r3);
(lat4,long4)=cart_latlong(r4);

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

#magnetic fields
#B_B1=Q_N_B*NED_cart(igrf12(2018,norm(r)*1000,lat1*pi/180,long1*pi/180)*1E-9,lat1,long1);
#B_B2=Q_N_B*NED_cart(igrf12(2018,norm(r2)*1000,lat2*pi/180,long2*pi/180)*1E-9,lat2,long2);
#B_B_dot=(B_B2-B_B1)/dt; #T/s
#use Vandemonde Matrix to find cubic regression
B0=-8E15*(1E-3)^3;
B_B1=Q_N_B*[3*r[1]*r[3]*(B0/norm(r)^5);3*r[2]*r[3]*(B0/norm(r)^5);(3*r[3]^2-norm(r)^2)*(B0/norm(r)^5)];
B_B2=Q_N_B*[3*r2[1]*r2[3]*(B0/norm(r2)^5);3*r2[2]*r2[3]*(B0/norm(r2)^5);(3*r2[3]^2-norm(r2)^2)*(B0/norm(r2)^5)];
B_B3=Q_N_B*[3*r3[1]*r3[3]*(B0/norm(r3)^5);3*r3[2]*r3[3]*(B0/norm(r3)^5);(3*r3[3]^2-norm(r3)^2)*(B0/norm(r3)^5)];
B_B4=Q_N_B*[3*r4[1]*r4[3]*(B0/norm(r4)^5);3*r4[2]*r4[3]*(B0/norm(r4)^5);(3*r4[3]^2-norm(r4)^2)*(B0/norm(r4)^5)];
B_B_dot=-11/6*B_B1+3*B_B2-3/2*B_B3+1/3*B_B4; #T/s

m_c=zeros(3);
for i=1:3
   temp=cross(-omega,B_B1);
   m_c[i]=-m_max*sign(temp[i]);
end
tau_c=cross(m_c,B_B1); #N-m

#rotation
M=tau_c;
omega_dot=inv(J)*(M-cross(omega,J*omega)); #rad/s^2

du[:] = [v;a;omega;omega_dot;q_N_B_dot];
#print(omega);


end
