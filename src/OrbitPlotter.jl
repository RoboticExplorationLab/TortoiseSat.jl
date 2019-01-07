function OrbitPlotter(u,p,t)
#Summary of this function goes here
r = u[1:3] #km
v = u[4:6] #km/s

#force balance
#gravity
#2-body gravity
f_grav=GM*mass/(norm(r)^2)*-r/(norm(r)) #km*kg/s^2

#J2 oblateness term
J2=0.0010826359;
#J2=0;
f_J2=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))];


#acceleration
a=(f_grav+f_J2)/mass; #km/s^2

return [v;a];
#print(omega);


end
