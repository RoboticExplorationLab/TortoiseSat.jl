function OrbitPlotter(x,p,t)
#this function plots an orbit of a satellite to get the position over a time period
r = x[1:3] #km
v = x[4:6] #km/s

#force balance
#gravity
#2-body gravity
f_grav=GM/(norm(r)^2)*-r/(norm(r)) #km/s^2

# atmospheric drag (assuming that time is invariant)
# date is new years on 2019
altitude=norm(r)-R_E #meters

#find the latitude and longitude in a geometric frame
GMST = (280.4606 + 360.9856473*(t/24/60/60 .+ (MJD_0.+t./(24*60*60)) .- 51544.5))./180*pi #calculates Greenwich Mean Standard Time
ROT = [cos(GMST) sin(GMST) 0;-sin(GMST) cos(GMST) 0; 0 0 1]
pos_ecef = ROT*r
lat= pos_ecef[3]/norm(pos_ecef)
long= atan(pos_ecef[2]/pos_ecef[1])

# out=SatelliteToolbox.nrlmsise00(DatetoJD(2019, 1, 1, 00, 00, 00), # Julian Day
    #              altitude,# Altitude [m]
    #              -lat*pi/180,# Latitude [rad]
    #              -long*pi/180,# Longitude [rad]
    #              83.7,# 81 day average F10.7 flux
    #              102.5,# Daily F10.7 for previous day
    #              15;# Magnetic index (daily)
    #              output_si = false# Output in cm^-3 and g/cm^-3
    # );
# ω_earth = [0;0;7.2921150E-5] #rad/s
# v_ecef = v-cross(ω_earth,r)
# f_drag=.5*out.den_Total*norm(v_ecef)^2*mass/BC*-v_ecef/norm(v_ecef)

#J2 oblateness term
J2=0.0010826359;
#J2=0;
f_J2=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))];


#acceleration
a=(f_grav+f_J2); #km/s^2

return [v;a]
#print(omega);


end
