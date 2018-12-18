function cart_latlong(r)
#this function converts between cartesian coordinates and a latitude and
#longitude in degrees for the Earth

lat=asind(r[3]./norm(r));
long=atan(r[2],r[1])*180/pi;

lat,long

end
