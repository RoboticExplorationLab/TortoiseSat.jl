function NED_cart(vector,lat,long)
#this function converts a NED vector at a known location to a cartesian
#vector

north=vector[1];
east=vector[2];
down=vector[3];

Q=[-sind(lat)*cosd(long) -sind(long) -cosd(lat)*cosd(long);
    -sind(lat)*sind(long) cosd(long) -cosd(lat)*sind(long);
    cosd(lat) 0 -sind(lat)];

output=Q*vector;

end
