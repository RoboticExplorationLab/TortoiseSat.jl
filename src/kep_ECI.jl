function kep_ECI(kep_elements,t0,GM)
#this function converts keplerian orbital elements and a mean anomaly to
#cartesian coordinates

    #assumes all keplerian coordinates are in degrees

    A_temp=kep_elements
    A_temp[6]=rem(kep_elements[6]+t0*sqrt(GM./kep_elements[2].^3),360)

    #Use Newton Rhapson to find Eccentric Anomaly
    E=zeros(1,101);
    E[1]=A_temp[6]/180*pi;
    for i=1:100
        E[i+1]=E[i]-(E[i]-A_temp[1].*sin(E[i])-A_temp[6]/180*pi)./
            (1-A_temp[1].*cos(E[i]));
    end

        #find true anomaly
        nu=2*atand(sqrt(1+A_temp[1]).*sin(E[end]/2),
        sqrt(1-A_temp[1]).*cos(E[end]/2))

        #find distance to central body
        r_c=A_temp[2].*(1-A_temp[1].*cos(E[end]))

        #find orbital position and velocity
        o=r_c*[cosd(nu) sind(nu) 0]
        o_dot=sqrt(GM*A_temp[2])/r_c*[-sin(E[end]) sqrt(1-A_temp[1]^2)*cos(E[end]) 0];

        #transform into cartesian
        r = R_z(-A_temp[4])*R_x(-A_temp[3])*R_z(-A_temp[5])*o'
        rdot = R_z(-A_temp[4])*R_x(-A_temp[3])*R_z(-A_temp[5])*o_dot'

    output=[r';rdot']

end

function R_z(angle)

    return [cosd(angle) sind(angle) 0;
            -sind(angle) cosd(angle) 0;
            0 0 1]
end

function R_x(angle)

    return [1 0 0;
            0 cosd(angle) sind(angle);
            0 -sind(angle) cosd(angle)]
end
