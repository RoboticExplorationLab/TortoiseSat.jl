function kep_cart(A,t,GM)
#this function converts keplerian orbital elements and a mean anomaly to
#cartesian coordinates

    A_temp=A;
    A_temp[6]=rem(A[6]+t*sqrt(GM./A[2].^3)*180/pi,360);

    #Use Newton Rhapson to find Eccentric Anomaly
    E=zeros(1,11);
    E[1]=A_temp[6];
    for i=1:10
        E[i+1]=E[i]-(E[i]-A_temp[1].*sind(E[i])-A_temp[6])./
            (1-A_temp[1].*cosd(E[i]));
    end

        #find true anomaly
        nu=2*atan(sqrt(1+A_temp[1]).*sind(E[end]/2),
        sqrt(1-A_temp[1]).*cosd(E[end]/2));

        #find distance to central body
        r_c=A_temp[2].*(1-A_temp[1].*cos(E[end]));

        #find orbital position and velocity
        o=r_c*[cos(nu) sin(nu) 0];
        o_dot=sqrt(GM*A_temp[2])/r_c*[-sin(E[end]) sqrt(1-A_temp[1]^2)*cos(E[end]) 0];

        #transform into cartesian
        r=[o[1].*(cosd(A[5]).*cosd(A[4])-sind(A[5]).*cosd(A[3]).*sind(A[4])) -
    o[2].*(sind(A[5]).*cosd(A[4])-cosd(A[5]).*cosd(A[3]).*sind(A[4]))   #r1
    o[1].*(cosd(A[5]).*sind(A[4])+sind(A[5]).*cosd(A[3]).*sind(A[4])) +
    o[2].*(cosd(A[5]).*cosd(A[3].*cosd(A[4])-sind(A[5]).*sind(A[4])))  #r2
    o[1].*(sind(A[5]).*sind(A[3]))+o[2].*(cosd(A[5]).*sind(A[3]))]';#r3
        rdot=[
    o_dot[1].*(cosd(A[5]).*cosd(A[4])-sind(A[5]).*cosd(A[3]).*sind(A[4])) -
    o_dot[2].*(sind(A[5]).*cosd(A[4])-cosd(A[5]).*cosd(A[3]).*sind(A[4]))  #rdot1
    o_dot[1].*(cosd(A[5]).*sind(A[4])+sind(A[5]).*cosd(A[3]).*sind(A[4])) +
    o_dot[2].*(cosd(A[5]).*cosd(A[3].*cosd(A[4])-sind(A[5]).*sind(A[4])))  #rdot2
    o_dot[1].*(sind(A[5]).*sind(A[3]))+o_dot[2].*(cosd(A[5]).*sind(A[3]))]'; #rdot3

    output=[r;rdot];

end
