function magnetic_gramian(B_N,dt)
    #this function takes in a magnetic field of size
    # Nx3 and returns the magnetic field gramian B^T*B
    # of size 3x3xN
    B_gram = zeros(3,3,size(B_N,1))
    B_gram[:,:,1] = hat(B_N[1,:])*hat(B_N[1,:])'
    for i = 2:size(B_N,1)
        B_gram[:,:,i] = B_gram[:,:,i-1]+
        hat(B_N[i,:])*hat(B_N[i,:])'*dt
    end
    return B_gram
end

function condition_based_time(B_gram,cutoff)
    #this function takes in a 3x3xN array of magnetic field
    #gramian and spits out the time index at which the cutoff
    # is hit for the condition number

    B_gram_cond = zeros(1,size(B_gram,3))
    for i = 1:size(B_gram,3)
        (B_gram_cond[i],) = cond(B_gram[:,:,i])
    end #generate condition number of the gramian
    tf_index = 0
    for i = 1:size(B_gram_cond,2)
        if (B_gram_cond[i]) < cutoff
            tf_index = i
            break
        end
    end
    return tf_index
end

function magnetic_simulation(p,t0,tf,N,mag_field)
    #this funciton takes in an input and projects the magnetic field
    #and returns a magnetic field over the entire simulation

    # INPUTS:
    # paramters (p), initial time, final time, knot points, mag_field
    # OUTPUTS:
    # magnetic field in ECI, position, velocity

    A = p.Kep
    GM = p.GM
    R_E = p.R_E

    #initial state
    pos_0=kep_ECI(p.Kep,t0,p.GM) #km and km/s
    u0=[pos_0[1,:];pos_0[2,:]] #km or km/s

    # find magnetic field along 1 orbit
    tspan=(t0,tf*2) #go a little over for the ForwardDiff downstream
    prob=DiffEqBase.ODEProblem(OrbitPlotter,u0,tspan,p = p)
    dt = (tf-t0)/N
    sol=DiffEqBase.solve(prob, dt = dt ,adaptive=false,RK4())
    pos=sol[1:3,:]
    vel=sol[4:6,:]

    #convert from ECI to ECEF for IGRF
    t = t0:(tf-t0)/N:2*tf
    GMST = (280.4606 .+ 360.9856473*(t/24/60/60 .+ p.MJD) .- 51544.5)./180*pi #calculates Greenwich Mean Standard Time
    ROT = zeros(3,3,2*N) #creates rotation matrix for every time step
    pos_ecef = zeros(3,2*N) #populates ecef position vector
    lat = zeros(2*N,1)
    long = zeros(2*N,1) # creates latitude and longitude
    for i = 1:2*N
        ROT[:,:,i] = Rz(GMST[i])
        pos_ecef[:,i] = ROT[:,:,i]*pos[:,i]
        lat[i] = asin(pos_ecef[3,i]/norm(pos_ecef[:,i]))
        long[i] = atan(pos_ecef[2,i],pos_ecef[1,i])
    end

    #find the IGRF magnetic field for the orbit!
    B_N_sim = zeros(N*2,3)
    NED_to_ENU = [0 1 0;1 0 0;0 0 -1]

    for i = 1:size(B_N_sim,1)-1

        # using SatelliteToolbox
        # maximum polynomial of IGRF
        n_max = 3
        B_N_sim[i,:] = SatelliteToolbox.igrf12(2019,(alt+R_E)*1000,lat[i],long[i])/1.e9

        # # old way of interpolation
        # B_N_sim[i,1] = mag_field((lat[i]+π/2)/(π)*size(mag_field,1),(long[i]+π)/(2*π)*size(mag_field,1),1)
        # B_N_sim[i,2] = mag_field((lat[i]+π/2)/(π)*size(mag_field,1),(long[i]+π)/(2*π)*size(mag_field,1),2)
        # B_N_sim[i,3] = mag_field((lat[i]+π/2)/(π)*size(mag_field,1),(long[i]+π)/(2*π)*size(mag_field,1),3)

        # new way of calling MATLAB
        # B_N_sim[i,:] = mat"magnetic_estimation_400(0,0)"

        # convert to ECI
        R_ENU_to_XYZ =
            [-sin(long[i]) -sin(lat[i])*cos(long[i]) cos(lat[i])*cos(long[i]);
            cos(long[i]) -sin(lat[i])*sin(long[i]) cos(lat[i])*sin(long[i]);
            0 cos(lat[i]) sin(lat[i])];
        B_N_sim[i,:] = Rz(GMST[i])'*R_ENU_to_XYZ*NED_to_ENU*B_N_sim[i,:]

    end
    B_N_sim = convert(Array{Float64},B_N_sim)
    # plot(sol.t[1:end-1],B_N_sim[:,1])
    # plot!(sol.t[1:end-1],B_N_sim[:,2])
    # plot!(sol.t[1:end-1],B_N_sim[:,3])

    return B_N_sim, pos, vel

end

function igrf_data(altitude,year::Int64)
    #This function takes in an altitude in km and a year and returns an array
    #of 1E6 points evenly divided along the Earth's latitude and longitude
    R_E = 6378 #Earth radius in km
    N = 1000 # number of sample points
    mag_field = zeros(N,N,3)
    lat = collect(range(-π/2,length = N,π/2))
    long = collect(range(-π,length = N,π))
    for i = 1:size(mag_field,1)
        for j = 1:size(mag_field,2)
            mag_field[i,j,:] =
                SatelliteToolbox.igrf12(year,(altitude+R_E)*1000,lat[i],long[j])/1.e9
        end
    end

    #smooth the magnetic field because satellite toolbox has some discontinuities
    mag_field = Interpolations.interpolate(mag_field,BSpline(Cubic(Reflect(OnCell()))))
    mag_field = Interpolations.extrapolate(mag_field,Periodic())
    return mag_field
end


function Rx(theta)

    return [1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]

end

function Rz(theta)

    return [cos(theta) sin(theta) 0 ; -sin(theta) cos(theta) 0;0 0 1]

end
