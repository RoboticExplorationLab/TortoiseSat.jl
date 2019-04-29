function magnetic_gramian(B_N)
    #this function takes in a magnetic field of size
    # Nx3 and returns the magnetic field gramian B^T*B
    # of size 3x3xN
    B_gram = zeros(3,3,size(B_N,1))
    B_gram[:,:,1] = hat(B_N[1,:])*hat(B_N[1,:])'
    for i = 2:size(B_N,1)
        B_gram[:,:,i] = B_gram[:,:,i-1]+
        hat(B_N[i,:])*hat(B_N[i,:])'
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
    end #generate eigenvalues of the gramian
    tf_index = 0
    for i = 1:size(B_gram_cond,2)
        if (B_gram_cond[i]) < cutoff
            tf_index = i
            break
        end
    end
    return tf_index
end

function magnetic_simulation(A,t0,tf,N)
    #this funciton takes in an input and projects the magnetic field
    #and returns a magnetic field over the entire simulation

    #initial state
    pos_0=kep_ECI(A,t0,GM) #km and km/s
    u0=[pos_0[1,:];pos_0[2,:]] #km or km/s

    # find magnetic field along 1 orbit
    dt=(tf-t0)/(N)
    tspan=(t0,tf*2) #go a little over for the ForwardDiff downstream
    prob=DiffEqBase.ODEProblem(OrbitPlotter,u0,tspan)
    sol=DiffEqBase.solve(prob,dt=dt,adaptive=false,RK4())
    pos=sol[1:3,:]
    vel=sol[4:6,:]
    # plot(sol.t,pos[1,:])
    # plot!(sol.t,pos[2,:])
    # plot!(sol.t,pos[3,:])

    #convert from ECI to ECEF for IGRF
    t = t0:dt:2*tf
    GMST = (280.4606 .+ 360.9856473*(t/24/60/60 .+ (MJD_0.+t./(24*60*60)) .- 51544.5))./180*pi #calculates Greenwich Mean Standard Time
    ROT = zeros(3,3,length(t)) #cretes rotation matrix for every time step
    pos_ecef = zeros(3,length(t)) #populates ecef position vector
    lat = zeros(length(t),1)
    long = zeros(length(t),1) # creates latitude and longitude
    for i = 1:length(t)-1
        ROT[:,:,i] = [cos(GMST[i]) sin(GMST[i]) 0;
        -sin(GMST[i]) cos(GMST[i]) 0;
        0 0 1]
        pos_ecef[:,i] = ROT[:,:,i]*pos[:,i]
        lat[i] = pos_ecef[3,i]/norm(pos_ecef[:,i])
        long[i] = atan(pos_ecef[2,i],pos_ecef[1,i])
    end

    #interpolate the magnetic field to remove any unsmoothe issues found
    #in satellite toolbox package
    mag_field = igrf_data(alt,2019)

    #find the IGRF magnetic field for the orbit!
    B_N_sim = zeros(size(pos,2)-1,3)
    for i = 1:size(pos,2)-1
        B_N_sim[i,1] = mag_field((lat[i]+π/2)/(π)*size(mag_field,1),(long[i]+π)/(2*π)*size(mag_field,1),1)
        B_N_sim[i,2] = mag_field((lat[i]+π/2)/(π)*size(mag_field,1),(long[i]+π)/(2*π)*size(mag_field,1),2)
        B_N_sim[i,3] = mag_field((lat[i]+π/2)/(π)*size(mag_field,1),(long[i]+π)/(2*π)*size(mag_field,1),3)
        # B_N_sim[i,:] = SatelliteToolbox.igrf12(2019,(alt+R_E)*1000,lat[i],long[i])/1.e9
    end
    B_N_sim = convert(Array{Float64},B_N_sim)
    # plot(sol.t[1:end-1],B_N_sim[:,1])
    # plot!(sol.t[1:end-1],B_N_sim[:,2])
    # plot!(sol.t[1:end-1],B_N_sim[:,3])

    return B_N_sim

end
