#does not include all of the packages required

#earth parameters
GM = 3.986004418E14*(1/1000)^3 #earth graviational constant (km^3 s^-2)""
R_E = 6371.0 #earth radius (km)

#satellite parameters (1U)
mass=.75 #kg
J=[00.00125 0 0;
    0 0.00125 0;
    0 0 0.00125] #kgm2
J_inv=inv(J)
BC = mass/2.2/(J[1,1]*J[2,2])

alt=400 #altitude in km
A=zeros(1,6)
#enter orbital parameters
A[1] = 0 #circular orbits
A[2]= alt+R_E #semi-major axis
# A[3] = rand(1)[1]*90 #inclination in degrees
A[3] = 90 #inclination in degrees
A[4] = 0 #assign ascending nodes
A[5] = 0 #argument of periapsis
A[6] = 0 #initial true anomoly
T=2*pi*sqrt(A[2].^3/GM)
ω_0 = sqrt(GM/A[2]^3) #rad/s

#create vector of final times
t_final = zeros(10)

#Modified Julia Date
MJD_0 = 58155.0 #modified julian day

#magnetic field generation for half orbit
t0 = 0.0 #seconds
tf = 60*90 #seconds
cutoff = collect(range(30,length = length(t_final),100))
N = 5000

#extrapolate over the magnetic field
alt = 400 # km
mag_field = igrf_data(alt,2019)

for i in 1:length(t_final)
    B_ECI_initial,a,b = magnetic_simulation(A,t0,tf,N,mag_field,GM,MJD_0)

    # Find magnetic control garmian
    B_gram = magnetic_gramian(B_ECI_initial,(tf-t0)/N)

    #determine at what point the condition cutoff is satisfied
    tf_index = condition_based_time(B_gram,cutoff[i])

    t_final[i] = tf_index*(tf-t0)/N
    println("Finished ",i," out of ",length(t_final)," iterations.")
end

using Plots.PlotMeasures
plt = plot(cutoff,t_final,xlabel = "Condition Number Requirement",
    ylabel = "Time Required (s)",left_margin = 20px,legend = false)
Plots.savefig("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/cutoff_requirements.png")
