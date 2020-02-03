# this file simplifies the command line by introducing
# the parameters all at once based on options 

struct params
	type::AbstractString
	mass::Float64
	J::Array{Float64}
	BC::Float64
	alt::Float64
	Kep::Array{Float64}
	MJD::Float64
	GM::Float64
	R_E::Float64
	T::Float64
	ω_0::Float64
end

# INPUTS:
# - type of satellite (1P,1U)
# - Keplerian orbital elements:
# 	- e, a, i, RAAN, ω, ν (true anomaly)
# - Modified Julian day

function input_parameters(type,Kep,MJD)
	#earth parameters
	GM = 3.986004418E14*(1/1000)^3 #earth graviational constant (km^3 s^-2)""
	R_E = 6371.0 #earth radius (km)

	if type == "1U"
		#satellite parameters (1U)
		mass=.75 #kg
		J=[00.00125 0 0;
		    0 0.00125 0;
		    0 0 0.00125] #kgm2
		BC = mass/2.2/(J[1,1]*J[2,2])

	elseif type == "1P"
		#satellite parameters (1P)
		mass = .25 #kg
		J =[0.0001041667 0 0;
		    0 0.0001041667 0;
		    0 0 0.0001041667] #kgm2
		BC = mass/2.2/(J[1,1]*J[2,2])

	else 
		println("Type not recognized")
	end

	#enter orbital parameters
	alt=400 #altitude in km
	T=2*pi*sqrt(Kep[2].^3/GM)
	ω_0 = sqrt(GM/Kep[2]^3) #rad/s

	#Modified Julia Date
	MJD = 58155.0 #modified julian day

	return params(type,mass,J,BC,alt,Kep,MJD,GM,R_E,T,ω_0)
end

