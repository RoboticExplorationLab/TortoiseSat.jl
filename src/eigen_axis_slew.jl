function eigen_axis_slew(x0,xf,t)
	# This function calculates an eigen axis slew for a given
	# maneuver and time vector 
	# INPUTS:
	# x0: starting attitude in [ω;q]
	# xf: final attitude in [ω;q]
	# t: time vector

	# OUTPUTS:
	# ω_guess - guess of omega at each time step
	# q_guess: guess of quaternion at each time step

	# first, given both quaternions, find the angle/axis
	q1 = x0[4:7]
	q2 = xf[4:7]
	q_e = qmult([q2;-q2[2:4]],q1)
	theta_f = 2*acos(q_e[1])
	axis = -q_e[2:4]/(sin(theta_f/2))

	# find theta over time using versine
	alpha = pi/t[end]
	theta = theta_f*1/2*(ones(length(t))-cos.(alpha*t))
	d_theta = diff(theta)/(t[2]-t[1])
	push!(d_theta,d_theta[end])

	ω_guess = zeros(length(t),3)
	for i in 1:length(t)
		ω_guess[i,:] = d_theta[i]*axis
	end

	# now find the quaternion 
	q_guess = zeros(length(t),4)
	for i in 1:length(t)
		q_guess[i,:] = qmult(q1,[cos(theta[i]/2);axis*sin(theta[i]/2)])
	end

	return ω_guess, q_guess
end


function qmult(q1,q2)
      [q1[1]*q2[1] - q1[2:4]'*q2[2:4]; q1[1]*q2[2:4] + q2[1]*q1[2:4] + cross(q1[2:4],q2[2:4])]
end

