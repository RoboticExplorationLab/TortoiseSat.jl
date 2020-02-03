function simulator(dx,x,u)
	#test function once the input has been determined by the trajectory optimization package
	#assign noise to rotation

	omega_noise = randn(3)*(.38*pi/180)^2
	omega = x[1:3]+omega_noise

	#assign noise to quaternion

	q_noise = randn(3,1)*(1*pi/180)^2
	θ_noise = norm(q_noise)
	r_noise = q_noise/θ_noise
	q = qmult(x[4:7]./norm(x[4:7]),[cos(θ_noise/2);r_noise*sin(θ_noise/2)])
	t=x[8]

	skew_omega=[0 omega[3] -omega[2];
	   -omega[3] 0 omega[1];
	   omega[2] -omega[1] 0]

	#find quaternion derivative
	q_dot=0.5*qmult(q,[0; omega])
	B_N_noise = rand(3)*(1E-5)^2 #adding 10 percent error
	temp=qrot(q,B_ECI[floor(Int,t*N+1),:]+B_N_noise)
	# temp=qrot(q,B_N[1,:])
	B_B=temp[1:3]

	#m_c=zeros(3);
	#for i=1:3
	#   temp=cross(-omega,B_B1);
	#   m_c[i]=-m_max*sign(temp[i]);
	#end

	tau_c=cross(u[1:3]/100,B_B) #N-m
	# tau_c=u[1:3]

	#rotation
	omega_dot=inv(p.J)*(tau_c-cross(omega,p.J*omega)) #rad/s^2

	#redo concatination
	dx[1:8] = [omega_dot;q_dot;1/(tf-t0)]

end

function qrot(q,r)
      r + 2*cross(q[2:4],cross(q[2:4],r) + q[1]*r)
end

function qmult(q1,q2)
      [q1[1]*q2[1] - q1[2:4]'*q2[2:4]; q1[1]*q2[2:4] + q2[1]*q1[2:4] + cross(q1[2:4],q2[2:4])]
end
