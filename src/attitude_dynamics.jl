
function attitude_dynamics(x,u,B_B,J)
      #this function returns ẋ given x, an input of magnetic moment and a BF mag vector

      omega=x[1:3]
      q=x[4:7]./norm(x[4:7])

      skew_omega=[0 omega[3] -omega[2];
         -omega[3] 0 omega[1];
         omega[2] -omega[1] 0]

      #find quaternion derivative
      q_dot=0.5*qmult(q,[0; omega])

      tau_c=cross(u[1:3],B_B) #N-m

      #rotation
      omega_dot=inv(J)*(tau_c-cross(omega,J*omega)) #rad/s^2

      #redo concatination
      ẋ = [omega_dot;q_dot]
      #print(omega);

end

function attitude_dynamics_linear(x,u,x_linear, B_B,J)
      #this function returns ẋ given x, an input of magnetic moment and a BF mag vector

      omega=x[1:3]
      q=x[4:7]./norm(x[4:7])

      skew_omega=[0 omega[3] -omega[2];
         -omega[3] 0 omega[1];
         omega[2] -omega[1] 0]

      #find quaternion derivative
      q_dot=0.5*qmult(q,[0; x_linear[4:6]])

      tau_c=cross(u[1:3],B_B) #N-m

      #rotation
      omega_dot=inv(J)*(tau_c-cross(omega,J*omega)) #rad/s^2

      #redo concatination
      ẋ = [omega_dot;q_dot]
      #print(omega);

end


function qrot(q,r)
      (r + 2*cross(q[2:4],cross(q[2:4],r) + q[1]*r))
end

function qmult(q1,q2)
      [q1[1]*q2[1] - q1[2:4]'*q2[2:4]; q1[1]*q2[2:4] + q2[1]*q1[2:4] + cross(q1[2:4],q2[2:4])]
end
