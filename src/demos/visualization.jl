using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using MeshCat
using GeometryTypes
using CoordinateTransformations
using FileIO
using MeshIO
using Random
using Rotations
using StaticArrays
using GeometryTypes: GLUVMesh # we need a mesh type that stores texture coordinates
using Colors

#assuming that you already have a trajectory you want to plot....
using HDF5
X_sim = h5read("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/visualize.h5", "X_sim")
B_ECI = h5read("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/visualize.h5", "B_ECI")
N = size(X_sim,2)
#################
# Visualization #
#################

vis = Visualizer()
open(vis)

# color options
green_ = MeshPhongMaterial(color=RGBA(0, 1, 0, 1.0))
# green_transparent = MeshPhongMaterial(color=RGBA(0, 1, 0, 0.1))
red_ = MeshPhongMaterial(color=RGBA(1, 0, 0, 1.0))
# red_transparent = MeshPhongMaterial(color=RGBA(1, 0, 0, 0.1))
blue_ = MeshPhongMaterial(color=RGBA(0, 0, 1, 1.0))
# blue_transparent = MeshPhongMaterial(color=RGBA(0, 0, 1, 0.1))
# blue_semi = MeshPhongMaterial(color=RGBA(0, 0, 1, 0.5))
# yellow_ = MeshPhongMaterial(color=RGBA(1, 1, 0, 1.0))
# yellow_transparent = MeshPhongMaterial(color=RGBA(1, 1, 0, 0.75))
#
orange_ = MeshPhongMaterial(color=RGBA(233/255, 164/255, 16/255, 1.0))
# orange_transparent = MeshPhongMaterial(color=RGBA(233/255, 164/255, 16/255, 0.1))
# black_ = MeshPhongMaterial(color=RGBA(0, 0, 0, 1.0))
# black_transparent = MeshPhongMaterial(color=RGBA(0, 0, 0, 0.1))
# black_semi = MeshPhongMaterial(color=RGBA(0, 0, 0, 0.5))

# create 1U cubesat object
image = PngImage(joinpath(MeshCat.VIEWER_ROOT, "..", "data", "cubesatside.png"))
texture = Texture(image=image)
antenna_image = PngImage(joinpath(MeshCat.VIEWER_ROOT, "..", "data", "antenna.png"))
antenna2_image = PngImage(joinpath(MeshCat.VIEWER_ROOT, "..", "data", "antenna_flipped.png"))
antenna_texture = Texture(image=antenna_image)
antenna2_texture = Texture(image=antenna2_image)
material = MeshLambertMaterial(map=texture)
antenna_material = MeshLambertMaterial(map=antenna_texture)
antenna2_material = MeshLambertMaterial(map=antenna2_texture)
satellite = vis["satellite"]
obj_sat = setobject!(vis["satellite"], HyperRectangle(Vec(-.5, -.5, -.5), Vec(1, 1, 1)),material)
antenna = vis["antenna"]
antenna2 = vis["antenna2"]
obj_antenna = setobject!(vis["antenna"], HyperRectangle(Vec(0,0,0),Vec(.125,1,0)),antenna_material)
obj_antenna = setobject!(vis["antenna2"], HyperRectangle(Vec(0,0,0),Vec(.125,1,0)),antenna_material)
settransform!(vis["antenna"], compose(Translation(0,.5,-.1),LinearMap(Quat(0,0,-cosd(45/2),sind(45/2)))))
settransform!(vis["antenna2"], compose(Translation(0,.5,.1),LinearMap(Quat(0,0,-cosd(-45/2),sind(-45/2)))))

#create axes
axis1 = vis["axis1"]
magnetic_vector = vis["axis1"]
obj_mag = setobject!(vis["axis1"],
    Cylinder{3, Float32}(
         Point3f0(0., 0., 0.),
         Point3f0(0., 0., 1.),
         Float64(.01)),red_)
axis2 = vis["axis2"]
magnetic_vector = vis["axis2"]
obj_mag = setobject!(vis["axis2"],
 Cylinder{3, Float32}(
      Point3f0(0., 0., 0.),
      Point3f0(0., 1., 0.),
      Float64(.01)),green_)
axis3 = vis["axis3"]
magnetic_vector = vis["axis3"]
obj_mag = setobject!(vis["axis3"],
  Cylinder{3, Float32}(
       Point3f0(0., 0., 0.),
       Point3f0(1., 0., 0.),
       Float64(.01)),blue_)




# create very long thin cylinder pointing for magnetic field
magnetic_vector = vis["magnetic_vector"]
using LinearAlgebra
obj_mag = setobject!(vis["magnetic_vector"],
    Cylinder{3, Float64}(
         Point3f0(0., 0., 0.),
         Point3f0(2*B_ECI[1,1]/norm(B_ECI[1,:]),
         2*B_ECI[1,2]/norm(B_ECI[1,:]),
         2*B_ECI[1,3]/norm(B_ECI[1,:])),
         Float64(.01)),orange_)
# setobject!(vis["magnetic_vector"],
#      Cylinder{3, Float64}(
#           Point3f0(0., 0., 0.),
#           Point3f0(B_ECI[N,1]/norm(B_ECI[N,:]),
#           B_ECI[N,2]/norm(B_ECI[N,:]),
#           B_ECI[N,3]/norm(B_ECI[N,:])),
#           Float64(.01)))



# Set camera location
settransform!(vis["/Cameras/default"],LinearMap(Quat(X_sim[4:7,[1]]...)))
setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", .1)
setprop!(vis["/Grid"], "visible", false)
# setprop!(vis["/Background"],"color",MeshPhongMaterial(color=RGBA(0, 0, 0, 0.5)))

# # Animate satellite
# time_multiplier = 100 #goes this times as fastb
# for i = 1:N
#     settransform!(vis["satellite"],compose(Translation(0,0,0),LinearMap(quat2rot(X[4:7,[i]]))))
#     obj_mag = setobject!(vis["magnetic_vector"],
#         Cylinder{3, Float64}(
#              Point3f0(0., 0., 0.),
#              Point3f0(B_ECI[i,1]/norm(B_ECI[i,:]),
#              B_ECI[i,2]/norm(B_ECI[i,:]),
#              B_ECI[i,3]/norm(B_ECI[i,:])),
#              Float64(.01)))
#      sleep((t_final-t0)/N/time_multiplier)
#
# end

#find the transforms between the magnetic field vectors
B_ECI_shift = zeros(4,N)
for i = 1:N-1
    axis = cross(B_ECI[1,:],B_ECI[i+1,:])
    axis *= norm(axis)^(-1)
    angle = acos(dot(B_ECI[1,:],B_ECI[i+1,:])/norm(B_ECI[1,:])/norm(B_ECI[i+1,:]))
    quat = [cos(angle/2);axis*sin(angle/2)]
    quat *= norm(quat)^(-1)
    B_ECI_shift[:,i] = quat
end

#save animation
anim = MeshCat.Animation(30)
for i = 1:size(X_sim,2)-1
    MeshCat.atframe(anim,vis,i) do frame
        #transform the satellite and the axes
        settransform!(frame["satellite"], compose(Translation(0,0,0),LinearMap(Quat(X_sim[4:7,[i]]...))))
        settransform!(frame["axis1"], compose(Translation(0,0,0),LinearMap(Quat(X_sim[4:7,[i]]...))))
        settransform!(frame["axis2"], compose(Translation(0,0,0),LinearMap(Quat(X_sim[4:7,[i]]...))))
        settransform!(frame["axis3"], compose(Translation(0,0,0),LinearMap(Quat(X_sim[4:7,[i]]...))))

        # antennas
        antenna_quat = Quat(qmult(X_sim[4:7,[i]],[0,0,-cosd(45/2),sind(45/2)])...)
        antenna_quat2 = Quat(qmult(X_sim[4:7,[i]],[0,0,-cosd(-45/2),sind(-45/2)])...)

        new_point = antenna_quat * SVector(0,.5,-.1)
        new_point2 = antenna_quat2 * SVector(0,.5,.1)

        settransform!(frame["antenna"], compose(Translation(new_point),LinearMap(antenna_quat)))
        settransform!(frame["antenna2"], compose(Translation(new_point2),LinearMap(antenna_quat2)))

        #transform the magnetic vector
        settransform!(frame["magnetic_vector"], compose(Translation(0,0,0),LinearMap(Quat(B_ECI_shift[:,i]...))))
    end
end
MeshCat.setanimation!(vis,anim)

# put earth in there
#create Earth
length_scale = 1e0
earth_image = PngImage(joinpath(MeshCat.VIEWER_ROOT, "..", "data", "earth.png"))
earth_texture = Texture(image=earth_image)
earth_material = MeshLambertMaterial(map=earth_texture)
earth = HyperSphere(Point(pos[:,1][1]*length_scale, pos[:,1][2]*length_scale, pos[:,1][3]*length_scale)
    , 6731*length_scale)
setobject!(vis["earth"], earth, earth_material)
settransform!(vis["earth"], compose(LinearMap(AngleAxis(2*pi, 0, 0, 1)),
    LinearMap(AngleAxis(pi, 1, 0, 0)), ))


pwd()
joinpath(pwd(),"sat_frames")
using FileIO
using ImageMagick
using Images, ImageView

anim = @animate for _img in readdir(_path)
    im_ = FileIO.load(joinpath(_path,_img))
    im_crop = im_[250:end-50,400:end-400]

    plot(im_crop,axis=false)
end
gif(anim, joinpath(_path,"maze_v2.gif"), fps = 20)


# settransform!(vis["/Cameras/default"], compose(Translation(0., 75., 50.),LinearMap(RotX(pi/10)*RotZ(pi/2))))
# settransform!(vis["/Cameras/default"], compose(Translation(0., 35., 65.),LinearMap(RotX(pi/3)*RotZ(pi/2))))

# Create and place trajectory
# for i = 1:N
#     setobject!(vis["traj_uncon"]["t$i"],sphere_small,blue_)
#     settransform!(vis["traj_uncon"]["t$i"], Translation(results_uncon.X[i][1], results_uncon.X[i][2], results_uncon.X[i][3]))
# end

# for i = 1:N
#     setobject!(vis["traj_x0"]["t$i"],sphere_small,blue_)
#     settransform!(vis["traj_x0"]["t$i"], Translation(X0[1,i], X0[2,i], X0[3,i]))
# end
# for i = 1:size(X_guess,2)
#     setobject!(vis["traj_x0"]["t$i"],sphere_medium,yellow_transparent)
#     settransform!(vis["traj_x0"]["t$i"], Translation(X_guess[1:3,i]...))
# end

# for i = 1:N
#     setobject!(vis["traj"]["t$i"],obj,black_)
#     settransform!(vis["traj"]["t$i"], Translation(results_con.X[i][1:3]...))
# end

# setobject!(vis["robot_uncon"]["ball"],sphere_medium,orange_transparent)
# setobject!(vis["robot_uncon"]["quad"],robot_obj,black_)

# setobject!(vis["robot"]["ball"],sphere_medium,green_transparent)
# setobject!(vis["robot"]["quad"],robot_obj,black_)

# Animate quadrotor
# for i = 1:N
#     settransform!(vis["robot_uncon"], compose(Translation(results_uncon.X[i][1], results_uncon.X[i][2], results_uncon.X[i][3]),LinearMap(Quat(results_uncon.X[i][4:7]...))))
#     sleep(solver_uncon.dt)
# end

# Ghose quadrotor scene
# traj_idx = [1;12;20;30;40;50;N]
# n_robots = length(traj_idx)
# for i = 1:n_robots
#     robot = vis["robot_$i"]
#     setobject!(vis["robot_$i"]["quad"],robot_obj,black_semi)
#     settransform!(vis["robot_$i"], compose(Translation(results_con.X[traj_idx[i]][1], results_con.X[traj_idx[i]][2], results_con.X[traj_idx[i]][3]),LinearMap(quat2rot(results_con.X[traj_idx[i]][4:7]))))
# end
#
# # Animation
# # (N-1)/tf
# anim = MeshCat.Animation(20)
# for i = 1:N
#     MeshCat.atframe(anim,vis,i) do frame
#         settransform!(frame["robot"], compose(Translation(results_con.X[i][1:3]...),LinearMap(Quat(results_con.X[i][4:7]...))))
#     end
# end
# MeshCat.setanimation!(vis,anim)
