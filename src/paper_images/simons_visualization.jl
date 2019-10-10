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

#################
# Visualization #
#################

vis = Visualizer()
open(vis)

# color options
green_ = MeshPhongMaterial(color=RGBA(0, 1, 0, 1.0))
red_ = MeshPhongMaterial(color=RGBA(1, 0, 0, 1.0))
blue_ = MeshPhongMaterial(color=RGBA(0, 0, 1, 1.0))
orange_ = MeshPhongMaterial(color=RGBA(233/255, 164/255, 16/255, 1.0))

#set length scale
length_scale = 1E-3

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
obj_sat = setobject!(vis["satellite"], HyperRectangle(Vec(-.5, -.5, -.5)*length_scale, Vec(1, 1, 1)*length_scale),material)
antenna = vis["antenna"]
antenna2 = vis["antenna2"]
obj_antenna = setobject!(vis["antenna"], HyperRectangle(Vec(0,0,0)*length_scale,Vec(.125,1,0)*length_scale),antenna_material)
obj_antenna = setobject!(vis["antenna2"], HyperRectangle(Vec(0,0,0)*length_scale,Vec(.125,1,0)*length_scale),antenna_material)
settransform!(vis["antenna"], compose(Translation([0,.5,-.1]*length_scale),LinearMap(Quat(0,0,-cosd(45/2),sind(45/2)))))
settransform!(vis["antenna2"], compose(Translation([0,.5,.1]*length_scale),LinearMap(Quat(0,0,-cosd(-45/2),sind(-45/2)))))

#create Earth
earth_image = PngImage(joinpath(MeshCat.VIEWER_ROOT, "..", "data", "earth.png"))
earth_texture = Texture(image=earth_image)
earth_material = MeshLambertMaterial(map=earth_texture)
earth = HyperSphere(Point(0.0,0.0,0.0), 6731. * length_scale)
setobject!(vis["earth"], earth, earth_material)
settransform!(vis["earth"], compose(LinearMap(AngleAxis(-0.8*pi, 0, 0, 1)),
    LinearMap(AngleAxis(pi/2, 1, 0, 0)), ))

#create axes
axis1 = vis["axis1"]
magnetic_vector = vis["axis1"]
obj_mag = setobject!(vis["axis1"],
    Cylinder{3, Float32}(
         Point3f0(0., 0., 0.)*length_scale,
         Point3f0(0., 0., 1.)*length_scale,
         Float64(.01)*length_scale),red_)
axis2 = vis["axis2"]
magnetic_vector = vis["axis2"]
obj_mag = setobject!(vis["axis2"],
 Cylinder{3, Float32}(
      Point3f0(0., 0., 0.)*length_scale,
      Point3f0(0., 1., 0.)*length_scale,
      Float64(.01)*length_scale),green_)
axis3 = vis["axis3"]
magnetic_vector = vis["axis3"]
obj_mag = setobject!(vis["axis3"],
  Cylinder{3, Float32}(
       Point3f0(0., 0., 0.)*length_scale,
       Point3f0(1., 0., 0.)*length_scale,
       Float64(.01)*length_scale),blue_)


# create very long thin cylinder pointing for magnetic field
magnetic_vector = vis["magnetic_vector"]
obj_mag = setobject!(vis["magnetic_vector"],
    Cylinder{3, Float64}(
         Point3f0(0., 0., 0.)*length_scale,
         Point3f0(2*B_ECI[1,1]/norm(B_ECI[1,:]),
         2*B_ECI[1,2]/norm(B_ECI[1,:]),
         2*B_ECI[1,3]/norm(B_ECI[1,:])*length_scale),
         Float64(.01)*length_scale),orange_)
# setobject!(vis["magnetic_vector"],
#      Cylinder{3, Float64}(
#           Point3f0(0., 0., 0.),
#           Point3f0(B_ECI[N,1]/norm(B_ECI[N,:]),
#           B_ECI[N,2]/norm(B_ECI[N,:]),
#           B_ECI[N,3]/norm(B_ECI[N,:])),
#           Float64(.01)))



# Set camera location
settransform!(vis["/Cameras/default"],compose(Translation(0,0,0),LinearMap(Quat(X[4:7,[1]]...))))
setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", .1)
setprop!(vis["/Grid"], "visible", false)
# setprop!(vis["/Background"],"color",MeshPhongMaterial(color=RGBA(0, 0, 0, 0.5)))

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
for i = 1:size(U_sim,2)-1
    MeshCat.atframe(anim,vis,i) do frame
        #transform the satellite and the axes
        sat_point = LinearMap(Quat(X_sim[4:7,[i]]...)) * SVector(pos[:,i])
        settransform!(frame["satellite"], compose(Translation(pos[:,i]),LinearMap(Quat(X_sim[4:7,[i]]...))))
        settransform!(frame["axis1"], compose(Translation(pos[:,i]),LinearMap(Quat(X_sim[4:7,[i]]...))))
        settransform!(frame["axis2"], compose(Translation(pos[:,i]),LinearMap(Quat(X_sim[4:7,[i]]...))))
        settransform!(frame["axis3"], compose(Translation(pos[:,i]),LinearMap(Quat(X_sim[4:7,[i]]...))))

        # # antennas
        # antenna_quat = compose(Translation(pos[:,i]),LinearMap(Quat(qmult(X_sim[4:7,[i]],[0,0,-cosd(45/2),sind(45/2)])...)))
        # antenna_quat2 = compose(Translation(pos[:,i]),LinearMap(Quat(qmult(X_sim[4:7,[i]],[0,0,-cosd(-45/2),sind(-45/2)])...)))
        # new_point = antenna_quat * SVector(0,.5,-.1)
        # new_point2 = antenna_quat2 * SVector(0,.5,.1)
        # settransform!(frame["antenna"], compose(Translation(new_point),LinearMap(antenna_quat)))
        # settransform!(frame["antenna2"], compose(Translation(new_point2),LinearMap(antenna_quat2)))

        #transform the magnetic vector
        settransform!(frame["magnetic_vector"], compose(Translation(0,0,0),LinearMap(Quat(B_ECI_shift[:,i]...))))

        #transform the camera
        camera_transformation = Translation(pos[:,i]*length_scale+[1,1,1]*length_scale)
        settransform!(vis["/Cameras/default"], camera_transformation)

    end
end
MeshCat.setanimation!(vis,anim)


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
