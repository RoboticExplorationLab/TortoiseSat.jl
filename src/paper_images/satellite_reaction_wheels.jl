# generate an image showing the size of
# reaction wheels versus the size of sat
using Pkg
using LinearAlgebra
using Plots

# max dimension (m), torque (Nm)
u = [.1 100*.1*40E-6*6/1000*.1^-3]
uuu = [.166667 100*.1*40E-6*.1*.3/(1/12*1000*.3*.1*.1*(.1*.3))]
uuuuuu = [(.2+.2+.3)/3 100*.1*40E-6*.2*.3/(1/12*1000*.3*.2*.2*(.2*.3))]
sizes = [u;uuu;uuuuuu]

# general curve
x = range(.05,.25; step = .001)
y = 100*.15*40E-6*6/1000*x.^-3

plt = plot(x,y,xflip = false,yscale = :log,legend = false)
xlabel!("Average Satellite Dimension (m)")
ylabel!("Angular Acceleration (1/s^2)")
scatter!(sizes[:,1],sizes[:,2],series_annotations = text.(["1U","3U","6U"], :top))

Plots.savefig("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/satellite_scaling.png")

# satellite mass fractions
masses = [1.33 3*1.33 6*1.33 110 3000 420000]
acs_fraction = [.366 .226 .164 .109 .02 .0095]

# general curve
x = range(1,500000; step = .1)
y = x.^-.8
plt = plot(x,y,xflip = false,xscale = :log,legend = false)
scatter(masses[:], acs_fraction[:], series_annotations = text.(["  1U","  3U","  6U","  LEO SmallSat","  GEO Sat","  ISS"], :left),
    xscale = :log, color = :red, legend = false, xlims = (1,10^6.5))
ylabel!("ACS Mass Fraction")
xlabel!("Satellite Mass (kg)")
Plots.savefig("/Users/agatherer/Desktop/Academic-Papers/Manchester/Data/acs_mass_fraction.png")
