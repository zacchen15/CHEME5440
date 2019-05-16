# constants
n = 1
a = 10

# create the data
minval = 0.0
maxval = 20
steps = 200
u = repeat(range(minval,stop=maxval,length=steps)',steps)
v = repeat(range(minval,stop=maxval,length=steps),1,steps)
f = (a ./ (1 .+ v.^n)) .- u
g = (a ./ (1 .+ u.^n)) .- v

using PyPlot
close("all")
figure(1)
streamplot(u,v,f,g,density=1)
title("Streamplot for n = 2")
xlabel("u")
ylabel("v")
