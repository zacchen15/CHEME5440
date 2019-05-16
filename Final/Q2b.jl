# constants
n = 2
a = 10

# create range for each repressor concentration
uStart = 0
uStep = 0.1
uStop = 20.0
uRange = collect(uStart:uStep:uStop)

vStart = 0
vStep = 0.1
vStop = 20.0
vRange = collect(vStart:vStep:vStop)

# initialize the array to store the repressor concentrations for the nullclines
u = Array{Float64}(undef,length(vRange),1)
v = Array{Float64}(undef,length(uRange),1)

# calculate the nullclines for repressor 1
for i in 1:length(vRange)
    u[i] = (a / (1 + vRange[i]^n))
end
# calculate the nullclines for repressor 2
for i in 1:length(uRange)
    v[i] = (a / (1 + uRange[i]^n))
end

using PyPlot
close("all")
figure(1)
plot(u,vRange,color="blue",linestyle="solid",label="Repressor 1")
plot(uRange,v,color="red",linestyle="solid",label="Repressor 2")
axvline(x=0,color="black",linestyle=":")
axhline(y=0,color="black",linestyle=":")
title("Nullclines for n = 2")
legend(loc="upper right")
xlabel("u")
ylabel("v")
