# calculate the eiganvalues of the Jacobian
# constants
n = 2
a = 10

# fixed point of interest
u = 2
v = 2

# evaluate each partial derivative
fu = -1
fv = (a * n * v^(n-1)) / (1 + v^n)^2
gu = (a * n * u^(n-1)) / (1 + u^n)^2
gv = -1

# Jacobian
J = [fu fv; gu gv]

using LinearAlgebra
# calculate eigenvalues
lambda1 = (tr(J) + sqrt(tr(J)^2 - 4*det(J))) / 2
lambda2 = (tr(J) - sqrt(tr(J)^2 - 4*det(J))) / 2
println("lambda1 = ",lambda1)
println("lambda2 = ",lambda2)

if (lambda1 < 0) && (lambda2 < 0)
    println("Fixed Point is Stable")
else
    println("Fixed Point is Unstable")
end
