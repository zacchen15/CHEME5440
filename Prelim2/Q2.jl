using LinearAlgebra

# Tunable parameter
DaArray = collect(0.0:1.0:10.0)
store = Array{Float64}(undef,6,length(DaArray))

for i in 1:length(DaArray)
        # Possible steady state solution
        Da = DaArray[i]
        A = 1.0
        B = 2.0-Da
        # Linearizations
        fA = -Da-B
        fB = -A
        gA = B
        gB = 2*A-1
        # Create Jacobian
        J = Array{Float64}(undef,2,2)
        J[1,:] = [fA fB]
        J[2,:] = [gA gB]
        # Solve for the things we want and store them
        store[1,i] = Da
        store[2,i],store[3,i] = eigvals(J)
        store[4,i] = tr(J)
        store[5,i] = det(J)
        store[6,i] = 0.25*(tr(J))^2
end
