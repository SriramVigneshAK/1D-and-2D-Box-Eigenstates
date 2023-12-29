#This is a code to calculate the eigenstates of a 1-Dimensional box
using QuantumOptics
using Plots

#Mass of particle and length of box
m = 1 #mass
l = 10 #length
n = 1 #Quantum number

#Defining a Basis to represent everything
xmin = -l
xmax = l
Npoints = 100
b_position = PositionBasis(xmin, xmax, Npoints) #Position Basis
b_momentum = MomentumBasis(b_position) #Momentum Basis

#Defining Operators 
x = position(b_position) #Position Operator in Position Basis
p = momentum(b_momentum) #Momentum Operator in Momentum Basis
px = momentum(b_position) #Momentum Operator in Position Basis

#Defining the Potential 
function V(x)
    if abs(x)>=5
        v = 1000
    else
        v=0
    end
    return v
end

#Plotting the Potential
PE = potentialoperator(b_position, V)
PE = dense(PE)
ptsx = samplepoints(b_position)
#plot(ptsx , V)


#Defining the Hamiltonian
H = LazySum(px^2/2m , PE)
# A more feasible way is to use Fourier transforms(FT), so defining FT operators
Txp = transform(b_position, b_momentum)
Tpx = transform(b_momentum, b_position)
Hkin = LazyProduct(Txp, p^2/2m, Tpx) #A faster method to do Hkin * psi
H = dense(LazySum(Hkin, PE))

#finding eigenstates
E,ψ = eigenstates(H) #Computes all  eigenstates
En,ψn = eigenstates(H,n) #Computes only the required eigenstate
ψn
#Plotting the States
#plot!(ptsx, 500*real(ψ[n].data)) #Can be plotted from all eigenstates 
plot(ptsx, [500*real(ψn[1].data) , V], title="1D Box n=$n eigenstate", label=["ψ(x)" "V(x)"])
xlabel!("x")
ylabel!("ψ(x)")
print("The Energy of the corresponding state is  : $(real(E[n])) ")
