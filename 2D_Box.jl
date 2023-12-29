#This is a code to calculate the eigenstates of a 2-Dimensional box
using QuantumOptics
using Plots
#using PyPlot


#Mass of particle and length of box
m = 1 #mass
l = 10 #length
nx = 1 #Quantum number
ny = 1

#Defining a Basis to represent everything
#We have to define two bases, one along x and one more along y
xmin = -l
xmax = l
Npointsx = 100
b_positionx = PositionBasis(xmin, xmax, Npointsx) #Position Basis
b_momentumx = MomentumBasis(b_positionx) #Momentum Basis

ymin = -l
ymax = l
Npointsy = 100
b_positiony = PositionBasis(ymin, ymax, Npointsy) #Position Basis
b_momentumy = MomentumBasis(b_positiony) #Momentum Basis

#We have to define a Composite Bases to Represent everything
#This is done by taking direct product of both bases 
bcomp_x = b_positionx ⊗ b_positiony
bcomp_p = b_momentumx ⊗ b_momentumy

# A more feasible way is to use Fourier transforms(FT), so defining FT operators
Txp = transform(bcomp_x, bcomp_p)
Tpx = transform(bcomp_p, bcomp_x)

#Defining Operators 
x = position(b_positionx) #Position Operator in Position Basis
p1 = momentum(b_momentumx) #Momentum Operator in Momentum Basis
px = momentum(b_positionx) #Momentum Operator in Position Basis

y = position(b_positiony) #Position Operator in Position Basis
p2 = momentum(b_momentumy) #Momentum Operator in Momentum Basis
py = momentum(b_positiony) #Momentum Operator in Position Basis

#Defining the Potential 
function V(x,y)
    if abs(x)>=5 && abs(y)>=5
        v = 1000
    else
        v=0
    end
    return v
end

#Plotting the Potential
PE = potentialoperator(bcomp_x, V)
PE = dense(PE)
ptsx, ptsy = samplepoints(b_positionx), samplepoints(b_positiony)
plot(ptsx , V , color = "red")
contour(ptsx, ptsy,V, fill = true)


#Defining the Hamiltonian
#H = LazySum(px^2/2m , PE)
p=p1⊗p2
Hkin = LazyProduct(Txp, p^2/2m, Tpx) #A faster method to do Hkin * psi
H = dense(LazySum(Hkin, PE))

#finding eigenstates
E,ψ = eigenstates(H) #Computes all  eigenstates
En,ψn = eigenstates(H,nx,ny) #Computes only the required eigenstate
ψn
#Plotting the States
#contour(ptsx,ptsy, 500*real(ψ[n].data)) #Can be plotted from all eigenstates 
contour(ptsx, ptsy, 100*real(ψ[3].data))
plot(ptsx, [500*real(ψn[1].data) , V], title="1D Box n=$n eigenstate", label=["ψ(x)" "V(x)"])
xlabel!("x")
ylabel!("ψ(x)")
print("The Energy of the corresponding state is  : $(real(E[1])) ")
