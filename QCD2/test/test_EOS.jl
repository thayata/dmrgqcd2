#test code
using ITensors
using Plots
include("../src/observables.jl")
include("../src/hamiltonian.jl")

######################
# demo to plot the μ dependence of observables  #
######################

# parameters
N = 40
g = 1.0
M = 1.
μ = 0
Λ = 400.
J = 0.5*g*g
w= 1 
Nc = 2
Nf = 1

sites = siteinds_qcd2(N, Nc, Nf, false)

sweeps = Sweeps(40) # number of sweeps
maxdim!(sweeps,20,20,40,40,60,60,80,80,100,100,200,200,400,400,400,400) # gradually increase states kept
cutoff!(sweeps,1E-8) # desired truncation error


H = hamiltonian(sites, w, J, M, μ, Λ, Nc, Nf, "fnc")
psi0 = randomMPS(sites,10)
E0, psi = dmrg(H,psi0,sweeps)

Δμ = 0.01
μseq = 1.1:Δμ:1.7
res_pressure = []
res_energy = []
res_Nq = []
res_sigma = []
for μ in μseq
    println(μ)
    H = hamiltonian(sites, w, J, M, μ, Λ, Nc, Nf, "fnc")
    energy,psi = dmrg(H,psi,sweeps)
    #pseq, distribution = quark_distribution(sites, psi, 1, 1, Nc, Nf) # this function is time-consuming
    Nq = quark_number(psi, Nc, Nf)
    push!(res_pressure,E0-energy)
    push!(res_energy,energy-E0+μ*sum(Nq))
    push!(res_Nq, Nq)
    push!(res_sigma,-chiral_condensate(psi, Nc, Nf))
end


plot(μseq, res_Nq,
label        ="quark number",
linecolor    =:red,       
)


plot(res_Nq)
plot(res_pressure)
plot(res_energy)
plot(res_sigma)
plot(res_sigma[30])
