#generate wave_functions
using ITensors
using Plots
using ITensors.HDF5
include("../src/observables.jl")
include("../src/hamiltonian.jl")
include("../src/io.jl")
let

######################
# generate ground-state wavefunction using dmrg #
######################

    
######################
# Initial parameters #
######################
N = 40
M = 1.0 #quark mass
Λ = 400. #cutoff
w = 1.   #lattice scale 1/(2a)
Nc = 2   #color
Nf = 1   #flavor
g = 1.   #coupling strength
J = 0.5*g*g #coeffcient of electric term
Δμ = 0.01
μseq = 0.01:Δμ:2.0 #quark chemical potential

target_dir = "../" #directory to save a h5 file in

sites = siteinds_qcd2(N, Nc, Nf, false)
sweeps0 = Sweeps(100) # number of sweeps
maxdim!(sweeps0,20,20,40,40,60,60,80,80,100,100,200,200) # gradually increase states kept
cutoff!(sweeps0,1E-8) # desired truncation error

resweeps = Sweeps(100) # number of sweeps
maxdim!(resweeps,100) # gradually increase states kept
cutoff!(resweeps,1E-8) # desired truncation error

sweeps = sweeps0

psi = randomMPS(sites,10)

# mu dependence
for μ in μseq
    H =  hamiltonian(sites, w, J, M, μ, Λ, Nc, Nf)
    energy, psi = dmrg(H,psi,sweeps)
    save_wavefunction(psi, target_dir, N, w, g, M, μ, Λ, Nc, Nf)
    println("E=$energy, μ=$μ")
    sweeps = resweeps
end

end
