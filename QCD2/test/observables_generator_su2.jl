using ITensors
using Plots
using ITensors.HDF5
include("../src/observables.jl")
include("../src/hamiltonian.jl")
include("../src/io.jl")
let

######################
# compute observables using a ground-state wavefunction generate by wave_function_generator_su2.jl #
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
g = 1.0   #coupling strength
J = 0.5*g*g #coeffcient of electric term

L = Int(N/2)
g=1.
J=g^2/2
Δμ = 0.01
μseq = 0.01:Δμ:2.0
cutoff = 0

target_dir = "../" #directory, which a h5 file of wavefunction is in
 
Nq_out = open("..//Nq_N=$N-w=$w-g=$g-M=$M-Lambda=$Λ-Nc=$Nc-Nf=$Nf.csv","w")
σ_out = open("..//sigma_N=$N-w=$w-g=$g-M=$M-Lambda=$Λ-Nc=$Nc-Nf=$Nf.csv","w")
fq_out= open("../fq_N=$N-w=$w-g=$g-M=$M-Lambda=$Λ-Nc=$Nc-Nf=$Nf-cutoff=$cutoff.csv","w")
thermo_out= open("../thermo_N=$N-w=$w-g=$g-M=$M-Lambda=$Λ-Nc=$Nc-Nf=$Nf.csv","w")

print(Nq_out,"μ")
print(σ_out,"μ")
for n in 1:L
    print(Nq_out,",$n")
    print(σ_out,",$n")
end
println(Nq_out,"")
println(σ_out,"")

print(fq_out,"μ")
for y in -L+cutoff+1:L-cutoff-1
    p = (2π/(2(L-cutoff)-1))*y
    print(fq_out,",$p")
end
println(fq_out,"")


println(thermo_out,"μ,pressure,E,Nq,σ")

for (n, μ) in enumerate(μseq)
    println("μ=$μ, N=$N, M=$M, Nc=$Nc, g=$g")
    sites, psi = load_wavefunction(target_dir, N, w, g, M, μ, Λ, Nc, Nf)
    H = hamiltonian(sites, w, J, M, μ, Λ, Nc, Nf)
    H0 = hamiltonian(sites, w, J, M, 0, Λ, Nc, Nf)
    energy = inner(psi',H0,psi)
    pressure = - inner(psi',H,psi)
    Nq = quark_number(psi, Nc, Nf)
    σ = chiral_condensate(psi, Nc, Nf)
    pseq, distribution = quark_distribution(psi, 1, 1, Nc, Nf, cutoff = cutoff)
    print(Nq_out,"$μ")
    for s in Nq
        print(Nq_out,",$s")
    end
    println(Nq_out,"")
    print(σ_out,"$μ")
    for s in σ
        print(σ_out,",$s")
    end
    println(σ_out,"")
    print(fq_out,"$μ")
    for s in distribution
        print(fq_out,",$s")
    end
    println(fq_out,"")
    totalNq = sum(Nq)
    totalσ = sum(σ)
    println(thermo_out,"$μ,$pressure,$energy,$totalNq,$totalσ")
end

close(thermo_out)
close(Nq_out)
close(σ_out)
close(fq_out)

end

