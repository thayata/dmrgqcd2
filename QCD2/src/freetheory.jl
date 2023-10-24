function particle_momentum(k, L)
    p = (2k-1)π/(2L+1)/2
    return p
end

function particle_energy_free(p, w, M)
    Ep = sqrt( (2w*sin(p))^2 + M^2 )
    return Ep
end

function wavefunction_antiparticle_free(n, s, p, w, M, N)
    res = 0
    Ep = particle_energy_free(p, w, M)
    c = sqrt(4/(N+1))*sqrt((Ep+M)/(2Ep))
    if s == 1
        res = -2im*w*sin(p)/(Ep+M)*c*sin(2p*n)
    elseif s == 2
        res = c*cos(p*(2n-1))
    else
        res = false
    end
    return res
end

function wavefunction_particle_free(n, s, p, w, M, N)
    res = 0
    Ep = particle_energy_free(p, w, M)
    c = sqrt(4/(N+1))*sqrt((Ep+M)/(2Ep))
    if s == 1
        res = c*sin(2p*n)
    elseif s == 2
        res = -2im*w*sin(p)/(Ep+M)*c*cos(p*(2n-1))
    else
        res = false
    end
    return res
end



function vacuum_energy_free(w, M, μ, N, Nc=1, Nf=1)
    E0 = 0
    L = Int(N/2)
    for k in 1:L
        p = particle_momentum(k,L)
        Ep = particle_energy_free(p,w,M)
        E0 -= Ep
        if Ep <= μ
            E0 += Ep
        end
    end
    E0 *= Nc*Nf
    return E0
end

function quark_number_free_total(w, M, μ, N, Nc=1, Nf=1)

    Nb = 0
    L = Int(N/2)
    for k in 1:L
        p = particle_momentum(k, L)
        Ep = particle_energy_free(p,w,M)
        if Ep <= μ
            Nb += 1
        end
    end
    Nb *= Nc*Nf
    return Nb
end


function pressure_free(w, M, μ, N, Nc=1, Nf=1)
    
    pressure = 0
    L = Int(N/2)
    for k in 1:L
        p = particle_momentum(k, L)
        Ep = particle_energy_free(p,w,M)
        pressure += Ep
        if Ep <= μ
            pressure += μ-Ep
        end
    end
    pressure *= Nc*Nf
    return pressure
end

function quark_number_free(w, M, μ, N)
    # \psi^\dag\psi
    L = Int(N/2)
    res = [sum( real(propagator_lessor_free(n, s, n, s, w, M, μ, N))-1/2 for s in 1:2) for n = 1:L]
    return res
end


function quark_number_free_semianalytic_local(n, w, M, μ, N)
    # \psi^\dag\psi
    L = Int(N/2)
    res = 0
    c = 4/(N+1)
    kf =0
    for k in 1:L
        p = particle_momentum(k, L)
        Ep = particle_energy_free(p,w,M)
        if Ep < μ
            kf += +1
            #res += (1+ sin(p)*sin(p*(4n-1)))/2
            res += -M*(cos(p)*cos(p*(4n-1)))/(2Ep)
        end
        res += M*(cos(p)*cos(p*(4n-1)))/(2Ep)
    end
    res += kf/2
    res += -sin((4n*(kf)*pi)/(2L+1))/sin((2n*pi)/(2L+1))/8
    res += sin(((4n-2)*(kf)*pi)/(2L+1))/sin(((2n-1)*pi)/(2L+1))/8
    res = res*c
    return res
end


function quark_number_free_semianalytic(w, M, μ, N)
    L = Int(N/2)
    res = [quark_number_free_semianalytic_local(n, w, M, μ, N) for n = 1:L]
    return res
end


function chiral_condensate_free_semianalytic_local(n, w, M, μ, N)
    # \psi^\dag\psi
    L = Int(N/2)
    res = 0
    c = 4/(N+1)
    kf =0
    for k in 1:L
        p = particle_momentum(k, L)
        Ep = particle_energy_free(p,w,M)
        if Ep < μ
            kf += +1
            res += M*(1+sin(p)*sin(p*(4n-1)))/(2Ep)
        end
        res += -M*(1+sin(p)*sin(p*(4n-1)))/(2Ep)
    end
    res += -sin((4n*(kf)*pi)/(2L+1))/sin((2n*pi)/(2L+1))/8
    res += -sin(((4n-2)*(kf)*pi)/(2L+1))/sin(((2n-1)*pi)/(2L+1))/8
    res = res*c
    return res
end


function chiral_condensate_free_semianalytic(w, M, μ, N)
    L = Int(N/2)
    res = [chiral_condensate_free_semianalytic_local(n, w, M, μ, N) for n = 1:L]
    return res
end


function chiral_condensate_free(w, M, μ, N)
    # \bar{\psi}\psi
    L = Int(N/2)
    res =  [sum(real(propagator_lessor_free(n, s, n, s, w, M, μ, N))*(3-2s) for s in 1:2) for n = 1:L]
    return res
end

function vector_condensate_free(w, M, μ, N) 
    # \bar{\psi}\gamma^1\psi
    L = Int(N/2)
    res = [real(propagator_lessor_free(n, 1, n, 2, w, M, μ, N) + propagator_lessor_free(n, 2, n, 1, w, M, μ, N))  for n = 1:L]
    return res
end

function pseudoscalar_condensate_free(w, M, μ, N)
    # \bar{\psi}i\gamma_5\psi
    L = Int(N/2)
    res = [real(-im*propagator_lessor_free(n, 1, n, 2, w, M, μ, N) + im*propagator_lessor_free(n, 2, n, 1, w, M, μ, N))  for n = 1:L]
    return res
end


function propagator_greater_free(n1, s1, n2, s2, w, M, μ, N)
    # <\psi_{s1}(n1)\psi_{s2}^\dag(n2)>
    res = 0
    L = Int(N/2)
    for k in 1:L
        p = particle_momentum(k, L)
        Ep = particle_energy_free(p,w,M)
        if Ep >= μ
            u = wavefunction_particle_free(n1, s1, p, w, M, N)
            udag = conj(wavefunction_particle_free(n2, s2, p, w, M, N))
            res += u*udag
        end
    end
    return res

end

function propagator_lessor_free(n1, s1, n2, s2, w, M, μ, N)
    # <\psi_{s2}^\dag(n2)\psi_{s1}(n1)>
    res = -propagator_greater_free(n1, s1, n2, s2, w, M, μ, N)
    if s1 == s2 && n1 == n2
        res += 1
    end
    return res
end


function quark_distribution_free_old(w, M, μ, N)
    res = 0
    L = Int(N/2)
    X = Int(ceil(L/2))
    pseq = 2π/(2X)*(-X:X-1)
    yseq = -Int(X)+1:Int(X)-1
    n1 = X.+yseq
    n2 = X.-yseq
    prop_free = propagator_greater_free.(n1, 1, n2, 1, w, M, μ, N)+ propagator_greater_free.(n1, 2, n2, 2, w, M, μ, N)
    fourier_kernel = [exp(im*p*y) for p in pseq, y in yseq ]
    res = 2 .- 2real(fourier_kernel*prop_free)
    return pseq./2, res
end


function quark_distribution_base_free(w, M, μ, N; cutoff=0)
    L = Int(N/2) - cutoff
    yseq = -L+1:L-1
    n1seq = Int.(floor.((yseq.+L.+1)./2))
    n2seq = Int.(floor.((L.-yseq.+1)./2))
    prop_free = propagator_greater_free.(n1seq, 1, n2seq, 1, w, M, μ, N)+ propagator_greater_free.(n1seq, 2, n2seq, 2, w, M, μ, N)
    prop_free = prop_free./(exp.(1*(yseq.^2 .-400) ) .+1)
    return yseq, prop_free
end


function quark_distribution_free(w, M, μ, N; cutoff=0)
    res = 0
    L = Int(N/2)    
#    L = Int(N/2) - cutoff
#    yseq = -L+1:L-1
#    pseq = (2π/(2L-1)).*yseq
#    n1seq = Int.(floor.((yseq.+L.+1)./2))
#    n2seq = Int.(floor.((L.-yseq.+1)./2))
#    prop_free = propagator_greater_free.(n1seq, 1, n2seq, 1, w, M, μ, N)+ propagator_greater_free.(n1seq, 2, n2seq, 2, w, M, μ, N)
    yseq, prop_free = quark_distribution_base_free(w, M, μ, N; cutoff)
    pseq = (2π/(2L-1)).*yseq
fourier_kernel = [exp(im*p*y) for p in pseq, y in yseq ]
    res = 1 .- real(fourier_kernel*prop_free)
    return pseq, res
end



function fermi_momentum(w, M, μ, N)
    pf = 0
    Ef = 0
    L = Int(N/2)
    for k in 1:L
        p = particle_momentum(k, L)
        Ep = particle_energy_free(p, w, M)
        if Ef < Ep < μ
            Ef = Ep
            pf = p
        end
    end
    return pf
end

function baryon_propagator_lessor_free(n1, s1vec, n2, s2vec, w, M, μ, N)
    res = reduce(*,[propagator_lessor_free(n1, s1, n2, s2, w, M, μ, N) for (s1, s2) in zip(s1vec, s2vec)])
    return res
end

function baryon_distribution_free(w, M, μ, N, Nc)
    res = 0
    L = Int(N/2)
    yseq = -L+1:L-1
    pseq = (2π/(2L-1)).*yseq
    n1seq = Int.(floor.((yseq.+L.+1)./2))
    n2seq = Int.(floor.((L.-yseq.+1)./2))
    prop_free = zeros(length(yseq))
    for s in 1:2
        svec = fill(s,Nc)
        prop_free .+= [ baryon_propagator_lessor_free(n1, svec, n2, svec, w, M, μ, N) for (n1, n2) in zip(n1seq, n2seq)] 
    end
    fourier_kernel = [exp(im*p*y) for p in pseq, y in yseq ]
    res = real(fourier_kernel*prop_free)
    return pseq, res
end


###continuum 


function fermi_momentum_continuum(μ, m)
    res = 0
    if μ > m
        res= sqrt(μ^2-m^2)
    end
    return res
end

function pressure_free_continuum(μ, m)
    res = 0
    if μ > m
        pf = fermi_momentum_continuum(μ, m)
        res= μ*pf/(2π)-m^2/(4π) * log((μ+pf)/(μ-pf))
    end
    return res
end

function number_density_free_continuum(μ, m)
    res = 0
    if μ > m
        pf = fermi_momentum_continuum(μ, m)
        res= pf/π
    end
    return res
end
function energy_density_free_continuum(μ, m)
    res = 0
    if μ > m
        pf = fermi_momentum_continuum(μ, m)
        res= μ*pf/(2π)+m^2/(4π)*log((μ+pf)/(μ-pf))
    end
    return res
end

function cs_square_free_continuum(μ, m)
    res = 0
    if μ > m
        res = 1-m^2/μ^2
    end
    return res
end