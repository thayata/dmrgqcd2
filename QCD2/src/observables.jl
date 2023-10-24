function mapping_coordinate(n, c, s, f, Nc, Nf, N, order="fnc")
    #mapping a physical site with position n, color c, spinor s, and flavor f
    # to the internal coordinate i
    if order == "fnc"
        i =  c + Nc*(2n-s) + N*Nc*(f-1)
    elseif order == "nfc"
        i = c + Nc*(f-1) + Nc*Nf*(2n-s)
    else
        i = -1
    end
    return i
end

function chiral_condensate_local(sz, n, Nc, Nf, N, order="fnc")
    # \bar{\psi}\psi
    res = 0
    for f in 1:Nf, c in 1:Nc
        i1 = mapping_coordinate(n, c, 1, f, Nc, Nf, N, order)
        i2 = mapping_coordinate(n, c, 2, f, Nc, Nf, N, order)
        res += sz[i1] - sz[i2]
    end
    return res
end

function quark_number_local(sz, n, Nc, Nf, N, order="fnc")
    # \psi^\dag\psi
    res = 0
    for f in 1:Nf, c in 1:Nc, s in 1:2
        i = mapping_coordinate(n, c, s, f, Nc, Nf, N, order)
        res += sz[i]
    end
    return res
end

function flavor_number_local(sz, n, f, Nc, Nf, N, order="fnc")
    # \psi^\dag_f \psi_f
    res = 0
    for c in 1:Nc, s in 1:2
        i = mapping_coordinate(n, c, s, f, Nc, Nf, N, order)
        res += sz[i]
    end
    return res
end

function quark_number(psi, Nc, Nf, order="fnc")
    # \psi^\dag\psi
    sz = expect(psi,"Sz")
    L = Int(length(psi)/(2*Nc*Nf))
    N = 2L
    res = [quark_number_local(sz, n, Nc, Nf, N, order) for n in 1:L]
    return res
end

function chiral_condensate(psi, Nc, Nf, order="fnc")
    # \bar{\psi}\psi
    sz = expect(psi,"Sz")
    L = Int(length(psi)/(2*Nc*Nf))
    N = 2L
    res = [chiral_condensate_local(sz, n, Nc, Nf, N, order) for n in 1:L]
    return res
end

function vector_condensate(psi, Nc, Nf, order="fnc")  #check it carefully
    # \bar{\psi}\gamma^1\psi
    L = Int(length(psi)/(2*Nc*Nf))
    res = [sum(fermion_propagator_lessor(psi, n, 1, f, n, 2, f, Nc, Nf) + fermion_propagator_lessor(psi, n, 2, f, n, 1, f, Nc, Nf) for f in 1:Nf) for n = 1:L]
    return res
end

function pseudoscalar_condensate(psi, Nc, Nf, order="fnc") #check it carefully
    # \bar{\psi}i\gamma_5\psi
    L = Int(length(psi)/(2*Nc*Nf))
    res = [sum(-im*fermion_propagator_lessor(psi, n, 1, f, n, 2, f, Nc, Nf,order) + im*fermion_propagator_lessor(psi, n, 2, f, n, 1, f, Nc, Nf,order) for f in 1:Nf) for n = 1:L]
    return res
end

function flavor_number(f, psi, Nc, Nf, order="fnc")
    # \psi^\dag_f \psi_f
    sz = expect(psi,"Sz")
    L = Int(length(psi)/(2*Nc*Nf))
    N = 2L
    res = [flavor_number_local(sz, n, f, Nc, Nf, N, order) for n in 1:L]
    return res
end

function fermion_propagator_internalcoordinate(psi, i, j)
    sites = siteinds(psi)
    if i<j
        A_i = op(sites,"S-",i)
        Adag_j = op(sites,"S+",j)

        orthogonalize!(psi, i)
        psidag=prime(dag(psi))

        if i==1
            Cij = psi[i]*A_i*psidag[i]
        else
            ##index linking i to i-1:
            lind = commonind(psi[i], psi[i-1])
            Cij = prime(psi[i],lind)*A_i*psidag[i]
        end

        for  k= i+1:j-1
            Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*psidag[k]
        end

        if j==length(psi)
            Cij *= psi[j]*Adag_j*psidag[j]
        else
            ##index linking j to j+1:
            rind = commonind(psi[j], psi[j+1])
            Cij *= prime(psi[j],rind)*Adag_j*psidag[j]
        end

        return -scalar(Cij)

    elseif j<i
        A_i = op(sites,"S-",i)
        Adag_j = op(sites,"S+",j)

        orthogonalize!(psi, j)
        psidag=prime(dag(psi))

        if j==1
            Cij = psi[j]*Adag_j*psidag[j]
        else
            ##index linking j to j-1:
            lind = commonind(psi[j], psi[j-1])
            Cij = prime(psi[j],lind)*Adag_j*psidag[j]
        end

        for  k= j+1:i-1
            Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*psidag[k]
        end

        if i==length(psi)
            Cij *= psi[i]*A_i*psidag[i]
        else
            ##index linking i to i+1:
            rind = commonind(psi[i], psi[i+1])
            Cij *= prime(psi[i],rind)*A_i*psidag[i]
        end
        return -scalar(Cij)

    else
        orthogonalize!(psi, i)
        return 0.5-scalar(psi[i]*op(sites,"Sz",i)*dag(prime(psi[i],"Site")))
    end
end

function fermion_propagator_internalcoordinate_dagger(psi, i, j)
    sites = siteinds(psi)
    if i<j
        Adag_i = op(sites,"S+",i)
        A_j = op(sites,"S-",j)

        orthogonalize!(psi, i)
        psidag=prime(dag(psi))

        if i==1
            Cij = psi[i]*Adag_i*psidag[i]
        else
            ##index linking i to i-1:
            lind = commonind(psi[i], psi[i-1])
            Cij = prime(psi[i],lind)*Adag_i*psidag[i]
        end

        for  k= i+1:j-1
            Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*psidag[k]
        end

        if j==length(psi)
            Cij *= psi[j]*A_j*psidag[j]
        else
            ##index linking j to j+1:
            rind = commonind(psi[j], psi[j+1])
            Cij *= prime(psi[j],rind)*A_j*psidag[j]
        end

        return scalar(Cij)

    elseif j<i
        Adag_i = op(sites,"S-",i)
        A_j = op(sites,"S+",j)

        orthogonalize!(psi, j)
        psidag=prime(dag(psi))

        if j==1
            Cij = psi[j]*A_j*psidag[j]
        else
            ##index linking j to j-1:
            lind = commonind(psi[j], psi[j-1])
            Cij = prime(psi[j],lind)*A_j*psidag[j]
        end

        for  k= j+1:i-1
            Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*psidag[k]
        end
        if i==length(psi)
            Cij *= psi[i]*Adag_i*psidag[i]
        else
            ##index linking i to i+1:
            rind = commonind(psi[i], psi[i+1])
            Cij *= prime(psi[i],rind)*Adag_i*psidag[i]
        end
        return scalar(Cij)

    else
        orthogonalize!(psi, i)
        return 0.5+scalar(psi[i]*op(sites,"Sz",i)*dag(prime(psi[i],"Site")))
    end
end

function fermion_propagator_greater(psi, n1, s1, f1, n2, s2, f2, Nc, Nf, order="fnc")
    # n1, n2 physical site
    # s1 s2 spinor
    # f1, f2 flavor
    L = Int(length(psi)/(2*Nc*Nf))
    N = 2L
    res = 0
    for c in 1:Nc
        i = mapping_coordinate(n1, c, s1, f1, Nc, Nf, N, order)
        j = mapping_coordinate(n2, c, s2, f2, Nc, Nf, N, order)
        res += fermion_propagator_internalcoordinate(psi, i, j)
    end
    return res*(-1)^(n1+n2)*(-im)^(s1-1)*(im)^(s2-1)
end

function fermion_propagator_lessor(psi, n1, s1, f1, n2, s2, f2, Nc, Nf, order="fnc")
    # n1, n2 physical site
    # s1 s2 spinor
    # f1, f2 flavor
    L = Int(length(psi)/(2*Nc*Nf))
    N = 2L
    res = 0
    for c in 1:Nc
        i = mapping_coordinate(n1, c, s1, f1, Nc, Nf, N, order)
        j = mapping_coordinate(n2, c, s2, f2, Nc, Nf, N, order)
        res += fermion_propagator_internalcoordinate_dagger(psi, j, i) # note ordering is j, i
    end
    return res*(-1)^(n1+n2)*(-im)^(s1-1)*(im)^(s2-1)
end


function quark_distribution_old(psi, f1, f2, Nc, Nf, order="fnc")
    res = 0
    L = Int(length(psi)/(2*Nc*Nf))
    N = 2L
    X = Int(ceil(L/2))
    pseq = 2π/(2X)*(-X:X-1)
    yseq = -Int(X)+1:Int(X)-1
    fourier_kernel = [exp(im*p*y) for p in pseq, y in yseq]
    prop = [sum(fermion_propagator_greater(psi, X+y, s, f1, X-y, s, f2, Nc, Nf, "fnc") for s in 1:2) for y in yseq]
    res = 2Nc .- 2real(fourier_kernel*prop)
    return pseq./2, res
end

function quark_distribution(psi, f1, f2, Nc, Nf, order="fnc";cutoff = 0)
    res = 0
    L = Int(length(psi)/(2*Nc*Nf))-cutoff
    yseq = -L+1:L-1
    prop = [sum(fermion_propagator_greater(psi, Int(floor((y+L+1)/2)), s, f1, Int(floor((L-y+1)/2)), s, f2, Nc, Nf, "fnc") for s in 1:2) for y in yseq]
    pseq = (2π/(2L-1)).*yseq
    fourier_kernel = [exp(im*p*y) for p in pseq, y in yseq ]
    res = Nc .- real(fourier_kernel*prop)
    return pseq, res
end

function fourier_transform(fx)
    L = length(fx)
    xseq = 1:L
    pseq = 2π/L .*xseq
    fourier_kernel = [exp(im*p*x) for p in pseq, x in xseq ]
    fp = fourier_kernel*fx
    return pseq, fp
end

function density_correlation_internalcoordinate(psi, op1, op2, i, j)

    sites = siteinds(psi)
    
    if i<j

        orthogonalize!(psi, i)
        rind = commonind(psi[i], psi[i+1])
        Cij=psi[i]*op(op1,sites,i)*dag(prime(prime(psi[i],"Site"),rind))

        for  k= i+1:j-1
            Cij *= psi[k]*dag(prime(psi[k],"Link"))
        end

        lind = commonind(psi[j], psi[j-1])
        Cij *= psi[j]*op(op2,sites,j)*dag(prime(prime(psi[j],"Site"),lind))

        return scalar(Cij)

    elseif j<i

        orthogonalize!(psi, j)
        rind = commonind(psi[j], psi[j+1])
        Cij=psi[j]*op(op2,sites,j)*dag(prime(prime(psi[j],"Site"),rind))

        for  k= j+1:i-1
            Cij *= psi[k]*dag(prime(psi[k],"Link"))
        end

        lind = commonind(psi[i], psi[i-1])
        Cij *= psi[i]*op(op1,sites,i)*dag(prime(prime(psi[i],"Site"),lind))

        return scalar(Cij)

    else
        orthogonalize!(psi, i)
        return scalar(psi[i]*op(op1*op2,sites,i)*dag(prime(psi[i],"Site")))
    end
end

function density_correlation(psi, op1, op2, n1, n2, Nc, Nf, order="fnc")
    # n1, n2 physical site
    # s1 s2 spinor
    # f1, f2 flavor
    L = Int(length(psi)/(2*Nc*Nf))
    N = 2L
    res = 0

    for s1 in 1:2, f1 in 1:Nf, c1 in 1:Nc
        for s2 in 1:2, f2 in 1:Nf, c2 in 1:Nc
            Op1=(op1=="density")*[1 0; 0 0]-(-1)^s1*(op1=="chiral")*[1/2 0; 0 -1/2]
            Op2=(op2=="density")*[1 0; 0 0]-(-1)^s2*(op2=="chiral")*[1/2 0; 0 -1/2]
            i = mapping_coordinate(n1, c1, s1, f1, Nc, Nf, N, order)
            j = mapping_coordinate(n2, c2, s2, f2, Nc, Nf, N, order)
            res += density_correlation_internalcoordinate(psi, Op1, Op2, i, j) 
        end
    end

    return res
end

function baryon_propagator_internalcoordinate(psi, x, s, Nc)
    sites = siteinds(psi)
    psidag = prime(dag(psi))
    Cij=1.0

    i=1
    orthogonalize!(psi, x[i])

    if x[i]==1
       if x[i]==x[i+1]
            Cij*=psi[x[i]]*(op(sites,"Sz",x[i])+0.5*op(sites,"Id",x[i]))*psidag[x[i]]
       else
            A_i = -s[i]*op(sites,"S-",x[i])+(1-s[i])*op(sites,"S+",x[i])
            Cij*= psi[x[i]]*A_i*psidag[x[i]]
       end
    else
       lind = commonind(psi[x[i]], psi[x[i]-1])
       if x[i]==x[i+1]
            Cij*=prime(psi[x[i]],lind)*(op(sites,"Sz",x[i])+0.5*op(sites,"Id",x[i]))*psidag[x[i]]
       else
            A_i = -s[i]*op(sites,"S-",x[i])+(1-s[i])*op(sites,"S+",x[i])
            Cij  *= prime(psi[x[i]],lind)*A_i*psidag[x[i]]
       end
    end

    for  k= x[i]+1:x[i+1]-1
         Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*psidag[k]
    end

    for i=2:2Nc-2

        if x[i]!=x[i-1]
            if x[i]==x[i+1]
                Cij*=(-1)^(i%2==0)*psi[x[i]]*(op(sites,"Sz",x[i])+0.5*op(sites,"Id",x[i]))*psidag[x[i]]
            else
                A_i = (-1)^(i%2==1)*s[i]*op(sites,"S-",x[i])+(1-s[i])*op(sites,"S+",x[i])
                Cij*= psi[x[i]]*A_i*psidag[x[i]]
            end
        end

        if i%2==1
            for  k= x[i]+1:x[i+1]-1
                Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*psidag[k]
            end
        else
            for  k= x[i]+1:x[i+1]-1
                Cij *= psi[k]*noprime(psidag[k],"Site")
            end
        end
    end

    i=2Nc-1

    if x[i]!=x[i-1]
       if x[i]==x[i+1]
          if x[i]==length(psi)
                Cij*=psi[x[i]]*(op(sites,"Sz",x[i])+0.5*op(sites,"Id",x[i]))*psidag[x[i]]
          else
                rind = commonind(psi[x[i]], psi[x[i]+1])
                Cij*=prime(psi[x[i]],rind)*(op(sites,"Sz",x[i])+0.5*op(sites,"Id",x[i]))*psidag[x[i]]
          end
       else
          A_i = -s[i]*op(sites,"S-",x[i])+(1-s[i])*op(sites,"S+",x[i])
          Cij*= psi[x[i]]*A_i*psidag[x[i]]

          for  k= x[i]+1:x[i+1]-1
            Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*psidag[k]
          end

          A_2Nc = s[2Nc]*op(sites,"S-",x[2Nc])+(1-s[2Nc])*op(sites,"S+",x[2Nc])

          if x[2Nc]==length(psi)
             Cij *= psi[x[2Nc]]*A_2Nc*psidag[x[2Nc]]
          else
            rind = commonind(psi[x[2Nc]], psi[x[2Nc]+1])
            Cij *= prime(psi[x[2Nc]],rind)*A_2Nc*psidag[x[2Nc]]
          end
       end
    else
        for  k= x[i]+1:x[i+1]-1
             Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*psidag[k]
        end

        A_2Nc = s[2Nc]*op(sites,"S-",x[2Nc])+(1-s[2Nc])*op(sites,"S+",x[2Nc])

        if x[2Nc]==length(psi)
             Cij *= psi[x[2Nc]]*A_2Nc*psidag[x[2Nc]]
        else
            rind = commonind(psi[x[2Nc]], psi[x[2Nc]+1])
            Cij *= prime(psi[x[2Nc]],rind)*A_2Nc*psidag[x[2Nc]]
        end
    end

   return scalar(Cij)
end

function baryon_propagator_internalcoordinate_debug(psi, x, s, Nc)
    ##here  (-1) phases from sorting in baryon_propagator_lessor_coodinate_mapping is not taken into account 
    ## (-1) phases in jordan wigner transformation is taken into account
    sites = siteinds(psi)
    Cij=1.0

    i=1
    orthogonalize!(psi, x[i])
    rind = commonind(psi[x[i]], psi[x[i]+1])

    if x[i]==x[i+1]
       Cij*=psi[x[i]]*(op(sites,"Sz",x[i])+0.5*op(sites,"Id",x[i]))*dag(prime(prime(psi[x[i]],"Site"),rind))
    else
       A_i = -s[i]*op(sites,"S-",x[i])+(1-s[i])*op(sites,"S+",x[i])
       Cij  *= psi[x[i]]*A_i*dag(prime(prime(psi[x[i]],"Site"),rind))
    end

    for  k= x[i]+1:x[i+1]-1
         Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*dag(prime(prime(psi[k],"Site"),"Link"))
    end

    for i=2:2Nc-2

        if x[i]!=x[i-1]
            if x[i]==x[i+1]
                Cij*=(-1)^(i%2==0)*psi[x[i]]*(op(sites,"Sz",x[i])+0.5*op(sites,"Id",x[i]))*dag(prime(prime(psi[x[i]],"Site"),"Link"))
            else
                A_i = (-1)^(i%2==1)*s[i]*op(sites,"S-",x[i])+(1-s[i])*op(sites,"S+",x[i])
                Cij*= psi[x[i]]*A_i*dag(prime(prime(psi[x[i]],"Site"),"Link"))
            end
        end

        if i%2==1
            for  k= x[i]+1:x[i+1]-1
                Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*dag(prime(prime(psi[k],"Site"),"Link"))
            end
        else
            for  k= x[i]+1:x[i+1]-1
                Cij *= psi[k]*dag(prime(psi[k],"Link"))
            end
        end
    end

    i=2Nc-1

    if x[i]!=x[i-1]
       if x[i]==x[i+1]
          lind = commonind(psi[x[i]], psi[x[i]-1])
          Cij*=psi[x[i]]*(op(sites,"Sz",x[i])+0.5*op(sites,"Id",x[i]))*dag(prime(prime(psi[x[i]],"Site"),lind))
       else
          A_i = -s[i]*op(sites,"S-",x[i])+(1-s[i])*op(sites,"S+",x[i])
          Cij*= psi[x[i]]*A_i*dag(prime(prime(psi[x[i]],"Site"),"Link"))

          for  k= x[i]+1:x[i+1]-1
            Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*dag(prime(prime(psi[k],"Site"),"Link"))
          end

          A_2Nc = s[2Nc]*op(sites,"S-",x[2Nc])+(1-s[2Nc])*op(sites,"S+",x[2Nc])
          lind = commonind(psi[x[2Nc]], psi[x[2Nc]-1])
          Cij *= psi[x[2Nc]]*A_2Nc*dag(prime(prime(psi[x[2Nc]],"Site"),lind))
       end
    else
        for  k= x[i]+1:x[i+1]-1
             Cij *= psi[k]*(-2.0)*op(sites,"Sz",k)*dag(prime(psi[k],"Link"))
        end

        A_2Nc = s[2Nc]*op(sites,"S-",x[2Nc])+(1-s[2Nc])*op(sites,"S+",x[2Nc])
        lind = commonind(psi[x[2Nc]], psi[x[2Nc]-1])
        Cij *= psi[x[2Nc]]*A_2Nc*dag(prime(prime(psi[x[2Nc]],"Site"),lind))                            
    end

   return scalar(Cij)
end

function baryon_coordinate_mapping(n, svec, fvec, Nc, Nf, N, order="fnc"; flg=1)
    if flg == 1 #psi
        cseq = 1:Nc
    elseif flg == 0 #psidag
        cseq = Nc:-1:1
    else
        return false
    end
    if length(svec) == Nc | length(fvec) == Nc
        baryon_internal_coordinate = mapping_coordinate.(n, cseq, svec, fvec, Nc, Nf, N, order)
    else
        println("length of s or f is not matched in baryon_coordinate_mapping")
        baryon_internal_coordinate =  false
    end
    return baryon_internal_coordinate
end

function baryon_propagator_lessor_coodinate_mapping(n1, s1, f1, n2, s2, f2, Nc, Nf, N, order="fnc")

    psidag_n = baryon_coordinate_mapping(n1, s1, f1, Nc, Nf, N, order)
    psi_n = baryon_coordinate_mapping(n2, s2, f2, Nc, Nf, N, order)
    lst = [[i, 1] for i in psi_n]
    append!(lst, [[i, 0] for i in psidag_n])
    xvec = []
    svec = []
        
    for zvec in sort(lst)
            push!(xvec,zvec[1])
            push!(svec,zvec[2])
    end

    return xvec, svec
end

function baryon_propagator_lessor(psi, n1, s1, f1, n2,s2, f2, Nc, Nf, order="fnc")
    N = Int(length(psi)/(Nc*Nf))
    xvec, svec = baryon_propagator_lessor_coodinate_mapping(n1, s1, f1, n2, s2, f2, Nc, Nf, N, order)
    res = baryon_propagator_internalcoordinate_debug(psi, xvec, svec, Nc)
    res = res*(-1)^(Nc*(n1-n2))*(im)^(sum(s2)-Nc)*(-im)^(sum(s1)-Nc)
    return res
end

function diquark_propagator_greator(psi, n1, n2, s, Nc, Nf, order="fnc")

    if s==1
        N = Int(length(psi)/(Nc*Nf))
        s0vec= fill(1,Nc)
        fvec= fill(1,Nc)
        xvec, svec = baryon_propagator_lessor_coodinate_mapping(n1, s0vec, fvec, n2, s0vec, fvec, Nc, Nf, N, order)
        res = baryon_propagator_internalcoordinate_debug(psi, xvec, svec, Nc)
        res = res*(-1)^(n1!=n2)
        return res
    else
        s1vec= fill(2,Nc)
        xvec, svec = baryon_propagator_lessor_coodinate_mapping(n1, s1vec, fvec, n2, s1vec, fvec, Nc, Nf, N, order)
        res2 = baryon_propagator_internalcoordinate_debug(psi, xvec, svec, Nc)
        res2 = res2*(-1)^(n1!=n2)
        return res2
    end
end

function apply_baryon_operator(psi, n, svec, fvec, Nc, Nf, order="fnc")
    N = Int(length(psi)/(Nc*Nf))
    psidag_n = baryon_coordinate_mapping(n, svec, fvec, Nc, Nf, N, order, flg=0)
    for i in psidag_n
        Adag_j = op(sites,"S+",i)
        newpsi = Adag_j*psi[i]
        noprime!(newpsi)
        psi[i] = newpsi
    end
    psi = psi*(-1)^(Nc*n)*(-im)^(sum(svec)-Nc)
end
