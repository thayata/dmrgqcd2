function save_wavefunction(psi, target_dir, N, w, g, M, μ, Λ, Nc, Nf; overwrite=false)
    if length(psi) != Nc*Nf*N
        println("size of wave function mismatch")
        return false
    end
    file_base = "psi_N=$N-w=$w-g=$g-M=$M-mu=$μ-Lambda=$Λ-Nc=$Nc-Nf=$Nf.h5"
    file_name = joinpath(target_dir, file_base)
    if (isfile(file_name) == false) | (overwrite == true)
        f = h5open(file_name,"w")
        write(f,"psi",psi)
        close(f)
        return true
    else
        println("file exists.")
        return false
    end
end

function load_wavefunction(target_dir, N, w, g, M, μ, Λ, Nc, Nf)
    psi = false
    file_base = "psi_N=$N-w=$w-g=$g-M=$M-mu=$μ-Lambda=$Λ-Nc=$Nc-Nf=$Nf.h5"
    file_name = joinpath(target_dir, file_base)
    try
        f = h5open(file_name,"r")
        psi = read(f,"psi",MPS)
        close(f)
    catch
        files = readdir(target_dir)
        for fname in files
            m = match(r"psi_N=([0-9]+)-w=([0-9]*\.*[0-9]*)-g=([0-9]*\.*[0-9]*)-M=([0-9]*\.*[0-9]*)-mu=([0-9]*\.*[0-9]*)-Lambda=([0-9]*\.*[0-9]*)-Nc=([0-9]+)-Nf=([0-9]+)\.h5", fname)
            if m !== nothing
                 (NN, ww, gg, MM, μμ, ΛΛ, NNc, NNf) = m.captures
                 if (parse(Int,NN)==N) && (parse(Float64,ww)==w) && (parse(Float64,gg)==g) && (parse(Float64,MM)==M) && (parse(Float64,μμ)==μ) && (parse(Float64,ΛΛ)==Λ) && (parse(Int,NNc)==Nc) && (parse(Int,NNf)==Nf) 
                    file_name = joinpath(target_dir, fname)
                    println(file_base)
                    f = h5open(file_name,"r")
                    psi = read(f,"psi",MPS)
                    close(f)
                 end
            else
                println("File not found:N=$N, w=$w, g=$g, M=$M, μ=$μ, Λ=$Λ, Nc=$Nc, Nf=$Nf")
            end
        end
    end
    sites = siteinds(psi)
    return sites, psi
end
