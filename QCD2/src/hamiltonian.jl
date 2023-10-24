function siteinds_qcd2(N, Nc, Nf, conserve_qns)
    sites = siteinds("S=1/2", Nc*Nf*N;conserve_qns = conserve_qns)
    return sites
end

function hamiltonian(sites, w, J, M, μ, Lambda, Nc, Nf, order="fnc")
    if order == "fnc"
        H = hamiltonian_fnc(sites, w, J, M, μ, Lambda, Nc, Nf)
    elseif order == "nfc"
        H = hamiltonian_nfc(sites, w, J, M, μ, Lambda, Nc, Nf)
    end
end


## e.g., Nc=3 Nf=2 is labeld as
## |n_ur(1)n_ug(1)n_ub(1)n_dr(1)n_dg(1)n_db(1)...n_ur(N)n_ug(N)n_ub(N)n_dr(N)n_dg(N)n_db(N)>
function hamiltonian_nfc(sites, w, J, M, μ, Lambda, Nc, Nf)
  N = Int(length(sites)/Nc/Nf)
  ampo = OpSum()

  t=w*(-2.0)^(Nc*Nf-1)
  for i=1:N-1
      ii=Nf*Nc*(i-1)

      for j=1:(Nc*Nf)
          jj=ii+j
          list=[t,"S+",jj]
          for k=1:(Nc*Nf-1)
              push!(list,"Sz")
              push!(list,jj+k)
          end
          push!(list,"S-")
          push!(list,jj+Nc*Nf)
          ampo += tuple(list...)
      end

      for j=1:(Nc*Nf)
          jj=ii+j
          list=[t,"S-",jj]
          for k=1:(Nc*Nf-1)
              push!(list,"Sz")
              push!(list,jj+k)
          end
          push!(list,"S+")
          push!(list,jj+Nc*Nf)
          ampo += tuple(list...)
      end
  end

  for i=1:N
      ii=Nc*Nf*(i-1)
      for j=1:(Nc*Nf)
          ampo += M*(-1)^i-μ,"Sz",ii+j;
      end
  end

## QaQa same sites

  for i=1:N
      ii=Nc*Nf*(i-1)

      for j=1:Nf
          jj=ii+Nc*(j-1)
          for k=1:Nc-1
              for m=k+1:Nc
                  ampo += 0.5*J*(N-i)+0.5*Lambda,"S+",jj+k,"S-",jj+k,"S-",jj+m,"S+",jj+m
                  ampo += 0.5*J*(N-i)+0.5*Lambda,"S-",jj+k,"S+",jj+k,"S+",jj+m,"S-",jj+m
              end
          end
      end

      for j1=1:Nf-1
         jj1=ii+Nc*(j1-1)
         for j2=j1+1:Nf
          jj2=ii+Nc*(j2-1)
          for k=1:Nc-1
              for m=k+1:Nc
                  list=[4^(m-1-k)*J*(N-i)+4^(m-1-k)*Lambda,"S+",jj1+k]
                  for u=k+1:m-1
                      push!(list,"Sz")
                      push!(list,jj1+u)
                  end
                  push!(list,"S-")
                  push!(list,jj1+m)
                  push!(list,"S-")
                  push!(list,jj2+k)
                  for u=k+1:m-1
                       push!(list,"Sz")
                       push!(list,jj2+u)
                   end
                   push!(list,"S+")
                   push!(list,jj2+m)
                   ampo += tuple(list...)
                  list=[4^(m-1-k)*J*(N-i)+4^(m-1-k)*Lambda,"S-",jj1+k]
                  for u=k+1:m-1
                      push!(list,"Sz")
                      push!(list,jj1+u)
                  end
                  push!(list,"S+")
                  push!(list,jj1+m)
                  push!(list,"S+")
                  push!(list,jj2+k)
                  for u=k+1:m-1
                       push!(list,"Sz")
                       push!(list,jj2+u)
                   end
                   push!(list,"S-")
                   push!(list,jj2+m)
                   ampo += tuple(list...)
              end
            end
         end
      end

        
      for j=1:Nf
          jj=ii+Nc*(j-1)
          for k=2:Nc
              for m=1:k-1
                  for n=1:k-1
                      ampo += 0.5*J*(N-i)/(k*(k-1))+0.5*Lambda/(k*(k-1)),"Sz",jj+m,"Sz",jj+n
                      ampo += -0.5*J*(N-i)/(k*(k-1))-0.5*Lambda/(k*(k-1)),"Sz",jj+m,"Sz",jj+k
                      ampo += -0.5*J*(N-i)/(k*(k-1))-0.5*Lambda/(k*(k-1)),"Sz",jj+k,"Sz",jj+n
                      ampo += 0.125*J*(N-i)/(k*(k-1))+0.125*Lambda/(k*(k-1)),"Id",jj+k
                  end
              end
          end
      end

      for j1=1:Nf-1
          jj1=ii+Nc*(j1-1)
          for j2=j1+1:Nf
             jj2=ii+Nc*(j2-1)
             for k=2:Nc
              for m=1:k-1
                  for n=1:k-1
                      ampo += J*(N-i)/(k*(k-1))+Lambda/(k*(k-1)),"Sz",jj1+m,"Sz",jj2+n
                      ampo += -J*(N-i)/(k*(k-1))-Lambda/(k*(k-1)),"Sz",jj1+m,"Sz",jj2+k
                      ampo += -J*(N-i)/(k*(k-1))-Lambda/(k*(k-1)),"Sz",jj1+k,"Sz",jj2+n
                      ampo += J*(N-i)/(k*(k-1))+Lambda/(k*(k-1)),"Sz",jj1+k,"Sz",jj2+k
                  end
              end
          end
       end
      end  
  end

## QaQa diffrent sites

  for i=1:N-1
      for j=i+1:N
          ii=Nc*Nf*(i-1)
          jj=Nc*Nf*(j-1)

          for k1=1:Nf
              kk1=ii+Nc*(k1-1)
              for k2=1:Nf
                kk2=jj+Nc*(k2-1)              
                  for m=1:Nc-1
                      for n=m+1:Nc
                          list=[4^(n-1-m)*J*(N-j)+4^(n-1-m)*Lambda,"S+",kk1+m]
                          for u=m+1:n-1
                              push!(list,"Sz")
                              push!(list,kk1+u)
                          end
                          push!(list,"S-")
                          push!(list,kk1+n)
                          push!(list,"S-")
                          push!(list,kk2+m)
                          for u=m+1:n-1
                              push!(list,"Sz")
                              push!(list,kk2+u)
                          end
                          push!(list,"S+")
                          push!(list,kk2+n)
                          ampo += tuple(list...)
                          list=[4^(n-1-m)*J*(N-j)+4^(n-1-m)*Lambda,"S-",kk1+m]
                          for u=m+1:n-1
                              push!(list,"Sz")
                              push!(list,kk1+u)
                          end
                          push!(list,"S+")
                          push!(list,kk1+n)
                          push!(list,"S+")
                          push!(list,kk2+m)
                          for u=m+1:n-1
                              push!(list,"Sz")
                            push!(list,kk2+u)
                          end
                          push!(list,"S-")
                          push!(list,kk2+n)
                          ampo += tuple(list...)
                    end
                  end
                 end
           end
                
          for k1=1:Nf
              kk1=ii+Nc*(k1-1)
              for k2=1:Nf
                kk2=jj+Nc*(k2-1)
                for m=2:Nc
                  for u=1:m-1
                      for v=1:m-1
                          ampo += J*(N-j)/(m*(m-1))+Lambda/(m*(m-1)),"Sz",kk1+u,"Sz",kk2+v
                          ampo +=-J*(N-j)/(m*(m-1))-Lambda/(m*(m-1)),"Sz",kk1+u,"Sz",kk2+m
                          ampo +=-J*(N-j)/(m*(m-1))-Lambda/(m*(m-1)),"Sz",kk1+m,"Sz",kk2+v
                          ampo += J*(N-j)/(m*(m-1))+Lambda/(m*(m-1)),"Sz",kk1+m,"Sz",kk2+m
                      end
                  end
              end
             end
           end
     
     end
  end

  H = MPO(ampo,sites)
  return H
end

## e.g., Nc=3 Nf=2 is labeld as
## |n_ur(1)n_ug(1)n_ub(1)...n_ur(N)n_ug(N)n_ub(N);n_dr(1)n_dg(1)n_db(1)...n_dr(N)n_dg(N)n_db(N);s flavor...>
function hamiltonian_fnc(sites, w, J, M, μ, Lambda, Nc, Nf)
  N = Int(length(sites)/Nc/Nf)
  ampo = OpSum()
  
  t=w*(-2.0)^(Nc-1)
  for f=1:Nf
    for i=1:N-1
      ii=N*Nc*(f-1)+Nc*(i-1)
      for j=1:Nc
          jj=ii+j
          list=[t,"S+",jj]
          for k=1:Nc-1
              push!(list,"Sz")
              push!(list,jj+k)
          end
          push!(list,"S-")
          push!(list,jj+Nc)
          ampo += tuple(list...)
      end

      for j=1:Nc
          jj=ii+j
          list=[t,"S-",jj]
          for k=1:Nc-1
              push!(list,"Sz")
              push!(list,jj+k)
          end
          push!(list,"S+")
          push!(list,jj+Nc)
          ampo += tuple(list...)
      end
      end
    end
    
  for f=1:Nf
    for i=1:N
      for j=1:Nc
          ii=N*Nc*(f-1)+Nc*(i-1)+j
          ampo += M*(-1)^i-μ,"Sz",ii ;
      end
    end
  end

## QaQa same sites

 for f=1:Nf
  for i=1:N
      ii=N*Nc*(f-1)+Nc*(i-1)
      for k=1:Nc-1
          for m=k+1:Nc
              ampo += 0.5*J*(N-i)+0.5*Lambda,"S+",ii+k,"S-",ii+k,"S-",ii+m,"S+",ii+m
              ampo += 0.5*J*(N-i)+0.5*Lambda,"S-",ii+k,"S+",ii+k,"S+",ii+m,"S-",ii+m
          end
      end
        
      for k=2:Nc
          for m=1:k-1
              for n=1:k-1
                  ampo += 0.5*J*(N-i)/(k*(k-1))+0.5*Lambda/(k*(k-1)),"Sz",ii+m,"Sz",ii+n
                  ampo += -0.5*J*(N-i)/(k*(k-1))-0.5*Lambda/(k*(k-1)),"Sz",ii+m,"Sz",ii+k
                  ampo += -0.5*J*(N-i)/(k*(k-1))-0.5*Lambda/(k*(k-1)),"Sz",ii+k,"Sz",ii+n
                  ampo += 0.125*J*(N-i)/(k*(k-1))+0.125*Lambda/(k*(k-1)),"Id",ii+k
              end
           end
      end
   end
  end      

 for f1=1:Nf-1
     for f2=f1+1:Nf
         for i=1:N
             ii1=N*Nc*(f1-1)+Nc*(i-1)
             ii2=N*Nc*(f2-1)+Nc*(i-1)
              
             for k=1:Nc-1
                 for m=k+1:Nc
                    list=[4^(m-1-k)*J*(N-i)+4^(m-1-k)*Lambda,"S+",ii1+k]
                    for u=k+1:m-1
                        push!(list,"Sz")
                        push!(list,ii1+u)
                    end
                    push!(list,"S-")
                    push!(list,ii1+m)
                    push!(list,"S-")
                    push!(list,ii2+k)
                    for u=k+1:m-1
                        push!(list,"Sz")
                        push!(list,ii2+u)
                    end
                    push!(list,"S+")
                    push!(list,ii2+m)
                    ampo += tuple(list...)
                    list=[4^(m-1-k)*J*(N-i)+4^(m-1-k)*Lambda,"S-",ii1+k]
                    for u=k+1:m-1
                        push!(list,"Sz")
                        push!(list,ii1+u)
                    end
                    push!(list,"S+")
                    push!(list,ii1+m)
                    push!(list,"S+")
                    push!(list,ii2+k)
                    for u=k+1:m-1
                        push!(list,"Sz")
                        push!(list,ii2+u)
                    end
                    push!(list,"S-")
                    push!(list,ii2+m)
                    ampo += tuple(list...)
                end
            end

            for k=2:Nc
                for m=1:k-1
                    for n=1:k-1
                        ampo += J*(N-i)/(k*(k-1))+Lambda/(k*(k-1)),"Sz",ii1+m,"Sz",ii2+n
                        ampo += -J*(N-i)/(k*(k-1))-Lambda/(k*(k-1)),"Sz",ii1+m,"Sz",ii2+k
                        ampo += -J*(N-i)/(k*(k-1))-Lambda/(k*(k-1)),"Sz",ii1+k,"Sz",ii2+n
                        ampo += J*(N-i)/(k*(k-1))+Lambda/(k*(k-1)),"Sz",ii1+k,"Sz",ii2+k
                    end
                end
            end
      end
   end  
end
        
## QaQa diffrent sites

for f1=1:Nf
    for f2=1:Nf
        for i=1:N-1
            for j=i+1:N
                ii=N*Nc*(f1-1)+Nc*(i-1)
                jj=N*Nc*(f2-1)+Nc*(j-1)

                for m=1:Nc-1
                    for n=m+1:Nc
                        list=[4^(n-1-m)*J*(N-j)+4^(n-1-m)*Lambda,"S+",ii+m]
                        for u=m+1:n-1
                              push!(list,"Sz")
                              push!(list,ii+u)
                         end
                         push!(list,"S-")
                         push!(list,ii+n)
                         push!(list,"S-")
                         push!(list,jj+m)
                         for u=m+1:n-1
                             push!(list,"Sz")
                             push!(list,jj+u)
                         end
                         push!(list,"S+")
                         push!(list,jj+n)
                         ampo += tuple(list...)
                         list=[4^(n-1-m)*J*(N-j)+4^(n-1-m)*Lambda,"S-",ii+m]
                         for u=m+1:n-1
                             push!(list,"Sz")
                             push!(list,ii+u)
                         end
                         push!(list,"S+")
                         push!(list,ii+n)
                         push!(list,"S+")
                         push!(list,jj+m)
                         for u=m+1:n-1
                              push!(list,"Sz")
                            push!(list,jj+u)
                          end
                          push!(list,"S-")
                          push!(list,jj+n)
                          ampo += tuple(list...)
                    end
                  end
                
                  for m=2:Nc
                      for u=1:m-1
                          for v=1:m-1
                              ampo += J*(N-j)/(m*(m-1))+Lambda/(m*(m-1)),"Sz",ii+u,"Sz",jj+v
                              ampo +=-J*(N-j)/(m*(m-1))-Lambda/(m*(m-1)),"Sz",ii+u,"Sz",jj+m
                              ampo +=-J*(N-j)/(m*(m-1))-Lambda/(m*(m-1)),"Sz",ii+m,"Sz",jj+v
                              ampo += J*(N-j)/(m*(m-1))+Lambda/(m*(m-1)),"Sz",ii+m,"Sz",jj+m
                          end
                      end
                  end
             
          end
        end     
     end
  end

  H = MPO(ampo,sites)
  return H
end

function densityMPO(sites, op1, op2, n1, n2, Nc, Nf,order)

    ampo = OpSum()

    for s1 in 1:2, f1 in 1:Nf, c1 in 1:Nc
        for s2 in 1:2, f2 in 1:Nf, c2 in 1:Nc
            i = mapping_coordinate(n1, c1, s1, f1, Nc, Nf, N, order)
            j = mapping_coordinate(n2, c2, s2, f2, Nc, Nf, N, order)

            #list=[0.25*(op1=="density")*(op2=="density"),"Id",i,"Id",j]
            #ampo += tuple(list...)
            #list=[0.5*(op1=="density")*(op2=="density"),"Sz",i,"Id",j]
            #ampo += tuple(list...)
            #list=[0.5*(op1=="density")*(op2=="density"),"Id",i,"Sz",j]
            #ampo += tuple(list...)
            list=[(op1=="density")*(op2=="density"),"Sz",i,"Sz",j]
            ampo += tuple(list...)

            #list=[-(-1)^s2*0.5*(op1=="density")*(op2=="chiral"),"Id",i,"Sz",j]
            #ampo += tuple(list...)
            list=[-(-1)^s2*(op1=="density")*(op2=="chiral"),"Sz",i,"Sz",j]
            ampo += tuple(list...)

            #list=[-(-1)^s1*0.5*(op1=="chiral")*(op2=="density"),"Sz",i,"Id",j]
            #ampo += tuple(list...)
            list=[-(-1)^s1*(op1=="chiral")*(op2=="density"),"Sz",i,"Sz",j]
            ampo += tuple(list...)

            list=[(-1)^s1*(-1)^s2*(op1=="chiral")*(op2=="chiral"),"Sz",i,"Sz",j]
            ampo += tuple(list...)

        end
    end
  
    H = MPO(ampo,sites)
    return H
  end
