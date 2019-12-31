function block_compact_POP(x,f,g,h,k,r)
    n=length(x)
    m=length(g)
    l=length(h)

    dg = Array{UInt16}(undef, m)
    dh = Array{UInt16}(undef,l)

    rf=binomial(k+n,n)
    rg=Array{UInt16}(undef, m)
    rh=Array{UInt16}(undef,l)

    v=monomials(x, 0:k)
    v=reverse(v)

    Uf=union(2*v,monomials(f))
    Ug=Array{Any}(undef, m)
    Uh=Array{Any}(undef, l)

    Tf=v*v'
    Tg=Array{Any}(undef, m)
    Th=Array{Any}(undef, l)

    length_ind_block_f=Array{UInt16}(undef, 1)
    length_ind_block_g=Array{UInt16}(undef, m)
    length_ind_block_h=Array{UInt16}(undef, l)

    ind_block_f=Array{Any}(undef, 1)
    ind_block_g=Array{Any}(undef, m)
    ind_block_h=Array{Any}(undef, l)


    for i = 1:m
        dg[i] = ceil(UInt16,maxdegree(g[i])/2)
        rg[i]=binomial(k-dg[i]+n,n)
        Uf=union(Uf,monomials(g[i]))
        Tg[i]=g[i].*v[1:rg[i]]*v[1:rg[i]]'
    end



    for i = 1:l
        dh[i] = ceil(UInt16,maxdegree(h[i])/2)
        rh[i] = binomial(k-dh[i]+n,n)
        Uf=union(Uf,monomials(h[i]))
        Th[i]=h[i].*v[1:rh[i]]*v[1:rh[i]]'
    end


    for t=1:r
        U=Uf
        for i=1:m
            if t>1 && isempty(Ug[i])==0
                U=union(U,Ug[i])
            end
        end
        for i=1:l
            if t>1 && isempty(Uh[i])==0
                U=union(U,Uh[i])
            end
        end

        Cf=zeros(rf,rf)
        Cf[findall(p -> p in U, Tf)].=1
        Cf[diagind(Cf)] .= 0

        Cg=Array{Any}(undef, m)
        Ch=Array{Any}(undef, l)

        for i=1:m
            Cg[i]=zeros(rg[i],rg[i])

            for p=1:rg[i]
                for q=1:rg[i]
                    if isintersecting(monomials(Tg[i][p,q]),U)==1
                        Cg[i][p,q]=1
                    end
                end
            end
            Cg[i][diagind(Cg[i])] .= 0
        end


        for i=1:l
            Ch[i]=zeros(rh[i],rh[i])

            for p=1:rh[i]
                for q=1:rh[i]
                    if isintersecting(monomials(Th[i][p,q]),U)==1
                        Ch[i][p,q]=1
                    end
                end
            end
            Ch[i][diagind(Ch[i])].= 0
        end

        ind_block_f=connected_components(SimpleGraph(Cf))
        length_ind_block_f=length(ind_block_f)



        for i=1:m
            ind_block_g[i]=connected_components(SimpleGraph(Cg[i]))
            length_ind_block_g[i]=length(ind_block_g[i])
        end




        for i=1:l        
            ind_block_h[i]=connected_components(SimpleGraph(Ch[i]))
            length_ind_block_h[i]=length(ind_block_h[i])
        end 

        if r >=2
            for j=1:length_ind_block_f
                Uf=union(Uf,vec(v[ind_block_f[j]]*v[ind_block_f[j]]'))
            end
            

            for i=1:m
                for j=1:length_ind_block_g[i]
                    mono=monomials(g[i])
                    Ug[i]=vec(mono[1]*v[ind_block_g[i][j]]*v[ind_block_g[i][j]]')
                    for p=2:length(mono)
                        Ug[i]=union(Ug[i],vec(mono[p]*v[ind_block_g[i][j]]*v[ind_block_g[i][j]]'))
                    end
                end
            end

            for i=1:l
                for j=1:length_ind_block_h[i]
                    mono=monomials(h[i])
                    Uh[i]=vec(mono[1]*v[ind_block_h[i][j]]*v[ind_block_h[i][j]]')
                    for p=2:length(mono)
                        Uh[i]=union(Uh[i],vec(mono[p]*v[ind_block_h[i][j]]*v[ind_block_h[i][j]]'))
                    end
                end
            end
        end

    end








    model = Model(with_optimizer(Mosek.Optimizer, QUIET=true))

    wSOS=0

    G=Array{Any}(undef, length_ind_block_f)
    for i=1:length_ind_block_f
        G[i]= @variable(model, [1:length(ind_block_f[i]), 1:length(ind_block_f[i])],PSD)
        wSOS +=v[ind_block_f[i]]'*G[i]*v[ind_block_f[i]]   
    end


    Gg=Array{Any}(undef, m)
    for j=1:m
        Gg[j]=Array{Any}(undef, length_ind_block_g[j])
        for i=1:length_ind_block_g[j]
            Gg[j][i]= @variable(model, [1:length(ind_block_g[j][i]), 1:length(ind_block_g[j][i])],PSD)
            wSOS += g[j]*v[ind_block_g[j][i]]'*Gg[j][i]*v[ind_block_g[j][i]]   
        end
    end


    Gh=Array{Any}(undef, l)
    length_Gh=Array{Any}(undef, l)
    vh=Array{Any}(undef, l)
    for j=1:l
        Gh[j]=Array{Any}(undef, length_ind_block_h[j])
        length_Gh[j]=Array{Any}(undef, length_ind_block_h[j])
        vh[j]=Array{Any}(undef, length_ind_block_h[j])
        for i=1:length_ind_block_h[j]
            vh[j][i]= unique(vec(v[ind_block_h[j][i]]*v[ind_block_h[j][i]]'))
            length_Gh[j][i]= length(vh[j][i])
            Gh[j][i]=  @variable(model, [1:length_Gh[j][i]])
            wSOS += h[j]*vh[j][i]'*Gh[j][i]   
        end
    end


    @variable(model, lambda)

    @constraint(model, coefficients(f-lambda - wSOS).== 0)

    @objective(model, Max, lambda)

    optimize!(model)

    opt_val = value(lambda)
    
    println("Index of block for f : ", ind_block_f)
    println("------------------------------------")
    println("Indices of block for g : ", ind_block_g)
    println("------------------------------------")
    println("Indices of block for h ; ", ind_block_h)
    println("------------------------------------")
    println("Termination status = ", termination_status(model))

    println("Optimal value = ",opt_val)
    
    have_sol=0
    
    println("------------------------------------")
    G0=zeros(rf,rf)
    for i=1:length_ind_block_f
        G0[ind_block_f[i],ind_block_f[i]]=value.(G[i]) 
    end


    
    rk=rf-rank(G0, 1e-4)

    println("Rank of moment matrix = ", rk)

    if rk >0
        # extraction of Henrion and Lasserre
        F = eigen(G0)
        V = F.vectors

        Ix=sortperm(F.values)

        V=V[:,Ix[1:rk]]
        V=Matrix(V')
        V= rref_with_pivots!(V,1e-3);
        U=V[1]

        U=Matrix(U')
        # Figure out multiplying matrices using YALMIP code
        w=v[V[2]];
        N=zeros(length(V[2]),rk,n)
        for i in 1:n
            xw=x[i]*w
            kk=indexin(xw,v)
            N[:,:,i]=U[kk,:]
        end



        # Create random convex combination
        rands = rand(n,1);rands = rands/sum(rands);
        M = zeros(length(V[2]),rk);
        for i in 1:n
            M+=rands[i]*N[:,:,i];
        end

        F= schur(M);
        L=F.Z
        # Extract solution
        for i in 1:rk
            solution=[]
            for j = 1:n
                solution=[solution;L[:,i]'*N[:,:,j]*L[:,i]];
            end
            println("------------------------------------")
            println("Solution = ",solution)
            flag=1
            check=opt_val-polynomial(f)(x => solution)
            println("check lower bound  = ",check)
            if abs(check)>1e-3
                flag=0
            end
            for i in 1:m
                check=polynomial(g[i])(x => solution)
                println("check inequality ",i," = ",check) 
                if check<-1e-3
                flag=0
                end
            end

            for i in 1:l
                check=polynomial(h[i])(x => solution)
                println("check equality ",i," = ",check)
                if abs(check)>1e-3
                flag=0
                end
            end
            if flag ==1
                have_sol=1
            end

        end
    end
    
return opt_val, have_sol
end


function block_noncompact_POP(x,f,g,h,eps,k,r)
    n=length(x)
    m=length(g)
    l=length(h)

    if m==0 && l==0
        d = ceil(UInt16,maxdegree(f)/2)
    else
        d = floor(UInt16,maxdegree(f)/2) + 1
    end
    dg = Array{UInt16}(undef, m)
    dh = Array{UInt16}(undef,l)

    rf=binomial(k+d+n,n)
    rg=Array{UInt16}(undef, m)
    rh=Array{UInt16}(undef,l)

    v=monomials(x, 0:d+k)
    v=reverse(v)


    theta=1+x'*x

    Uf=union(2*v,monomials(theta^k*f))
    Ug=Array{Any}(undef, m)
    Uh=Array{Any}(undef, l)

    Tf=v*v'
    Tg=Array{Any}(undef, m)
    Th=Array{Any}(undef, l)

    length_ind_block_f=Array{UInt16}(undef, 1)
    length_ind_block_g=Array{UInt16}(undef, m)
    length_ind_block_h=Array{UInt16}(undef, l)

    ind_block_f=Array{Any}(undef, 1)
    ind_block_g=Array{Any}(undef, m)
    ind_block_h=Array{Any}(undef, l)


    for i = 1:m
        dg[i] = ceil(UInt16,maxdegree(g[i])/2)
        rg[i]=binomial(k+d-dg[i]+n,n)
        Uf=union(Uf,monomials(g[i]))
        Tg[i]=g[i].*v[1:rg[i]]*v[1:rg[i]]'
    end



    for i = 1:l
        dh[i] = ceil(UInt16,maxdegree(h[i])/2)
        rh[i] = binomial(k+d-dh[i]+n,n)
        Uf=union(Uf,monomials(h[i]))
        Th[i]=h[i].*v[1:rh[i]]*v[1:rh[i]]'
    end


    for t=1:r
        U=Uf
        for i=1:m
            if t>1 && isempty(Ug[i])==0
                U=union(U,Ug[i])
            end
        end
        for i=1:l
            if t>1 && isempty(Uh[i])==0
                U=union(U,Uh[i])
            end
        end

        Cf=zeros(rf,rf)
        Cf[findall(p -> p in U, Tf)].=1
        Cf[diagind(Cf)] .= 0

        Cg=Array{Any}(undef, m)
        Ch=Array{Any}(undef, l)

        for i=1:m
            Cg[i]=zeros(rg[i],rg[i])

            for p=1:rg[i]
                for q=1:rg[i]
                    if isintersecting(monomials(Tg[i][p,q]),U)==1
                        Cg[i][p,q]=1
                    end
                end
            end
            Cg[i][diagind(Cg[i])] .= 0
        end


        for i=1:l
            Ch[i]=zeros(rh[i],rh[i])

            for p=1:rh[i]
                for q=1:rh[i]
                    if isintersecting(monomials(Th[i][p,q]),U)==1
                        Ch[i][p,q]=1
                    end
                end
            end
            Ch[i][diagind(Ch[i])].= 0
        end

        ind_block_f=connected_components(SimpleGraph(Cf))
        length_ind_block_f=length(ind_block_f)



        for i=1:m
            ind_block_g[i]=connected_components(SimpleGraph(Cg[i]))
            length_ind_block_g[i]=length(ind_block_g[i])
        end




        for i=1:l        
            ind_block_h[i]=connected_components(SimpleGraph(Ch[i]))
            length_ind_block_h[i]=length(ind_block_h[i])
        end 

        if r >=2
            for j=1:length_ind_block_f
                Uf=union(Uf,vec(v[ind_block_f[j]]*v[ind_block_f[j]]'))
            end

            for i=1:m
                for j=1:length_ind_block_g[i]
                    mono=monomials(g[i])
                    Ug[i]=vec(mono[1]*v[ind_block_g[i][j]]*v[ind_block_g[i][j]]')
                    for p=2:length(mono)
                        Ug[i]=union(Ug[i],vec(mono[p]*v[ind_block_g[i][j]]*v[ind_block_g[i][j]]'))
                    end
                end
            end

            for i=1:l
                for j=1:length_ind_block_h[i]
                    mono=monomials(h[i])
                    Uh[i]=vec(mono[1]*v[ind_block_h[i][j]]*v[ind_block_h[i][j]]')
                    for p=2:length(mono)
                        Uh[i]=union(Uh[i],vec(mono[p]*v[ind_block_h[i][j]]*v[ind_block_h[i][j]]'))
                    end
                end
            end
        end

    end








    model = Model(with_optimizer(Mosek.Optimizer, QUIET=true))

    wSOS=0

    G=Array{Any}(undef, length_ind_block_f)
    for i=1:length_ind_block_f
        G[i]= @variable(model, [1:length(ind_block_f[i]), 1:length(ind_block_f[i])],PSD)
        wSOS += v[ind_block_f[i]]'*G[i]*v[ind_block_f[i]]   
    end


    Gg=Array{Any}(undef, m)
    for j=1:m
        Gg[j]=Array{Any}(undef, length_ind_block_g[j])
        for i=1:length_ind_block_g[j]
            Gg[j][i]= @variable(model, [1:length(ind_block_g[j][i]), 1:length(ind_block_g[j][i])],PSD)
            wSOS += g[j]*v[ind_block_g[j][i]]'*Gg[j][i]*v[ind_block_g[j][i]]   
        end
    end


    Gh=Array{Any}(undef, l)
    length_Gh=Array{Any}(undef, l)
    vh=Array{Any}(undef, l)
    for j=1:l
        Gh[j]=Array{Any}(undef, length_ind_block_h[j])
        length_Gh[j]=Array{Any}(undef, length_ind_block_h[j])
        vh[j]=Array{Any}(undef, length_ind_block_h[j])
        for i=1:length_ind_block_h[j]
            vh[j][i]= unique(vec(v[ind_block_h[j][i]]*v[ind_block_h[j][i]]'))
            length_Gh[j][i]= length(vh[j][i])
            Gh[j][i]=  @variable(model, [1:length_Gh[j][i]])
            wSOS +=h[j]*vh[j][i]'*Gh[j][i]   
        end
    end


    @variable(model, lambda)

    @constraint(model, coefficients(theta^k*(f-lambda+eps*theta^d) - wSOS).== 0)

    @objective(model, Max, lambda)

    optimize!(model)

    opt_val = value(lambda)
    
    println("Index of block for f : ", ind_block_f)
    println("------------------------------------")
    println("Indices of block for g : ", ind_block_g)
    println("------------------------------------")
    println("Indices of block for h ; ", ind_block_h)
    println("------------------------------------")
    println("Termination status = ", termination_status(model))

    println("Optimal value = ",opt_val)
    
return opt_val
end


function block_noncompact_POP_with_lower_bound(x,f,g,h,eps,k,r)
	
    n=length(x)
    m=length(g)
    l=length(h)

    d = ceil(UInt16,maxdegree(f)/2)

    dg = Array{UInt16}(undef, m)
    dh = Array{UInt16}(undef,l)

    rf=binomial(k+d+n,n)
    rg=Array{UInt16}(undef, m)
    rh=Array{UInt16}(undef,l)

    v=monomials(x, 0:d+k)
    v=reverse(v)


    theta=1+x'*x

    Uf=union(2*v,monomials(theta^k*f))
    Ug=Array{Any}(undef, m)
    Uh=Array{Any}(undef, l)

    Tf=v*v'
    Tg=Array{Any}(undef, m)
    Th=Array{Any}(undef, l)

    length_ind_block_f=Array{UInt16}(undef, 1)
    length_ind_block_g=Array{UInt16}(undef, m)
    length_ind_block_h=Array{UInt16}(undef, l)

    ind_block_f=Array{Any}(undef, 1)
    ind_block_g=Array{Any}(undef, m)
    ind_block_h=Array{Any}(undef, l)


    for i = 1:m
        dg[i] = ceil(UInt16,maxdegree(g[i])/2)
        rg[i]=binomial(k+d-dg[i]+n,n)
        Uf=union(Uf,monomials(g[i]))
        Tg[i]=g[i].*v[1:rg[i]]*v[1:rg[i]]'
    end



    for i = 1:l
        dh[i] = ceil(UInt16,maxdegree(h[i])/2)
        rh[i] = binomial(k+d-dh[i]+n,n)
        Uf=union(Uf,monomials(h[i]))
        Th[i]=h[i].*v[1:rh[i]]*v[1:rh[i]]'
    end


    for t=1:r
        U=Uf
        for i=1:m
            if t>1 && isempty(Ug[i])==0
                U=union(U,Ug[i])
            end
        end
        for i=1:l
            if t>1 && isempty(Uh[i])==0
                U=union(U,Uh[i])
            end
        end

        Cf=zeros(rf,rf)
        Cf[findall(p -> p in U, Tf)].=1
        Cf[diagind(Cf)] .= 0

        Cg=Array{Any}(undef, m)
        Ch=Array{Any}(undef, l)

        for i=1:m
            Cg[i]=zeros(rg[i],rg[i])

            for p=1:rg[i]
                for q=1:rg[i]
                    if isintersecting(monomials(Tg[i][p,q]),U)==1
                        Cg[i][p,q]=1
                    end
                end
            end
            Cg[i][diagind(Cg[i])] .= 0
        end


        for i=1:l
            Ch[i]=zeros(rh[i],rh[i])

            for p=1:rh[i]
                for q=1:rh[i]
                    if isintersecting(monomials(Th[i][p,q]),U)==1
                        Ch[i][p,q]=1
                    end
                end
            end
            Ch[i][diagind(Ch[i])].= 0
        end

        ind_block_f=connected_components(SimpleGraph(Cf))
        length_ind_block_f=length(ind_block_f)



        for i=1:m
            ind_block_g[i]=connected_components(SimpleGraph(Cg[i]))
            length_ind_block_g[i]=length(ind_block_g[i])
        end




        for i=1:l        
            ind_block_h[i]=connected_components(SimpleGraph(Ch[i]))
            length_ind_block_h[i]=length(ind_block_h[i])
        end 

        if r >=2
            for j=1:length_ind_block_f
                Uf=union(Uf,vec(v[ind_block_f[j]]*v[ind_block_f[j]]'))
            end

            for i=1:m
                for j=1:length_ind_block_g[i]
                    mono=monomials(g[i])
                    Ug[i]=vec(mono[1]*v[ind_block_g[i][j]]*v[ind_block_g[i][j]]')
                    for p=2:length(mono)
                        Ug[i]=union(Ug[i],vec(mono[p]*v[ind_block_g[i][j]]*v[ind_block_g[i][j]]'))
                    end
                end
            end

            for i=1:l
                for j=1:length_ind_block_h[i]
                    mono=monomials(h[i])
                    Uh[i]=vec(mono[1]*v[ind_block_h[i][j]]*v[ind_block_h[i][j]]')
                    for p=2:length(mono)
                        Uh[i]=union(Uh[i],vec(mono[p]*v[ind_block_h[i][j]]*v[ind_block_h[i][j]]'))
                    end
                end
            end
        end

    end








    model = Model(with_optimizer(Mosek.Optimizer, QUIET=true))

    wSOS=0

    G=Array{Any}(undef, length_ind_block_f)
    for i=1:length_ind_block_f
        G[i]= @variable(model, [1:length(ind_block_f[i]), 1:length(ind_block_f[i])],PSD)
        wSOS += v[ind_block_f[i]]'*G[i]*v[ind_block_f[i]]   
    end


    Gg=Array{Any}(undef, m)
    for j=1:m
        Gg[j]=Array{Any}(undef, length_ind_block_g[j])
        for i=1:length_ind_block_g[j]
            Gg[j][i]= @variable(model, [1:length(ind_block_g[j][i]), 1:length(ind_block_g[j][i])],PSD)
            wSOS += g[j]*v[ind_block_g[j][i]]'*Gg[j][i]*v[ind_block_g[j][i]]   
        end
    end


    Gh=Array{Any}(undef, l)
    length_Gh=Array{Any}(undef, l)
    vh=Array{Any}(undef, l)
    for j=1:l
        Gh[j]=Array{Any}(undef, length_ind_block_h[j])
        length_Gh[j]=Array{Any}(undef, length_ind_block_h[j])
        vh[j]=Array{Any}(undef, length_ind_block_h[j])
        for i=1:length_ind_block_h[j]
            vh[j][i]= unique(vec(v[ind_block_h[j][i]]*v[ind_block_h[j][i]]'))
            length_Gh[j][i]= length(vh[j][i])
            Gh[j][i]=  @variable(model, [1:length_Gh[j][i]])
            wSOS +=h[j]*vh[j][i]'*Gh[j][i]   
        end
    end


    @variable(model, lambda)

    @constraint(model, coefficients(theta^k*(f-lambda+eps*theta^d) - wSOS).== 0)

    @objective(model, Max, lambda)

    optimize!(model)

    opt_val = value(lambda)
    
    println("Index of block for f : ", ind_block_f)
    println("------------------------------------")
    println("Indices of block for g : ", ind_block_g)
    println("------------------------------------")
    println("Indices of block for h ; ", ind_block_h)
    println("------------------------------------")
    println("Termination status = ", termination_status(model))

    println("Optimal value = ",opt_val)
    
return opt_val
end


function adding_spherical_constraints(x,g,h,k,r)

    # small parameter
    eps = 1e-1

    n=length(x)
    m=length(g)
    l=length(h)

    dg = Array{UInt16}(undef, m)
    dh = Array{UInt16}(undef,l)

    for i = 1:m
        dg[i] = ceil(UInt16,maxdegree(g[i])/2)
    end

    for i = 1:l
        dh[i] = ceil(UInt16,maxdegree(h[i])/2)
    end

    d_max=maximum(union(union(dg,dh),1))

    # Define centers and square of radius
    a0=zeros(Float64,(n, 1)); a = Matrix{Float64}(I, n, n)

    println("Determine L0:")

    println("------------------------------------")

    # Polynomial to optimize 
    f = [x-a0]'*[x-a0]

    if m==0
        g=[f]
    else
        g=[g;f]
    end

    L0=PV_TSSOS.block_noncompact_POP_with_lower_bound(x,f,g,h,eps,k,r)
    
    println("------------------------------------")

    println("L0 = ", L0)

    println("====================================")

    # Define omegat, t=0,...,n
    omega0 = 0; omega = zeros(n)

    #inequalities polynomial
    g[end]=L0-f

    println("Determine omega",0,":")

    println("------------------------------------")

    obtain_data=PV_TSSOS.block_compact_POP(x,f,g,h,k+d_max,r)
    
    omega0=obtain_data.opt_val
    
    println("------------------------------------")

    println("omega",0," = ", omega0)
    

    println("====================================")

    #equalities polynomial
    if l==0
        h=[omega0-f]
    else
        h=[h;omega0-f]
    end


    #inequalities polynomial
    g=g[1:end-1]
    
    if obtain_data.have_sol==0
        for t=1:n

            println("Determine omega",t,":")

            println("------------------------------------")



            if t>1
                #equalities polynomial
                h = [h;omega[t-1]-f] ; l=length(h)
            end        

            f=[x-a[:,t]]'*[x-a[:,t]]




            obtain_data=PV_TSSOS.block_compact_POP(x,f,g,h,k+d_max,r)
            omega[t]=obtain_data.opt_val
            println("------------------------------------")
            println("omega",t," = ", omega[t])



            println("====================================")
            if obtain_data.have_sol==0
                break
            end
        end
    end
end


function isintersecting(v1, v2)
    flag=0
    for v in v1
        if v in v2 
            flag=1
            break
        end
    end
return flag 
end
