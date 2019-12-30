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

    d_max=maximum([dg;dh;1])

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

    L0=block_noncompact_POP_with_lower_bound(x,f,g,h,eps,k,r)
    
    println("------------------------------------")

    println("L0 = ", L0)

    println("====================================")

    # Define omegat, t=0,...,n
    omega0 = 0; omega = zeros(n)

    #inequalities polynomial
    g[end]=L0-f

    println("Determine omega",0,":")

    println("------------------------------------")

    omega0=block_compact_POP(x,f,g,h,k+d_max,r)

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

    for t=1:n

        println("Determine omega",t,":")

        println("------------------------------------")



        if t>1
            #equalities polynomial
            h = [h;omega[t-1]-f] ; l=length(h)
        end        

        f=[x-a[:,t]]'*[x-a[:,t]]




        omega[t]=block_compact_POP(x,f,g,h,k+d_max,r)
        println("------------------------------------")
        println("omega",t," = ", omega[t])



        println("====================================")
    end

end