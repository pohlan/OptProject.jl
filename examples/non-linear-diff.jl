using OptProject, PyPlot, Optim, LaTeXStrings
close("all") # close all open figures
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 34

function make_f(u,a,b,p,Δx,g)
    du_1   = abs((u[1]-a)/Δx)^(p+1)
    du_end = abs((b-u[end])/Δx)^(p+1)
    du_mid = sum(abs.((u[2:end]-u[1:end-1])/Δx).^(p+1))
    uf     = sum(u.*g[2:end-1]) + (a*g[1] + b*g[end])/2
    return  Δx/(p+1) * (du_1 + du_mid + du_end) - uf * Δx
end

function make_df(u,a,b,p,g,Δx)
    ϵ            = eps()
    ddu          = zeros(length(u))
    ddu[2:end-1] = abs.(u[2:end-1] .- u[1:end-2] .+ ϵ).^(p-1) .* (u[2:end-1] .- u[1:end-2]) .-
                   abs.(u[3:end]   .- u[2:end-1] .+ ϵ).^(p-1) .* (u[3:end]   .- u[2:end-1])
    ddu[1]       = abs(u[1]-a+ϵ)^(p-1) * (u[1]-a) - abs(u[2]-u[1]+ϵ)^(p-1) * (u[2]-u[1])
    ddu[end]     = abs(u[end]-u[end-1]+ϵ)^(p-1) * (u[end]-u[end-1]) - abs(b-u[end]+ϵ)^(p-1) * (b-u[end])
    df           = ddu/Δx^p - g[2:end-1]*Δx
    return df
end

function make_ddf(u,a,b,p,Δx)
    ϵ = eps()
    n = length(u)
    A = zeros(n,n)
    for j = 1:n
        if j == 1
            ul = a
            ur = u[2]
        elseif j == n
            ul = u[n-1]
            ur = b
        else
            ul = u[j-1]
            ur = u[j+1]
        end
        if j<n
            A[j,j+1] = -p/Δx^p * abs(ur-u[j]+ϵ)^(p-1)
        end
        if j>1
            A[j,j-1] = -p/Δx^p * abs(u[j]-ul+ϵ)^(p-1)
        end
        A[j,j] =  p/Δx^p * (abs(u[j]-ul+ϵ)^(p-1) + abs(ur-u[j]+ϵ)^(p-1))
    end
    return A
end

a  = 1.0
b  = 0.0

ps = [1.0,1.3]
ns = 41:40:361     #31:20:161
it     = Dict(:newton => zeros(length(ps),length(ns)),
              :sr1    => zeros(length(ps),length(ns)),
              :BFGS   => zeros(length(ps),length(ns)),
              :LBFGS  => zeros(length(ps),length(ns)))
t_wall = Dict(:newton => zeros(length(ps),length(ns)),
              :sr1    => zeros(length(ps),length(ns)),
              :BFGS   => zeros(length(ps),length(ns)),
              :LBFGS  => zeros(length(ps),length(ns)))

figure(figsize=(18,12))
for (i,p) in enumerate(ps)
    println("p=$p")
    for (j,n) in enumerate(ns)
        println("Running for n=$n...")
        Δx = 1.0 / (n+1)   # grid spacing
        x0 = zeros(n,1);
        g  = 10 * ones(n+2,1)
        f = u -> make_f(u,a,b,p,Δx,g)
        df = u -> make_df(u,a,b,p,g,Δx)
        ddf = u-> make_ddf(u,a,b,p,Δx)

        # "warmup" for measuring time
        if i == 1 && j ==1
            newton(f, df, ddf, x0, maxiters=5*10^4)
            quasi_Newton(x0, f, df, method=:sr1, maxiters=2*10^3)
            optimize(f, df, x0, LBFGS(), inplace=false, Optim.Options(g_tol=1e-8,iterations=10^4))
        end
        # Newton
        tic = Base.time()
        ustar_newton, it[:newton][i,j] = newton(f, df, ddf, x0, maxiters=3*10^3)
        t_wall[:newton][i,j] = Base.time() - tic

        # SR1
        ρ = 0.4; μ = 1e-4
        ust = zeros(n,1)
        while all(ust .== 0.0) && ρ < 1.0
            try
                tic = Base.time()
                ust,   it[:sr1][i,j]   = quasi_Newton(x0, f, df, method=:sr1, tol=1e-7, maxiters=3*10^3; ρ,μ)
                t_wall[:sr1][i,j] = Base.time() - tic
            catch
                ρ += 0.1; μ *= 10
                println("Re-trying sr1 for ρ=$ρ, μ=$μ")
            end
        end

        # BFGS
        ρ = 0.4; μ = 1e-4
        ust = zeros(n,1)
        while all(ust .== 0.0) && ρ < 1.0
            try
                tic = Base.time()
                ust,   it[:BFGS][i,j]  = quasi_Newton(x0, f, df, method=:BFGS, tol=1e-7, maxiters=3*10^3; ρ,μ)
                t_wall[:BFGS][i,j] = Base.time() - tic
            catch
                ρ += 0.1; μ *= 10
                println("Re-trying BFGS for ρ=$ρ, μ=$μ")
            end
        end

        # LBFGS from Optim package
        tic = Base.time()
        s = optimize(f, df, x0, LBFGS(), inplace=false, Optim.Options(g_tol=1e-7,iterations=10^4))
        if s.ls_success == false
            @warn "no success for LBFGS p=$p, n=$n."
        end
        t_wall[:LBFGS][i,j]  = Base.time() - tic
        it[:LBFGS][i,j]      = s.iterations

        if j==length(ns)
            plot(0:Δx:1,[a;ustar_newton;b],label="p=$p")
        end
    end
end
xlabel("x")
ylabel("u")
legend()
savefig("ux_newton_n-"*string(ns[end])*".jpeg")


for (i, p_i) in enumerate(ps)
    figure(1,figsize=(39,13))
    for k in keys(it)
        subplot(1,length(ps),i)
        plot(ns,it[k][i,:],"-o",label=k,lw=3,ms=10)
    end
    xlabel("n")
    ylabel("iterations")
    legend()
    title("p = $p_i")
    savefig("iters_vs_n.jpeg")

    figure(2,figsize=(38,14))
    for k in keys(it)
        subplot(1,length(ps),i)
        plot(ns,t_wall[k][i,:],"-o",label=k,lw=3,ms=10)
    end
    xlabel("n")
    ylabel("elapsed time (s)")
    legend()
    title("p = $p_i")
    savefig("twall_vs_n.jpeg")
end
