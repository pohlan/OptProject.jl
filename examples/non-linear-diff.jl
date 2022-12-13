using OptProject, PyPlot, Optim

function make_ϕ(u,a,b,p,Δx,f)
    du_1   = abs((u[1]-a)/Δx)^(p+1)
    du_end = abs((b-u[end])/Δx)^(p+1)
    du_mid = sum(abs.((u[2:end]-u[1:end-1])/Δx).^(p+1))
    uf     = sum(u.*f[2:end-1]) + (a*f[1] + b*f[end])/2
    return  Δx/(p+1) * (du_1 + du_mid + du_end) - uf * Δx
end
    # gradient
    #function dϕ(u)
function make_dϕ(u,a,b,p,f,Δx)
    ϵ            = eps()
    ddu          = zeros(length(u))
    ddu[2:end-1] = abs.(u[2:end-1] .- u[1:end-2] .+ ϵ).^(p-1) .* (u[2:end-1] .- u[1:end-2]) .-
                   abs.(u[3:end]   .- u[2:end-1] .+ ϵ).^(p-1) .* (u[3:end]   .- u[2:end-1])
    ddu[1]       = abs(u[1]-a+ϵ)^(p-1) * (u[1]-a) - abs(u[2]-u[1]+ϵ)^(p-1) * (u[2]-u[1])
    ddu[end]     = abs(u[end]-u[end-1]+ϵ)^(p-1) * (u[end]-u[end-1]) - abs(b-u[end]+ϵ)^(p-1) * (b-u[end])
    dϕ           = ddu/Δx^p - f[2:end-1]*Δx
    return dϕ
end

function make_ddϕ(u,a,b,p)
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

p = 0.5
a = 1.0
b = 0.0
n  = 21
Δx = 1.0 / (n+1)   # grid spacing
f  = zeros(n+2); f[1:end] .= 10 # source term
ϕ = u -> make_ϕ(u,a,b,p,Δx,f)
dϕ = u -> make_dϕ(u,a,b,p,f,Δx)
ddϕ = u-> make_ddϕ(u,a,b,p)

x0 = zeros(n,1); #x0[3:5] .= 0.3

ustar, it = newton(ϕ, dϕ, ddϕ, x0, maxiters=10^4)

# ustar, nit = quasi_Newton(x0, ϕ, dϕ, method=:sr1, maxiters=10^5)
# ustar, nit = quasi_Newton(x0, ϕ, dϕ, method=:BFGS, tol=3e-8,maxiters=1000,μ=0.1,ρ=0.5)

# Q, c, d = build_Q_c(;a,b,n,f)
# fct(u) = only(0.5 * transpose(u) * Q * u - transpose(c)*u) + d
# dfct(u) = Q * u - c
# ustar = sr1bt(x0, fct, dfct)

# s = optimize(ϕ, dϕ, x0, LBFGS(), inplace=false, Optim.Options(g_tol=1e-8,iterations=10^5))
# ustar = s.minimizer


plot(0:Δx:1,[a;ustar;b])