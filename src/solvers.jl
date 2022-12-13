function newton(f, df, ddf, x0; tol=1e-8, maxiters=10^3)
    xk = x0
    for it = 1:maxiters
        xk = xk - ddf(xk) \ df(xk)
        if norm(df(xk)) < tol
            println("Solution found after $it iterations.")
            if all(eigen(ddf(xk)).values .< 0)
                warning("Solution is a maximum.")
            elseif any(eigen(ddf(xk)).values .<= 0)
                warning("Solution is neither a maximum nor a minimum.")
            end
            return xk, it
        end
    end
    error("No solution after $maxiters iterations")
end # function

"""
Quasi-Newton: Symmetric rank-one or BFGS
"""
function quasi_Newton(x0, f, df=nothing;
                      method=:BFGS,      # either :sr1 or :BFGS
                      tol=1e-8, maxiters=10^4,
                      # for backtracking:
                      μ=1e-4,       # sufficient decrease
                      ρ=0.5)        # backtracking by halving)
    xk = x0
    n  = length(xk)

    if df === nothing
        df = xk -> fdgrad(xk,f)   # finite difference approx. of derivative
    end

    B = I(n)
    dfxk = df(xk)
    for k = 1:maxiters
        if norm(dfxk) < tol
            println("Solution found after $k iterations.")
            return xk, k
        end
        pk = - B \ dfxk

        if only(transpose(dfxk) * pk) >= 0.0
            @warn "Non-descent direction at step $k .. revert to steepest-descent."
            pk = -dfxk
        end
        α = bt(xk,pk,dfxk,f;μ,ρ)
        oldxk = xk
        olddfxk = dfxk
        xk = xk + α * pk

        dfxk = df(xk)
        sk = xk - oldxk
        yk = dfxk - olddfxk

        if method == :sr1               # symmetric rank 1 update
            v  = yk - B * sk
            denom = transpose(v) * sk
            if denom == 0
                println("Solution found after $k iterations.")
                return xk, k
            end
            B = B + v * transpose(v / denom)
        elseif method == :BFGS          # BFGS update
            w = B*sk
            den = only(transpose(sk) * w)
            if den == 0 || only(transpose(yk)*sk) == 0
                println("Solution found after $k iterations.")
                return xk, k
            end
            B = B - w * transpose(w) ./ den + (yk*transpose(yk)) ./ (transpose(yk)*sk)
        end
    end
    error("No solution found after $maxiters iterations.")
end # function

"""
Backtracking
"""
function bt(xk,pk,dfxk,f;
            μ=1e-4,       # sufficient decrease
            ρ=0.5)        # backtracking by halving
    Dk = only(transpose(dfxk)*pk)
    fxk = f(xk)
    α   = 1.0             # α_0 = 1 because of Newton and by Thm 11.7
    while f(xk+α*pk) > fxk + μ*α*Dk
        α = ρ * α
    end
    return α
end # function

"""
Finite difference approximation to the gradient
"""
function fdgrad(xk,f)
    fxk = f(xk)
    n = length(xk)
    dfxk = zeros(n,1)
    h = sqrt(eps())      # 'simple' formula for h
    for i = 1:n
        xkh     = xk
        xkh[i]  = xk[i] + h
        dfxk[i] = (f(xkh) - fxk) / h
    end
    return dfxk
end # function
