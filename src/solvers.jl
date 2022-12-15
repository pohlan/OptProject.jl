function newton(f, df, ddf, x0; tol=1e-8, maxiters=10^3)
    xk = x0
    for it = 1:maxiters
        xk = xk - ddf(xk) \ df(xk)
        if norm(df(xk)) < tol
            println("Solution found after $it iterations.")
            if !all(eigen(ddf(xk)).values .>= 0)
                warning("Solution is not a minimum.")
            end
            return xk, it
        end
    end
    error("No solution found after $maxiters iterations")
end  # function

"""
Quasi-Newton: Symmetric rank-one or BFGS
"""
function quasi_Newton(x0, f, df;
                      method=:BFGS,      # either :sr1 or :BFGS
                      tol=1e-8, maxiters=10^4,
                      # for backtracking:
                      μ=1e-4,       # sufficient decrease
                      ρ=0.5)        # backtracking by halving)
    n  = length(x0)
    xk = x0
    dfxk = df(xk)
    B = I(n)
    for k = 1:maxiters
        if norm(dfxk) < tol
            return xk, k
        end
        pk = - B \ dfxk

        if only(transpose(dfxk) * pk) >= 0.0
            pk = -dfxk
        end
        # update xk
        α = bt(xk,pk,dfxk,f;μ,ρ)
        oldxk = xk
        olddfxk = dfxk
        xk = xk + α * pk

        # update Bk
        dfxk = df(xk)
        sk = xk - oldxk
        yk = dfxk - olddfxk
        if method == :sr1               # symmetric rank 1 update
            v  = yk - B * sk
            denom = transpose(v) * sk
            if denom == 0
                return xk, k
            end
            B = B + v * transpose(v / denom)
        elseif method == :BFGS          # BFGS update
            w = B*sk
            denom = transpose(sk) * w
            if denom == 0
                return xk, k
            end
            B = B - w * transpose(w) ./ denom + (yk*transpose(yk)) ./ (transpose(yk)*sk)
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
