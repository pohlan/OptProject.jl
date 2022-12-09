function newton(f, df, ddf, xo; tol=1e-8, maxiters=10^4)
    xk = xo
    xs = [xk]
    for it = 1:maxiters
        xk = xk - ddf(xk) \ df(xk)
        push!(xs, xk)
        if abs(df(xk)) < tol
            if all(eig(ddf(xk)) .< 0)
                warning("Solution is a maximum.")
            elseif any(eig(ddf(xk)) .<= 0)
                warning("Solution is neither a maximum nor a minimum.")
            end
            return xk, xs
        end
    end
    error("No solution after $maxiters iterations")
end # function
