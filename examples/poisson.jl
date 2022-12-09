using OptProject, PyPlot
tr(A) = transpose(A)

u1 = 0.0
un = 0.0
n  = 41
f  = ones(n); f[18:22] .= 1 # source term

Q, c, d = build_Q_c(;u1,un,n,f)

# solve
ustar = Q \ c
fstar = 0.5 * tr(ustar) * Q * ustar .- tr(c) * ustar .+ d

# plot
figure(1)
plot(1:n, ustar)

