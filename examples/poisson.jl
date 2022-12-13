using OptProject, PyPlot
tr(A) = transpose(A)

a = 0.0
b = 0.0
n  = 41
f  = ones(n+2); f[18:22] .= 1 # source term

Q, c, d = build_Q_c(;a,b,n,f)

# solve
ustar = Q \ c
fstar = 0.5 * tr(ustar) * Q * ustar .- tr(c) * ustar .+ d

# plot
figure(1)
plot(1:n, ustar)

