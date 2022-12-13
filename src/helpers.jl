"""
    make_Q(;u1,un,n)

Setting up a linear system
Modified from straight.m
"""
function build_Q_c(;a,b,n,f)
    h = 1.0 / (n+1)   # grid spacing
    d = (a^2 + b^2) / (2 * h) - (f[1]*a+f[end]*b)/2 * h
    Q = zeros(n,n);  c = zeros(n,1)
    for j = 1:n
        if j==1
            Q[j,1:2] = [2 -1]
            c[j]     = a / h
        elseif j==n
            Q[j,n-1:n] = [-1 2]
            c[j]       = b / h
        else
            Q[j,j-1:j+1] = [-1 2 -1]
        end
    end
    c += f[2:end-1]*h
    Q  = Q / h
    return Q, c, d
end # function
