"""
    make_Q(;u1,un,n)

Setting up a linear system
Modified from straight.m
"""
function build_Q_c(;u1,un,n,f)
    h = 1.0 / (n+1)   # grid spacing
    d = (u1^2 + un^2) / (2 * h)
    Q = zeros(n,n);  c = zeros(n,1)
    for j = 1:n
        if j==1
            Q[j,1:2] = [2 -1]
            c[j]     = u1 / h + h^2 * f[1]/2
        elseif j==n
            Q[j,n-1:n] = [-1 2]
            c[j]       = un / h + h^2 * f[end]/2
        else
            Q[j,j-1:j+1] = [-1 2 -1]
            c[j]         = f[j]
        end
    end
    Q = Q / h
    return Q, c, d
end # function
