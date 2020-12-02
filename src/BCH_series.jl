module BCH_series

export BCH_denominator, BCH_coefficient

using Primes 

BCH_denominator(n::Int; T::Type{I}=Int)  where I<:Integer = 
    factorial(I(n))*prod([I(p)^(floor(Int, log(p, sum(digits(n, base=p)))))
        for p in 2:n if isprime(p)])

            
function BCH_coefficient(q::Vector{Int}; Afirst::Bool=true, T::Type{I}=Int) where I<:Integer
    m = length(q)
    N = sum(q)
    d = BCH_denominator(N, T=I)
    C = zeros(I, N, N)
    Acurrent = Afirst
    if iseven(m)
        Acurrent = !Acurrent
    end
    n = 0
    for i = m:-1:1
        for r = 1:q[i]
            n += 1
            h = zero(I)
            if i==m
                h = div(d, factorial(I(n)))
            elseif Acurrent && i==m-1
                h = div(d, factorial(I(r))*factorial(I(q[i+1])))
            end
            C[1, n] = h
            for k = 2:n-1
                h = zero(I)
                for j=1:r
                    if n>j && C[k-1,n-j]!=0
                        h += div(C[k-1, n-j], factorial(I(j)))
                    end
                end                     
                if Acurrent && i<=m-1
                    for j=1:q[i+1]
                        if n>r+j && C[k-1, n-r-j]!=0
                            h += div(C[k-1, n-r-j], factorial(I(r))*factorial(I(j)))
                        end
                    end
                end
                C[k, n] = h
            end
            C[n, n] = d 
        end
        Acurrent = !Acurrent
    end
    sum([div((-1)^(k+1)*C[k, N], k) for k=1:N])//d
end


function BCH_coefficient(b::Vector{Int}, q::Vector{Int}; T::Type{I}=Int) where I<:Integer
    m = length(q)
    @assert m==length(b)
    N = sum(q)
    d = BCH_denominator(N, T=I)
    S = zeros(Int, m)
    S[m] = 1
    for i = m-1:-1:1
        @assert b[i]!=b[i+1]
        S[i] = b[i]<b[i+1] ? S[i+1]+1 : 1
    end
    C = zeros(I, N, N)
    n = 0
    for i = m:-1:1
        for r = 1:q[i]
            n += 1
            h = zero(I)
            if m-i+1==S[i]
                x = factorial(I(r))
                for s = 1:S[i]-1
                    x *= factorial(I(q[i+s]))
                end
                h = div(d, x)
            end
            C[1, n] = h
            for k = 2:n-1
                h = zero(I)
                for j=1:r
                    if n>j && C[k-1,n-j]!=0
                        h += div(C[k-1, n-j], factorial(I(j)))
                    end
                end                     
                x = factorial(I(r))
                t = r
                for s = 1:S[i]-1
                    if s>=2
                        t += q[i+s-1]
                        x *= factorial(I(q[i+s-1]))
                    end
                    for j = 1:q[i+s]
                        if n>t+j && C[k-1, n-t-j]!=0
                            h += div(C[k-1, n-t-j], x*factorial(I(j)))
                        end
                    end
                end
                C[k, n] = h
            end
            C[n,n] = d 
        end
    end
    sum([div((-1)^(k+1)*C[k, N], k) for k=1:N])//d
end


function word2be(w::Vector{Int})
    k = 1
    for j=2:length(w) 
        if w[j]!=w[j-1]
            k += 1
        end
    end
    b = zeros(Int, k)
    e = zeros(Int, k)
    b[1] = w[1]
    k = 1
    f = 1
    for j=2:length(w)
        if w[j]==w[j-1]
            f += 1
        else
            e[k] = f
            k += 1
            b[k] = w[j]
            f = 1
        end
    end
    e[k] = f                
    b, e
end


function BCH_coefficient(W::String)
    if length(W) <= 19
        T = Int
    else
        T = BigInt
    end
    w = [x for x in W]
    b, e = word2be(w .- minimum(w))
    BCH_coefficient(b, e, T=T)
end


end # module
