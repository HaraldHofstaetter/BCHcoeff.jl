module BCH_series

export BCH_denominator, BCH_coefficient, BCH_extremal_p_valuation

using Primes 

BCH_denominator(n::Int; T::Type{I}=Int)  where I<:Integer = 
    factorial(I(n))*prod([I(p)^(floor(Int, log(p, sum(digits(n, base=p)))))
        for p in 2:n if isprime(p)])

            
function BCH_coefficient(q::Vector{Int}; Afirst::Bool=true, T::Type{I}=Int) where I<:Integer
    m = length(q)
    N = sum(q)
    d = BCH_denominator(N, T=I)
    F = [factorial(I(n)) for n=1:N]
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
                h = div(d, F[n])
            elseif Acurrent && i==m-1
                h = div(d, F[r]*F[q[i+1]])
            end
            C[1, n] = h
            for k = 2:n-1
                h = zero(I)
                for j=1:r
                    if n>j && C[k-1,n-j]!=0
                        h += div(C[k-1, n-j], F[j])
                    end
                end                     
                if Acurrent && i<=m-1
                    for j=1:q[i+1]
                        if n>r+j && C[k-1, n-r-j]!=0
                            h += div(C[k-1, n-r-j], F[r]*F[j])
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
    F = [factorial(I(n)) for n=1:N]
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
                x = F[r]
                for s = 1:S[i]-1
                    x *= F[q[i+s]]
                end
                h = div(d, x)
            end
            C[1, n] = h
            for k = 2:n-1
                h = zero(I)
                for j=1:r
                    if n>j && C[k-1,n-j]!=0
                        h += div(C[k-1, n-j], F[j])
                    end
                end                     
                x = F[r]
                t = r
                for s = 1:S[i]-1
                    if s>=2
                        t += q[i+s-1]
                        x *= F[q[i+s-1]]
                    end
                    for j = 1:q[i+s]
                        if n>t+j && C[k-1, n-t-j]!=0
                            h += div(C[k-1, n-t-j], x*F[j])
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


function BCH_extremal_p_valuation(n::Int, p::Int)
    α = digits(n, base=p)
    r = floor(Int, log(p, n))
    l = floor(Int, log(p, sum(α)))
    if l==0
        if n<p
            if n==1
                return [1]
            elseif n==2 || isodd(n)
                return [n-1, 1]
            else
                return [n-2, 2]
            end
        else
            k = p^(r-1)*(p-1)
            return [n-k, k]
        end
    elseif l==1
        h = p-1
        i = 0
        β = zeros(Int, r+1)
        while h>0
            β[i+1] = min(h, α[i+1])
            h -= β[i+1]
            i += 1
        end
        k = sum([β[j+1]*p^j for j=0:(i-1)])
        return [n-k, k]
    else
        if p==2 || isodd(n)
            m = p^l
        else
            m = p^l+1
        end
        b = vcat([Int[p^j for i=1:α[j+1]] for  j=0:length(α)-1]...)
        return vcat(sum(b[m:end]), reverse(b[1:m-1]))      
    end
end


end # module
