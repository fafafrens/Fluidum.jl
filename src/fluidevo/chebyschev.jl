using BenchmarkTools
using Fluidum
function cheb_coeff(j;N=3,f=abs)
    @assert isless(j,N)
    arg1(k) = π*(k+0.5)/N
    arg2(k,j) = π*j*(k+0.5)/N
    sum = 0
    [sum+=f(cos(arg1(k)))*cos(arg2(k,j)) for k in 0:N-1]
    return 2/N*sum
end

[cheb_coeff(i,N=4) for i in 0:3]

function cheb_pol(n,x)
    @assert isless(-1,n)
    if n == 0
        return one(x)
    elseif n==1
        return x
    else
        return 2*x*cheb_pol(n-1,x)-cheb_pol(n-2,x)
    end
end

[cheb_pol(i,-0.5) for i in 0:3]

function cheb_approx(x;N=3, f= abs)
    sum = zero(x)
    [sum+=cheb_coeff(j,N=N,f=f)*cheb_pol(j,x) for j in 0:N-1]
    return sum -one(x)*0.5*cheb_coeff(0,N=N,f=f)
end

cheb_approx([0.5 0.5; 0 0],N=6)
cheb_approx([0.5 0.5; 0 0],N=10)
cheb_approx([0.5 0.5; 0 0],N=20)
cheb_approx([0.8 -0.8; 0 0.5],N=30)

poly(x) = 0.5*one(x)+0.5*x*x
A=[0.5 -0.5; 0.3 0.5]
poly([0.5 -0.5; 0 0.5])
@benchmark cheb_coeff($2)

temp = [0 0; 0 0]
diff = [1 1]
Fluidum.upwindflux!(diff,[1 1],A,temp)
diff
temp


temp = [0 0; 0 0]
diff = [1 1]
Fluidum.upwindflux_Chebyschev!(diff,[1 1],A,temp)

diff
temp

