


struct InverseFunction{N,F}
    fun::F
end 

function InverseFunction(f)
    InverseFunction{1,typeof(f)}(f)
end

function (func::InverseFunction{1,F})(x; autodiff=AutoForwardDiff(),u0=0.1*one(x)) where F 
    
    function f(u,p)
        func.fun(u)-x
    end 
    prob = NonlinearProblem{false}(f,u0;autodiff=autodiff) 
    
    return solve(prob, SimpleNewtonRaphson()).u
end

#this simple function define a decoration to a callble onject like 
#prova(x) = x^19 +x +5
#the inverse of this function can then be obtained just decorating the function call 
#InverseFunction(prova)(10.)


#prova3 = x-> x^2
#vect = (0.:1.:10.)
#InverseFunction.(Ref(prova)).(vect)

#more complicated examples are the following. 
#let consider the a two d function 

#prova2(x,y)= x^2+y^2

#let invert this function in the positive quadrant
#InverseFunction(x->prova2(x,10))(1.)




