using Fluidum
using Test
using LinearAlgebra
using ForwardDiff
#=
@testset "equation of state" begin
    @test IdealQCD()==IdealQCD(3,2)
    @test isfinite(thermodynamic(1,IdealQCD()).pressure[1])
    @test isfinite(thermodynamic(0,IdealQCD()).pressure_derivative[1])
  
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure)==1
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure_derivative)==1
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure_hessian)==1

    @test isfinite(free_charm(1,0,Heavy_Quark())[1])
    @test round(pressure(1,Heavy_Quark()), sigdigits=10)==round(pressure(1,FluiduMEoS()), sigdigits=10)
    @test round(pressure(1,Heavy_Quark()), sigdigits=10)==round(thermodynamic(1,FluiduMEoS()).pressure[1], sigdigits=10)
    @test pressure(1,Heavy_Quark()) ≈ pressure(1,FluiduMEoS()) atol=0.01   
end

@testset "evolution" begin
    eos = Heavy_Quark()
    @test det(Fluidum.one_d_viscous_matrix([0.2,0.1,0,0,0,-0.1,0],2.,2.,0.5,0.5,0,0,0,1.,1.,0.1,0.1,0.1,1.,0)[1])!=0 
    @test typeof(Fluidum.runFluidum_fo(eos,fug_pars=Fluidum.FugacityPars(rdrop=8),DsT=0.24,maxtime=20.).fo)<:FreezeOutResult
    @test typeof(Fluidum.initialize_fields(Fluidum.TabulatedTrento(pwd()*"/../examples/ic_data/only_shift_BKG_full_order_changed.txt"),0,10).initial_field)<:Matrix{Float64}
    @test typeof(Fluidum.runFluidum(eos,DsT=0,maxtime=5.)(1)[1,:])<:Vector
    @test isfinite(Fluidum.runFluidum(eos,maxtime=2.)(2)[1,1])
end


@testset "observables" begin 
begin
    eos = Heavy_Quark()
    obs=Fluidum.compute_observables(eos,1.5,Tfo=0.156,save = true)
    @test isfile(Fluidum.get_filename(obs))
    @test typeof(obs.yield_th)<:Float64
end

@testset "plots" begin
    #Fluidum.plot_params(gui=true)
    eos = Heavy_Quark()
    obs=Fluidum.compute_observables(eos,1.5,Tfo=0.156,save = true)
    println("ok...")
    @test isfile(Fluidum.get_filename(obs))
    println("ok...")
    Fluidum.plot_spectra(obs,save=true)
    println("ok...")
    begin
    result=Fluidum.runFluidum(eos,DsT=0,maxtime=5.)
    Fluidum.plot_field(result,:temperature,tspan=(0.4,5.),save=true)
    #Fluidum.plot_field(result,tspan=(0.4,5.),save=true)
    println("ok...")
    end
    #rm(Fluidum.get_filename(obs)) 
end
=#
@testset "Cheb" begin
    function abs_matrix(diff, A)
        λ=eigen(A).values
        vec = eigen(A).vectors
        return vec*Diagonal(abs.(λ))*inv(vec)*diff
    end
    A = rand(2,2)
    A = A/max(eigen(A).values...)
    diff = rand(2)

    B = copy(A)
    diff2 = copy(diff)

    Fluidum.cheb_flux!(diff,A,15)
    @show diff, abs_matrix(diff2,B)
    @test isapprox(diff,abs_matrix(diff2,B),rtol=0.5)
    
end
#end


