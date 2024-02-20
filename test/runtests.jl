using REIM
using Test
using Aqua

using CSV
using DataFrames
using Distributions
using JSON
using LinearAlgebra
using NLsolve
using Parameters
using Random
using Revise

include("test_model.jl")
#include("test_parameter_estimation.jl")



@testset "REIM" begin
    @test test_sqrt_p0_eps0_div_q0()

    test_∂qs_∂Δϕ_nd()
    test_∂qs_∂Δγ_nd()

    @test test_surface_tension_PF6K()
    @test test_surface_coverage_vacancy_PF6K()
    @test test_surface_coverage_species_PF6K()
    @test test_capacitance_NaF()

    @test test_surface_tension_PF6K_nd()
    @test test_surface_coverage_vacancy_PF6K_nd()
    @test test_surface_coverage_species_PF6K_nd()
    @test test_capacitance_NaF_nd()

    #@test test_parameter_estimation_capacitance_data()
end

@testset "Aqua" begin
    Aqua.test_all(
        REIM;
        ambiguities=false,#(exclude=[], broken=false),
        stale_deps=false,#(ignore=[],),
        deps_compat=false,#(ignore=[],),
        piracies=true,
      )

end
