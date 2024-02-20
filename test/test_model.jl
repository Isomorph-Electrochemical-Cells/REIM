function test_sqrt_p0_eps0_div_q0()
    # check consistency of nondimensional parameter
    return REIM.sqrt_p0_eps0_div_q0 ≈ 1.0 &&
           REIM.sqrt_p0_eps0_div_q0 ≈ REIM.V0/REIM.l0 * sqrt(REIM.ϵ0/REIM.p0)
end

function test_∂qs_∂Δϕ_nd()
    # load and preprocess model data
    dict_input_params = read_model_parameters("data/input/PF6K2.json")
    params = preprocess_model_parameters(dict_input_params["model_parameters"])
    params = expected_model_params(params)
    nondimensionalize!(params)

    Δϕ = 0.25
    Δγ = REIM.surface_tension_difference_nd(Δϕ, params)

    fd = REIM.∂qs_∂Δϕ_nd_fd(Δϕ, Δγ, params)
    ex = REIM.∂qs_∂Δϕ_nd(Δϕ, Δγ, params)

    @test fd ≈ ex rtol=1e-5
end

function test_∂qs_∂Δγ_nd()
    # load and preprocess model data
    dict_input_params = read_model_parameters("data/input/PF6K2.json")
    params = preprocess_model_parameters(dict_input_params["model_parameters"])
    params = expected_model_params(params)
    nondimensionalize!(params)

    Δϕ = 0.25
    Δγ = REIM.surface_tension_difference_nd(Δϕ, params)

    fd = REIM.∂qs_∂Δγ_nd_fd(Δϕ, Δγ, params)
    ex = REIM.∂qs_∂Δγ_nd(Δϕ, Δγ, params)

    @test fd ≈ ex rtol=1e-5
end


function test_surface_tension_PF6K()
    # load and preprocess model data
    dict_input_params = read_model_parameters("data/input/PF6K2.json")
    params = preprocess_model_parameters(dict_input_params["model_parameters"])
    params = expected_model_params(params)

    # load reference data
    df_surface_tension = CSV.read("data/reference/landstorfer_2016_fig9_surface_tension.csv", DataFrame; header =["Δϕ", "Δγ"])

    # evaluate surface tension at reference Δϕ values
    Δγ_values = [ REIM.surface_tension_difference(Δϕ, params).*10^3 for Δϕ in df_surface_tension.Δϕ]
    # check that the maximum relative error is below a small threshold value
    return maximum(abs.((df_surface_tension.Δγ-Δγ_values)./df_surface_tension.Δγ)) < 1e-2
end

function test_surface_tension_PF6K_nd()
    # load and preprocess model data
    dict_input_params = read_model_parameters("data/input/PF6K2.json")
    params = preprocess_model_parameters(dict_input_params["model_parameters"])
    params = expected_model_params(params)
    nondimensionalize!(params)

    # load reference data
    df_surface_tension = CSV.read("data/reference/landstorfer_2016_fig9_surface_tension.csv", DataFrame; header =["Δϕ", "Δγ"])
    # nondimensionalize the reference data
    df_surface_tension.Δϕ ./= REIM.V0
    df_surface_tension.Δγ ./= REIM.γ0

    # evaluate surface tension at reference Δϕ values
    Δγ_values = [ REIM.surface_tension_difference_nd(Δϕ, params).*10^3 for Δϕ in df_surface_tension.Δϕ]
    # check that the maximum relative error is below a small threshold value
    return maximum(abs.((df_surface_tension.Δγ-Δγ_values)./df_surface_tension.Δγ)) < 1e-2
end

function test_surface_coverage_vacancy_PF6K()
    # load and preprocess model data
    dict_input_params = read_model_parameters("data/input/PF6K2.json")
    params = preprocess_model_parameters(dict_input_params["model_parameters"])
    params = expected_model_params(params)

    # load reference data
    df_coverage_vacancy = CSV.read("data/reference/landstorfer_2016_fig8_coverage_Vacancies.csv", DataFrame; header =["Δϕ", "coverage"])

    # evaluate surface vacancy coverage at reference Δϕ values
    coverage_values = [ REIM.surface_coverage_vacancy(Δϕ, params) for Δϕ in df_coverage_vacancy.Δϕ]
    return maximum(abs.((df_coverage_vacancy.coverage.-coverage_values))) < 1e-2
end

function test_surface_coverage_vacancy_PF6K_nd()
    # load and preprocess model data
    dict_input_params = read_model_parameters("data/input/PF6K2.json")
    params = preprocess_model_parameters(dict_input_params["model_parameters"])
    params = expected_model_params(params)
    nondimensionalize!(params)

    # load reference data
    df_coverage_vacancy = CSV.read("data/reference/landstorfer_2016_fig8_coverage_Vacancies.csv", DataFrame; header =["Δϕ", "coverage"])
    # nondimensionalize the reference data
    df_coverage_vacancy.Δϕ ./= REIM.V0

    # evaluate surface vacancy coverage at reference Δϕ values
    coverage_values = [ REIM.surface_coverage_vacancy_nd(Δϕ, params) for Δϕ in df_coverage_vacancy.Δϕ]
    return maximum(abs.((df_coverage_vacancy.coverage.-coverage_values))) < 1e-2
end


function test_surface_coverage_species_PF6K()
    # load and preprocess model data
    dict_input_params = read_model_parameters("data/input/PF6K2.json")
    params = preprocess_model_parameters(dict_input_params["model_parameters"])
    params = expected_model_params(params)

    # load reference data
    df_coverage_species_h = CSV.read("data/reference/landstorfer_2016_fig8_coverage_H.csv",
    DataFrame; header =["Δϕ", "coverage"])
    df_coverage_species_oh = CSV.read("data/reference/landstorfer_2016_fig8_coverage_OH.csv",
    DataFrame; header =["Δϕ", "coverage"])

    # evaluate surface species coverage at reference Δϕ values
    coverage_oh = [ REIM.surface_coverage_species(Δϕ, params)[2] for Δϕ in df_coverage_species_oh.Δϕ]
    coverage_h = [ REIM.surface_coverage_species(Δϕ, params)[3] for Δϕ in df_coverage_species_h.Δϕ]

    return (maximum(abs.((df_coverage_species_oh.coverage .- coverage_oh)))) < 4e-2 &&
           (maximum(abs.((df_coverage_species_h.coverage .- coverage_h)))) < 4e-2
end

function test_surface_coverage_species_PF6K_nd()
    # load and preprocess model data
    dict_input_params = read_model_parameters("data/input/PF6K2.json")
    params = preprocess_model_parameters(dict_input_params["model_parameters"])
    params = expected_model_params(params)
    nondimensionalize!(params)

    # load reference data
    df_coverage_species_h = CSV.read("data/reference/landstorfer_2016_fig8_coverage_H.csv",
    DataFrame; header =["Δϕ", "coverage"])
    df_coverage_species_oh = CSV.read("data/reference/landstorfer_2016_fig8_coverage_OH.csv",
    DataFrame; header =["Δϕ", "coverage"])
    # nondimensionalize reference data
    df_coverage_species_h.Δϕ ./= REIM.V0
    df_coverage_species_oh.Δϕ ./= REIM.V0

    # evaluate surface species coverage at reference Δϕ values
    coverage_oh = [ REIM.surface_coverage_species_nd(Δϕ, params)[2] for Δϕ in df_coverage_species_oh.Δϕ]
    coverage_h = [ REIM.surface_coverage_species_nd(Δϕ, params)[3] for Δϕ in df_coverage_species_h.Δϕ]

    return (maximum(abs.((df_coverage_species_oh.coverage .- coverage_oh)))) < 4e-2 &&
           (maximum(abs.((df_coverage_species_h.coverage .- coverage_h)))) < 4e-2
end

function test_capacitance_NaF()
    # load and preprocess model data
    dict_input_params = read_model_parameters("data/input/NaF2.json")
    params = preprocess_model_parameters(dict_input_params["model_parameters"])
    params = expected_model_params(params)

    # load reference data
    df_capacitance_c01 = CSV.read("data/reference/landstorfer_2016_fig15a_capacitance_vs_EvsSC_c01.csv", DataFrame; header =["EvsSCE", "capacitance"])
    # reference voltage of a calomel electrode
    reference_voltage = -0.97

    capacitance_c01 = [10^2*REIM.capacitance(EvsSCE-reference_voltage, params) for EvsSCE in df_capacitance_c01.EvsSCE]

    return maximum(abs.((df_capacitance_c01.capacitance .- capacitance_c01)./df_capacitance_c01.capacitance)) < 5e-2
end

function test_capacitance_NaF_nd()
    # load and preprocess model data
    dict_input_params = read_model_parameters("data/input/NaF2.json")
    params = preprocess_model_parameters(dict_input_params["model_parameters"])
    params = expected_model_params(params)
    nondimensionalize!(params)

    # load reference data
    df_capacitance_c01 = CSV.read("data/reference/landstorfer_2016_fig15a_capacitance_vs_EvsSC_c01.csv", DataFrame; header =["EvsSCE", "capacitance"])
    # reference voltage of a calomel electrode
    reference_voltage = -0.97

    # nondimensionalize reference data
    reference_voltage /= REIM.V0
    df_capacitance_c01.EvsSCE ./= REIM.V0
    df_capacitance_c01.capacitance ./= REIM.C0

    capacitance_c01 = [10^2*REIM.capacitance_nd(EvsSCE-reference_voltage, params) for EvsSCE in df_capacitance_c01.EvsSCE]

    return maximum(abs.((df_capacitance_c01.capacitance .- capacitance_c01)./df_capacitance_c01.capacitance)) < 5e-2
end
