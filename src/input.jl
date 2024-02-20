function convert_to_normal(mean::T, std=zero(T)) where T<:Real
    return Normal(mean, std)
end

function convert_to_normal(input::Dict{String,T}) where {T}
    return Normal(input["mean"], input["std"])
end

function convert_to_truncated_normal(mean, std, lower, upper)
    return truncated(Normal(mean, std; lower, upper))
end

# function convert_to_truncated_normal(input::Dict{String,T}) where {T}
#     lb_value = haskey(input, "lower") && !isnan(input["lower"]) ? input["lower"] : -Inf
#     ub_value = haskey(input, "upper") && !isnan(input["lower"]) ? input["upper"] : Inf
#     return truncated(Normal(input["mean"], input["std"]); lower=lb_value, upper=ub_value)
# end

# function convert_to_truncated_normal(input::T) where {T<:AbstractFloat}
#      return truncated(Normal(input, 0.0); lower=input, upper=input)
# end

function convert_to_measurement(dict::Dict{String,T}) where {T}
    measurement(dict["mean"], dict["std"])
end

function convert_to_measurement(value::T) where {T<:Real}
    measurement(value, zero(value))
end


function convert_to_normal(measurement::Measurement)
    return Normal(Measurements.value(measurement), Measurements.uncertainty(measurement))
end

function read_model_parameters(input_file_path)
    dict_input_parameters = JSON.parsefile(input_file_path; dicttype=Dict, inttype=Int64, use_mmap=true)
    return dict_input_parameters
end

function expected_axis_array(axis_array::AxisArrays.AxisArray{Distr, 1, Vector{Normal{T}}}) where
    {T,Distr<:UnivariateDistribution}
    num_cols = size(axis_array)[1]
    names_col = [AxisArrays.axes(axis_array, 1)[i] for i in 1:length(AxisArrays.axes(axis_array, 1))]

    expected_axis_array = AxisArray(Vector{Float64}(undef, num_cols); col=names_col)
    [expected_axis_array[col] = mean(axis_array[col]) for col in 1:num_cols]

    return expected_axis_array
end

function expected_axis_array(axis_array::AxisArrays.AxisArray{Distr, 2, Matrix{Normal{T}}}) where
    {T,Distr<:UnivariateDistribution}

    (num_rows, num_cols) = size(axis_array)

    names_row = [AxisArrays.axes(axis_array, 1)[i] for i in 1:length(AxisArrays.axes(axis_array, 1))]
    names_col = [AxisArrays.axes(axis_array, 2)[i] for i in 1:length(AxisArrays.axes(axis_array, 2))]

    expected_axis_array = AxisArray(Matrix{Float64}(undef, num_rows, num_cols);
                                                    row=names_row,
                                                    col=names_col)
    [expected_axis_array[row, col] = mean(axis_array[row, col]) for
                                        row in 1:num_rows, col in 1:num_cols]

    return expected_axis_array
end

function expected_model_params(model_params_distr)
    ModelParameters{Float64,Float64,Float64}(
        electrode=expected_axis_array(model_params_distr.electrode),
        electrolyte=expected_axis_array(model_params_distr.electrolyte),
        species=expected_axis_array(model_params_distr.species))
end

function expected_reaction_params(model_params_distr)
    ReactionParameters{Float64}(reactions=expected_axis_array(model_params_distr.reactions))
end

function preprocess_model_parameters(params)
    ParameterDistribution = Normal{Float64} # Truncated{Normal{Float64},Continuous,Float64}
    ModelParameterDistribution = ModelParameters{ParameterDistribution,ParameterDistribution,ParameterDistribution}

    processed_params_distr = ModelParameterDistribution()

    species_property_names = [property for property in
                              AxisArrays.axes(processed_params_distr.species, 2)]


    # Electrode parameters
    for name in AxisArrays.axes(processed_params_distr.electrode, 1)
        processed_params_distr.electrode[col=name] = convert_to_normal(
                                                    params["electrode"][name])
    end

    # General electrolyte paremters
    for name in AxisArrays.axes(processed_params_distr.electrolyte, 1)
        processed_params_distr.electrolyte[col=name] = convert_to_normal(
                                                        params["electrolyte"][name])
    end

    # Determine total number of species (including the solvent)
    num_ionic_species = length(params["electrolyte"]["ionic_species"])
    num_species = num_ionic_species + 1
    # Store names of the columns and rows of the species
    species_col_names = [property_name for property_name in
                         AxisArrays.axes(processed_params_distr.species,2)]
    species_row_names = vcat("solvent",
                        [species["name"] for species in
                        params["electrolyte"]["ionic_species"]])

    # Convert species property values to measurements required to carry out the error
    # propagation to the parameters used in the model
    dict_solvent_params = params["electrolyte"]["solvent"]
    dict_ionic_species_params = params["electrolyte"]["ionic_species"]
    solvent_col_names_measurements = ["partial_molar_area", "partial_molar_volume",
                                      "adsorption_energy"]
    num_solvent_col_names_measurements = length(solvent_col_names_measurements)
    solvent_measurements = AxisArray(Vector{Measurement}(
        undef, num_solvent_col_names_measurements);
        col=solvent_col_names_measurements)
    # Set measurements of solvent and ionic species parameters
    [solvent_measurements[col=col_name] = convert_to_measurement(
        dict_solvent_params[col_name]) for col_name in solvent_col_names_measurements]

    ionic_species_col_names_measurements = ["partial_molar_area_central_ion",
                                            "partial_molar_volume_central_ion",
                                            "adsorption_energy_central_ion",
                                            "solvation_number", "surface_solvation_number",
                                            "concentration_bulk"]
    ionic_species_measurements = AxisArray(Matrix{Measurement}(undef,
            num_ionic_species, length(ionic_species_col_names_measurements));
            row=species_row_names[2:end], col=ionic_species_col_names_measurements)
    [ionic_species_measurements[row=idx_row, col=col_name] = convert_to_measurement(
        dict_ionic_species_params[idx_row][col_name]) for
        idx_row in 1:num_ionic_species,
        col_name in ionic_species_col_names_measurements]


    # Resize species matrix
    processed_params_distr.species = AxisArray(Matrix{ParameterDistribution}(
        undef, num_species, length(species_col_names));
        row=species_row_names, col=species_col_names)

    # Set distributions of solvent parameters
    solvent_col_names = ["charge_number",
                        "partial_molar_volume",
                        "partial_molar_area", "adsorption_energy"]
    ionic_species_col_names = ["charge_number", "concentration_bulk"]
    [processed_params_distr.species[row=1, col=col_name] = convert_to_normal(
        dict_solvent_params[col_name]) for col_name in solvent_col_names]
    [processed_params_distr.species[row=idx_row+1, col=col_name] = convert_to_normal(
        dict_ionic_species_params[idx_row][col_name]) for
        idx_row in 1:num_ionic_species, col_name in ionic_species_col_names]

    calculated_ionic_species_col_names = ["partial_molar_volume",
    "partial_molar_area", "adsorption_energy"]
    calculated_ionic_species_measurements = AxisArray(Matrix{Measurement}(
        undef, num_ionic_species, length(calculated_ionic_species_col_names));
        row=species_row_names[2:end], col=calculated_ionic_species_col_names)

    for (idx, species) in enumerate(params["electrolyte"]["ionic_species"])

        # evaluate statistics of the partial molar volume of species
        calculated_ionic_species_measurements[row=idx,col="partial_molar_volume"] =
            ionic_species_measurements[row=idx, col="partial_molar_volume_central_ion"] +
            ionic_species_measurements[row=idx, col="solvation_number"] *
            solvent_measurements[col="partial_molar_volume"]

        processed_params_distr.species[row=idx+1, col="partial_molar_volume"] =
        convert_to_normal(calculated_ionic_species_measurements[row=idx,col="partial_molar_volume"])

        # evaluate statistics of the molar area of species
        calculated_ionic_species_measurements[row=idx,col="partial_molar_area"] =
            ionic_species_measurements[row=idx, col="partial_molar_area_central_ion"] +
            ionic_species_measurements[row=idx, col="surface_solvation_number"] *
            solvent_measurements[col="partial_molar_area"]

        # update axis_array_species_stochastic
        processed_params_distr.species[row=idx+1, col="partial_molar_area"] =
            convert_to_normal(calculated_ionic_species_measurements[row=idx,col="partial_molar_area"])

        # evaluate statistics of the species adsorption energy
        calculated_ionic_species_measurements[row=idx,col="adsorption_energy"] =
            ionic_species_measurements[row=idx, col="adsorption_energy_central_ion"] +
            ionic_species_measurements[row=idx, col="surface_solvation_number"] *
            solvent_measurements[col="adsorption_energy"]

        # update axis_array_species_stochastic
        processed_params_distr.species[row=idx+1, col="adsorption_energy"] =
            convert_to_normal(calculated_ionic_species_measurements[row=idx,col="adsorption_energy"])
    end

    m_electrode_molar_area = convert_to_measurement(params["electrode"]["molar_area"])
    processed_params_distr.electrode["molar_area"] =
    convert_to_normal(m_electrode_molar_area / params["electrode"]["adsorption_sites_per_metal_ion"]) # FIXME: Molar area vacancy?

    solvent_concentration_bulk = (measurement(1.0) -
        sum(ionic_species_measurements[row=1:end, col="concentration_bulk"].*
        (calculated_ionic_species_measurements[row=1:end,col="partial_molar_volume"])))/
        (solvent_measurements[col="partial_molar_volume"])
    processed_params_distr.species[row=1, col="concentration_bulk"] = convert_to_normal(solvent_concentration_bulk)

    total_concentration_bulk = sum(ionic_species_measurements[col="concentration_bulk"]) + solvent_concentration_bulk
    processed_params_distr.species[row=1,col="molar_fraction_bulk"] = convert_to_normal(solvent_concentration_bulk/total_concentration_bulk)
    processed_params_distr.species[row=2:end,col="molar_fraction_bulk"] = convert_to_normal.(ionic_species_measurements[col="concentration_bulk"]/total_concentration_bulk)

    return processed_params_distr
end


function preprocess_reaction_parameters(params)

    ParameterDistribution = Normal{Float64} # Truncated{Normal{Float64},Continuous,Float64}
    processed_params_distr = ReactionParameters{ParameterDistribution}()

    reaction_params = params["heterogeneous_reactions"][1]
    processed_params_distr.reactions[row=1,col="k0"] =
        convert_to_normal(reaction_params["reaction_constant"])
    processed_params_distr.reactions[row=1,col="βs"] =
        convert_to_normal(reaction_params["coeff_beta"])
    processed_params_distr.reactions[row=1,col="as"] =
        convert_to_normal(reaction_params["coeff_a"])
    processed_params_distr.reactions[row=1,col="charge"] =
        convert_to_normal(reaction_params["charge"])

    return processed_params_distr
end
