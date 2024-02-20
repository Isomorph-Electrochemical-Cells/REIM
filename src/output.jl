function distr_to_number_or_dict(distr)
    distr_mean = mean(distr)
    distr_std = std(distr)

    if distr_std <= (abs(distr_mean)+eps(distr_mean))/1e+3
        return distr_mean
    else
        dict = Dict()
        dict["mean"] = mean(distr)
        dict["std"] = std(distr)
        return dict
    end
end

function parameter_distr_to_dict(parameter_distr)
    dict = Dict()
    dict["electrode"] = Dict()
    for key in AxisArrays.axes(parameter_distr.electrode, 1)
        dict["electrode"][key] = distr_to_number_or_dict(parameter_distr.electrode[col=key])
    end
    dict["electrolyte"] = Dict()
    for key in AxisArrays.axes(parameter_distr.electrolyte, 1)
        dict["electrolyte"][key] = distr_to_number_or_dict(parameter_distr.electrolyte[col=key])
    end
    solvent_dict = Dict()
    for species_param in AxisArrays.axes(parameter_distr.species, 2)
        solvent_dict[species_param] = distr_to_number_or_dict(parameter_distr.species[row=1,col=species_param])
    end
    dict["electrolyte"]["solvent"] = solvent_dict

    species = dict["electrolyte"]["ionic_species"] = []
    ionic_species_names = AxisArrays.axes(parameter_distr.species, 1)[2:end]
    for ionic_species_name in ionic_species_names
        dict_species = Dict()
        dict_species["name"] = ionic_species_name
        for species_param in AxisArrays.axes(parameter_distr.species, 2)
            dict_species[species_param] = distr_to_number_or_dict(parameter_distr.species[row=ionic_species_name,col=species_param])
        end
        push!(species, dict_species)
    end
    dict["electrolyte"]["species"] = species
    return dict
end

function write_parameter_distr_to_file(parameter_distr, file_path)

    # convert model parameter distributions to a dictionary
    parameter_dict = parameter_distr_to_dict(parameter_distr)

    indentation = 4 # pretty-print the JSON content with a specific indentation
    touch(file_path)
    open(file_path, "w") do file
        JSON.print(file, parameter_dict, indentation)
    end
end
