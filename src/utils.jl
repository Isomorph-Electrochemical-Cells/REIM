
function nondimensionalize!(params)
    params.electrode[col="molar_area"] /= a0
    params.electrolyte[col="pressure_bulk"] /= p0

    species = params.species
    species[col="partial_molar_area"] /= a0
    species[col="partial_molar_volume"] /= ν0
    species[col="adsorption_energy"] *= eV_div_kb_T0
    species[col="concentration_bulk"] /= c0

    return params
end

function dimenisonalize!(params)
    params.electrode[col="molar_area"] *= a0
    params.electrolyte[col="pressure_bulk"] *= p0

    species = params.species
    species[col="partial_molar_area"] *= a0
    species[col="partial_molar_volume"] *= ν0
    species[col="adsorption_energy"] /= eV_div_kb_T0
    species[col="concentration_bulk"] *= c0

    return params
end
