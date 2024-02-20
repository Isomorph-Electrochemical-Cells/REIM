@with_kw mutable struct ModelParameters{T1, T2, T3}
    electrode::AxisArray{T1, 1, Vector{T1}, Tuple{AxisArrays.Axis{:col, Vector{String}}}} =
        AxisArray(Vector{T1}(undef,3);
            col=["charge_number_metal_ion",
            "molar_area",
            "adsorption_sites_per_metal_ion"])
    electrolyte::AxisArray{T2, 1, Vector{T2}, Tuple{AxisArrays.Axis{:col, Vector{String}}}} =
        AxisArray(Vector{T2}(undef,2);
            col=["susceptibility", "pressure_bulk"])
    species::AxisArray{T3, 2, Matrix{T3},
        Tuple{AxisArrays.Axis{:row, Vector{String}}, AxisArrays.Axis{:col, Vector{String}}}} =
            AxisArray(Matrix{T3}(undef,1,6);
                    row=["solvent"],
                    col=["charge_number",
                        "partial_molar_area", "partial_molar_volume", "adsorption_energy",
                        "concentration_bulk", "molar_fraction_bulk"])
end


@with_kw mutable struct ReactionParameters{T1}
    reactions::AxisArray{T1, 2, Matrix{T1},
        Tuple{AxisArrays.Axis{:row, Vector{String}}, AxisArrays.Axis{:col, Vector{String}}}} =
            AxisArray(Matrix{T1}(undef,1,4);
                    row=["0"],
                    col=["k0", "βs", "as", "charge"])
end

function check_params(params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    if (params.electrode[col="molar_area"] <= 0.0 ||
        params.electrode[col="adsorption_sites_per_metal_ion"] <= 0.0 ||
        params.electrolyte[col="susceptibility"] <= 0.0 ||
        any(params.species[col="partial_molar_volume"].<=0.0) ||
        any(params.species[col="partial_molar_area"].<=0.0) ||
        any(params.species[col="concentration_bulk"].<=0.0) ||
        any(params.species[col="molar_fraction_bulk"].<=0.0))
        return false
    else
        return true
    end
end


function check_reaction_params(params::ReactionParameters{T1}) where {T1}
    if (any(params.reactions[col="k0"] .<= 0.0) ||
        any(params.reactions[col="βs"] .<= 0.0) ||
        any(params.reactions[col="as"] .<= 0.0))
        return false
    else
        return true
    end
end

function total_concentration(Δϕ, Δp, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    partial_molar_volume = @view params.species[col="partial_molar_volume"]
    one(Δϕ)/(sum(partial_molar_volume.*molar_fraction(Δϕ, Δp, params)))
end

function total_concentration_nd(Δϕ, Δp, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    partial_molar_volume = @view params.species[col="partial_molar_volume"]
    one(Δϕ)/(sum(partial_molar_volume.*molar_fraction_nd(Δϕ, Δp, params)))
end

function molar_fraction(idx_species, Δϕ, Δp, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    species = @view params.species[row=idx_species]
    return species[col="molar_fraction_bulk"]*exp(-species[col="charge_number"]*f0*Δϕ - species[col="partial_molar_volume"]/(RT0)*Δp)
end

function molar_fraction_nd(idx_species, Δϕ, Δp, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    species = @view params.species[row=idx_species]
    return species[col="molar_fraction_bulk"]*exp(-species[col="charge_number"]*Δϕ - species[col="partial_molar_volume"]*Δp)
end


function molar_fraction(Δϕ, Δp, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_fraction_bulk = @view params.species[col="molar_fraction_bulk"]
    charge_number = @view params.species[col="charge_number"]
    partial_molar_volume = @view params.species[col="partial_molar_volume"]

    return molar_fraction_bulk.*exp.(-charge_number*f0.*Δϕ .- partial_molar_volume/(RT0).*Δp)
end

function molar_fraction_nd(Δϕ, Δp, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_fraction_bulk = @view params.species[col="molar_fraction_bulk"]
    charge_number = @view params.species[col="charge_number"]
    partial_molar_volume = @view params.species[col="partial_molar_volume"]

    return molar_fraction_bulk.*exp.(-charge_number.*Δϕ .- partial_molar_volume.*Δp)
end

function free_charge_density(Δϕ, Δp, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    charge_number = @view species.species[col="charge_number"]
    return F*total_concentration(Δϕ, Δp, params)*sum(molar_fraction(idx_species, Δϕ, Δp, params).*charge_number)
end

function free_charge_density_nd(Δϕ, Δp, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    charge_number = @view species.species[col="charge_number"]
    return total_concentration(Δϕ, Δp, params)*sum(molar_fraction_nd(idx_species, Δϕ, Δp, params).*charge_number)
end

function stored_charge_boundary_layer(Δϕ, Δp, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    χ = params.electrolyte[col="susceptibility"]
    sign(Δϕ)*sqrt(2.0*ϵ0*(1.0+χ)*Δp)
end

function stored_charge_boundary_layer_nd(Δϕ, Δp, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    χ = params.electrolyte[col="susceptibility"]
    sign(Δϕ)*sqrt(2.0*(1.0+χ)*Δp) #*sqrt_p0_eps0_div_q0;
end

function volumetric_concentration_constraint(Δp, tuple_Δϕ_params)
    (Δϕ, params) = tuple_Δϕ_params
    sum(molar_fraction(Δϕ, Δp, params)) - one(typeof(Δp))
end

function volumetric_concentration_constraint_nd(Δp, tuple_Δϕ_params)
    (Δϕ, params) = tuple_Δϕ_params
    sum(molar_fraction_nd(Δϕ, Δp, params)) - one(typeof(Δp))
end

function pressure_difference(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    zero_problem = ZeroProblem(volumetric_concentration_constraint, zero(eltype(params.species)))
    solve(zero_problem, Order1(), p=(Δϕ, params), maxevals=2000)
end

function pressure_difference_nd(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    zero_problem = ZeroProblem(volumetric_concentration_constraint_nd, zero(eltype(params.species)))
    solve(zero_problem, Order1(), p=(Δϕ, params), maxevals=2000)
end

# Evaluate capacitance resulting from stored charges in the boundary layer
function capacitance_bl(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    Δp = pressure_difference(Δϕ, params)

    charge_number = @view params.species[col="charge_number"]
    partial_molar_volume = @view params.species[col="partial_molar_volume"]

    χ = params.electrolyte[col="susceptibility"]
    y_species = molar_fraction(Δϕ, Δp, params)
    ∂_Δp_∂Δϕ = -F*sum(charge_number.*y_species) / sum(partial_molar_volume.*y_species)

    sqrt.(ϵ0*(1.0+χ))*sign(Δϕ)*∂_Δp_∂Δϕ / sqrt(max(2.0*Δp, eps(Δp)))
end

# Evaluate capacitance resulting from stored charges in the boundary layer
function capacitance_bl_nd(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    Δp = pressure_difference_nd(Δϕ, params)

    charge_number = @view params.species[col="charge_number"]
    partial_molar_volume = @view params.species[col="partial_molar_volume"]

    χ = params.electrolyte[col="susceptibility"]
    y_species = molar_fraction_nd(Δϕ, Δp, params)
    ∂_Δp_∂Δϕ = -dot(charge_number, y_species) / dot(partial_molar_volume, y_species)

    sqrt.(1.0+χ)*sign(Δϕ)*∂_Δp_∂Δϕ / sqrt(max(2.0*Δp, eps(Δp)))#*sqrt_p0_eps0_div_q0
end

function ionic_strength(species_concentrations, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    charge_number = @view params.electrolyte["species"][col="charge_number"]
    0.5*sum((charge_number.^2).*species_concentrations)
end

function ionic_strength_nd(species_concentrations, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    charge_number = @view params.electrolyte["species"][col="charge_number"]
    0.5*sum((charge_number.^2).*species_concentrations)
end

function debye_length(ionic_strength_value, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3} # CHECK THIS!
    ϵr = one(eltype(params.electrolyte)) + params.electrolyte[col="susceptibility"]
    sqrt(ϵ0*ϵr*/(2.0*F*f0*ionic_strength_value))
end

function debye_length_nd(ionic_strength_value, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3} # CHECK THIS!
    ϵr = one(eltype(params.electrolyte)) + params.electrolyte[col="susceptibility"]
    sqrt(ϵ0*ϵr*/(2.0*F*f0*l0^2*c0*ionic_strength_value))
end

function molar_fraction_surface_vacancy(Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    partial_molar_area = params.electrode[col="molar_area"]
    ωm = params.electrode[col="adsorption_sites_per_metal_ion"]
    exp.((partial_molar_area/(ωm*RT0)) * Δγ)
end

function molar_fraction_surface_vacancy_nd(Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    partial_molar_area = params.electrode[col="molar_area"]
    ωm = params.electrode[col="adsorption_sites_per_metal_ion"]
    exp.((partial_molar_area/(ωm)) * Δγ)
end

function molar_fraction_surface_species(idx_species, Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    y_bulk = @view params.species[row=idx_species, col="molar_fraction_bulk"]
    charge_number = @view params.species[row=idx_species, col="charge_number"]
    adsorption_energy = @view params.species[row=idx_species, col="adsorption_energy"]
    partial_molar_area = @view params.species[row=idx_species, col="partial_molar_area"]

    y_bulk.*exp(-adsorption_energy*eV_div_kb_T0 .- charge_number*f0.*Δϕ .+
        partial_molar_area./(RT0).*Δγ)
end

function molar_fraction_surface_species_nd(idx_species, Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    y_bulk = @view params.species[row=idx_species, col="molar_fraction_bulk"]
    charge_number = @view params.species[row=idx_species, col="charge_number"]
    adsorption_energy = @view params.species[row=idx_species, col="adsorption_energy"]
    partial_molar_area = @view params.species[row=idx_species, col="partial_molar_area"]

    y_bulk.*exp(-adsorption_energy .- charge_number.*Δϕ .+ partial_molar_area.*Δγ)
end

function molar_fraction_surface_species(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    y_bulk = @view params.species[col="molar_fraction_bulk"]
    charge_number = @view params.species[col="charge_number"]
    adsorption_energy = @view params.species[col="adsorption_energy"]
    partial_molar_area = @view params.species[col="partial_molar_area"]

    y_bulk.*exp.(-adsorption_energy*eV/(kbT0) .- charge_number*f0.*Δϕ .+
        partial_molar_area./(RT0).*Δγ)
end

function molar_fraction_surface_species_nd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    y_bulk = @view params.species[col="molar_fraction_bulk"]
    charge_number = @view params.species[col="charge_number"]
    adsorption_energy = @view params.species[col="adsorption_energy"]
    partial_molar_area = @view params.species[col="partial_molar_area"]

    y_bulk.*exp.(-adsorption_energy .- charge_number.*Δϕ .+ partial_molar_area.*Δγ)
end

function surface_concentration_constraint(Δγ, tuple_Δϕ_params)
    (Δϕ, params) = tuple_Δϕ_params
    molar_fraction_surface_vacancy(Δγ, params) .+
        sum(molar_fraction_surface_species(Δϕ, Δγ, params)) .- one(Δγ)
end

function surface_concentration_constraint_nd(Δγ, tuple_Δϕ_params)
    (Δϕ, params) = tuple_Δϕ_params
    molar_fraction_surface_vacancy_nd(Δγ, params) .+
        sum(molar_fraction_surface_species_nd(Δϕ, Δγ, params)) .- one(Δγ)
end


function surface_tension_difference(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    zero_problem = ZeroProblem(surface_concentration_constraint, zero(eltype(params.species)))
    return solve(zero_problem, Order1(), p=(Δϕ,params), maxevals=2000)
end

function surface_tension_difference_nd(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    zero_problem = ZeroProblem(surface_concentration_constraint_nd, zero(eltype(params.species)))
    return solve(zero_problem, Order1(), p=(Δϕ,params), maxevals=2000)
end

function total_surface_concentration(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_area_surface = @view params.species[col="partial_molar_area"]
    molar_area_vacancy = params.electrode["molar_area"] # FIXME: molar_area_vacancy?
    y_species_surface = molar_fraction_surface_species(Δϕ, Δγ, params)
    y_surface_vacancy = molar_fraction_surface_vacancy(Δγ, params)

    one(T1)/(sum(molar_area_surface.*y_species_surface) + molar_area_vacancy*y_surface_vacancy)
end

function total_surface_concentration_nd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_area_surface = @view params.species[col="partial_molar_area"]
    molar_area_vacancy = params.electrode["molar_area"] # FIXME: molar_area_vacancy?
    y_species_surface = molar_fraction_surface_species_nd(Δϕ, Δγ, params)
    y_surface_vacancy = molar_fraction_surface_vacancy_nd(Δγ, params)

    one(T1)/(sum(molar_area_surface.*y_species_surface) + molar_area_vacancy*y_surface_vacancy)
end

function surface_coverage_species(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_area_surface = @view params.species[col="partial_molar_area"]
    y_species_surface = molar_fraction_surface_species(Δϕ, Δγ, params)

    y_species_surface.*molar_area_surface*total_surface_concentration(Δϕ, Δγ, params)
end

function surface_coverage_species_nd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_area_surface = @view params.species[col="partial_molar_area"]
    y_species_surface = molar_fraction_surface_species_nd(Δϕ, Δγ, params)

    y_species_surface.*molar_area_surface*total_surface_concentration_nd(Δϕ, Δγ, params)
end

function surface_coverage_species(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_area_surface = @view params.species[col="partial_molar_area"]
    Δγ = surface_tension_difference(Δϕ, params)
    y_species_surface = molar_fraction_surface_species(Δϕ, Δγ, params)

    y_species_surface.*molar_area_surface*total_surface_concentration(Δϕ, Δγ, params)
end

function surface_coverage_species_nd(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_area_surface = @view params.species[col="partial_molar_area"]
    Δγ = surface_tension_difference_nd(Δϕ, params)
    y_species_surface = molar_fraction_surface_species_nd(Δϕ, Δγ, params)

    y_species_surface.*molar_area_surface*total_surface_concentration_nd(Δϕ, Δγ, params)
end

function surface_coverage_vacancy(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_area_vacancy = @view params.electrode[col="molar_area"] # molar area vacancy?
    molar_fraction_surface_vacancy(Δγ, params)*total_surface_concentration(Δϕ, Δγ, params)*molar_area_vacancy
end

function surface_coverage_vacancy_nd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_area_vacancy = @view params.electrode[col="molar_area"] # molar area vacancy?
    molar_fraction_surface_vacancy_nd(Δγ, params)*total_surface_concentration_nd(Δϕ, Δγ, params)*molar_area_vacancy
end

function surface_coverage_vacancy(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_area_vacancy = params.electrode["molar_area"] # molar area vacancy?
    Δγ = surface_tension_difference(Δϕ, params)
    molar_fraction_surface_vacancy(Δγ, params)*total_surface_concentration(Δϕ, Δγ, params)*molar_area_vacancy
end

function surface_coverage_vacancy_nd(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    molar_area_vacancy = params.electrode["molar_area"] # molar area vacancy?
    Δγ = surface_tension_difference_nd(Δϕ, params)
    molar_fraction_surface_vacancy_nd(Δγ, params)*total_surface_concentration_nd(Δϕ, Δγ, params)*molar_area_vacancy
end

function surface_terms_nd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    y_species_surface = molar_fraction_surface_species_nd(Δϕ, Δγ, params)
    charge_numbers = @view params.species[col="charge_number"]
    molar_area_species = @view params.species[col="partial_molar_area"]
    y_vacancy = molar_fraction_surface_vacancy_nd(Δγ, params::ModelParameters{T1,T2,T3})
    molar_area_vacancy = params.electrode[col="molar_area"]

    f1s = dot(charge_numbers, y_species_surface)
    f2s = molar_area_vacancy * y_vacancy + dot(molar_area_species, y_species_surface)
    f3s = sum(charge_numbers.^2 .* y_species_surface)
    f4s = sum(molar_area_species .* y_species_surface .* charge_numbers)
    f5s = molar_area_vacancy^2 * y_vacancy + sum(molar_area_species.^2 .* y_species_surface)

    return @SVector [f1s, f2s, f3s, f4s, f5s]
end

function qs_nd(surface_terms)
    -surface_terms[1] / surface_terms[2]
end

function qs_nd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    y_species_surface = molar_fraction_surface_species_nd(Δϕ, Δγ, params)
    charge_numbers = @view params.species[col="charge_number"]
    molar_area_species = @view params.species[col="partial_molar_area"]
    y_vacancy = molar_fraction_surface_vacancy_nd(Δγ, params::ModelParameters{T1,T2,T3})
    molar_area_vacancy = params.electrode[col="molar_area"]

    -dot(charge_numbers, y_species_surface) /
        (molar_area_vacancy * y_vacancy + dot(molar_area_species, y_species_surface))
end

function ∂qs_∂Δϕ_nd(surface_terms)
    (surface_terms[3] * surface_terms[2] - surface_terms[1] * surface_terms[4]) /
        surface_terms[2]^2
end

function ∂qs_∂Δϕ_nd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    y_species_surface = molar_fraction_surface_species_nd(Δϕ, Δγ, params)
    charge_numbers = @view params.species[col="charge_number"]
    molar_area_species = @view params.species[col="partial_molar_area"]
    y_vacancy = molar_fraction_surface_vacancy_nd(Δγ, params::ModelParameters{T1,T2,T3})
    molar_area_vacancy = params.electrode[col="molar_area"]

    (sum(charge_numbers.^2 .* y_species_surface) * (molar_area_vacancy * y_vacancy +
    dot(molar_area_species, y_species_surface)) - dot(charge_numbers, y_species_surface) *
    sum(molar_area_species .* y_species_surface .* charge_numbers) ) /
        (molar_area_vacancy * y_vacancy + dot(molar_area_species, y_species_surface))^2
end

function ∂qs_∂Δγ_nd(surface_terms)
    (-surface_terms[4] * surface_terms[2] + surface_terms[1] * surface_terms[5]) /
        surface_terms[2]^2
end

function ∂qs_∂Δγ_nd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    y_species_surface = molar_fraction_surface_species_nd(Δϕ, Δγ, params)
    charge_numbers = @view params.species[col="charge_number"]
    molar_area_species = @view params.species[col="partial_molar_area"]
    y_vacancy = molar_fraction_surface_vacancy_nd(Δγ, params::ModelParameters{T1,T2,T3})
    molar_area_vacancy = params.electrode[col="molar_area"]

    (-sum(charge_numbers .* molar_area_species .* y_species_surface) * (molar_area_vacancy * y_vacancy + dot(molar_area_species, y_species_surface)) + dot(charge_numbers, y_species_surface) * (molar_area_vacancy^2 * y_vacancy + sum(molar_area_species.^2 .* y_species_surface))) / (molar_area_vacancy * y_vacancy + dot(molar_area_species, y_species_surface))^2
end

function dqs_dΔγ_nd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    #∂qs_∂Δϕ_nd(Δϕ, Δγ, params) - qs_nd(Δϕ, Δγ, params) * ∂qs_∂Δγ_nd(Δϕ, Δγ, params)
    surface_terms = surface_terms_nd(Δϕ, Δγ, params)
    ∂qs_∂Δϕ_nd(surface_terms) - qs_nd(surface_terms) * ∂qs_∂Δγ_nd(surface_terms)
end

function ∂qs_∂Δϕ_nd_fd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    dΔϕ = sqrt(eps(Δϕ))*abs(Δϕ)
    (qs_nd(Δϕ + dΔϕ, Δγ, params) - qs_nd(Δϕ - dΔϕ, Δγ, params)) / (2*dΔϕ)
end

function ∂qs_∂Δγ_nd_fd(Δϕ, Δγ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    dΔγ = sqrt(eps(Δγ))*abs(Δγ)
    (qs_nd(Δϕ, Δγ + dΔγ, params) - qs_nd(Δϕ, Δγ - dΔγ, params)) / (2*dΔγ)
end

function capacitance_surface(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    dΔϕ = 1e+2*sqrt(eps(Δϕ))
    Δγ_l = surface_tension_difference(Δϕ - dΔϕ, params)
    Δγ_c = surface_tension_difference(Δϕ, params)
    Δγ_r = surface_tension_difference(Δϕ + dΔϕ, params)

    ∂2_Δγ_∂x2 = (Δγ_l - 2.0*Δγ_c + Δγ_r) / (dΔϕ^2)
    return -∂2_Δγ_∂x2
end

function capacitance_surface_nd(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    Δγ = surface_tension_difference_nd(Δϕ, params)
    dqs_dΔγ_nd(Δϕ, Δγ, params)
end

function capacitance(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    return capacitance_bl(Δϕ, params) + capacitance_surface(Δϕ, params)
end

function capacitance_nd(Δϕ, params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    return capacitance_bl_nd(Δϕ, params) + capacitance_surface_nd(Δϕ, params)
end

function init_model!(params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    # re-evaluate solvent concentration and molar fractions
    bulk_solvent_concentration = one(T1) -
        (sum(params.species[col="concentration_bulk"][2:end] .*
        params.species[col="partial_molar_volume"][2:end])) /
        params.species[col="partial_molar_volume"][1]
    params.species[col="concentration_bulk"][1] = bulk_solvent_concentration
    total_concentration = sum(params.species[col="concentration_bulk"])
    params.species[col="molar_fraction_bulk"] =
        params.species[col="concentration_bulk"]./total_concentration
end

function init_model_nd!(params::ModelParameters{T1,T2,T3}) where {T1,T2,T3}
    # re-evaluate solvent concentration and molar fractions
    bulk_solvent_concentration = one(T1) -
        (sum(params.species[col="concentration_bulk"][2:end] .*
        params.species[col="partial_molar_volume"][2:end])) /
        params.species[col="partial_molar_volume"][1]
    params.species[col="concentration_bulk"][1] = bulk_solvent_concentration
    total_concentration = sum(params.species[col="concentration_bulk"])
    params.species[col="molar_fraction_bulk"] =
        params.species[col="concentration_bulk"]./total_concentration
end

# The following function evaluates the reaction rate for the simple
# redox reaction O + n e ⇋ R with no exlusive surface species.
# The overpotential is defined with respect to the electrolyte concentrations
# just outside the double layer
function reaction_redox_no_exlusive_surface_species(η, n, params)
    βs = params.reactions[1][col="βs"]
    As = params.reactions[1][col="As"]

    # forward reaction rate
    rf = rs0 * exp(βs*As*n*f0*η)
    # backward reaction rate
    rb = rs0 * exp((1-βs)*As*n*f0*η)
    # overall reaction rate
    return rf - rb
end


function current_electrode_plating(ηs, nc, nc_eq, params)
    rp = params.reactions[row=1]
    k0 = rp[col="k0"]/1e+5
    βs = rp[col="βs"]
    as = rp[col="as"]
    zo = rp[col="charge"]
    current =  zo * F * k0 * ((nc / nc_eq)^βs * exp(-βs*as*zo*f0*ηs) -
        (nc_eq / nc)^(1-βs) * exp((1-βs)*as*zo*f0*ηs))
    return current
end
