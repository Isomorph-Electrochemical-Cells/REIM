# Define physical constants
const R  = 8.31446261815324 # Ideal gas constant [J/(K mol)]
const F  = 96485.33212 # Faraday constant [C/mol]*)
const ϵ0 = 8.8541878176*10^-12 # Permittivity of free space [F/m=C/(V*m)]
const nA = 6.02214076*10^23 # Avogadro constant 1/Mol
const eV = 1.602176634*10^-19 # 1 electron volt in joule
const kb = 1.380649*10^-23 # Boltzmann constant [J/K]
const e0 = F/nA # Electron charge [C]

# Fix independent reference constants
const T0 = 273.15+25.0 # Constant reference temperature [K]
const c0 = 1000 # reference concentration [mol / m^3]
const ϵr0 = 2.0 # reference relative permittivity
const l0 = sqrt((ϵ0*ϵr0*R*T0)/(2*F^2*c0)) # reference Debye length [m]
const D0 = 1e-9 # reference diffusion coefficient [m^2 / s]
const V0 = R*T0/F # Reference voltage [V]

# Derived reference constants
const p0 = c0*R*T0 # reference pressure [J / m^3]
const cF0 = F*c0 # reference volumetric charge density [Coulomb / m^3]
const q0 = l0*cF0 # reference charge stored in double layer [Coulomb / m^2]
const E0 = p0*l0^3 # reference energy [J]
const γ0 = E0/l0^2 # reference surface tension [J/m^2]
const ν0 = 1/c0 # reference molar volume [m^3 / mol]
const a0 = 1/(c0*l0) # reference molar area [m^2 / mol]
const C0 = q0/V0 # reference capacitance [Coulomb / m^2 / V]

# Useful often-used combinations of constants to speed-up the computations
const f0 = F/(R*T0) # Reference inverse thermal voltage [1/V]
const eV_div_kb_T0 = eV/(kb*T0)
const RT0 = R*T0 # Thermal energy [J/mol]
const kbT0 = kb*T0 # Thermal energy [J]
const sqrt_p0_eps0_div_q0 = sqrt(p0*ϵ0)/q0 # equivalent to V0/l0 * sqrt(ϵ0/p0) = \sqrt(2/ϵr0) = 1

# Useful relationships: (R*T)/nA= kb*T, nA=R/kb, f0 = R/(R*T0) = e0/(kb*T0)
