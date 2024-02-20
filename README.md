# REIM (Reduced Electrochemical Interface Model)

The electrochemical interface region at an electrode-electrolyte is of great importance for surface processes, such as heterogeneous electrochemical reactions or sorption processes at the electrode surface. The model is based on a thermodynamic continuum-based metal-electrolyte interface theory developed in a series of papers [1,2]. This package allows for the determination of impedance and surface species coverage.

Model inputs:
- Bulk electrolyte composition
- Bulk electrolyte pressure
- Potential drop across the double layer

Model outputs:
- Double layer capacitance
- Surface fractions of adsorpted species

*Note*: This package is currently under development, so that the user interface is subject to change.

## License

The REIM.jl package is released under the MIT license (see https://opensource.org/license/MIT/ for additional information).

## Core developers

R. P. Schaerer

## Acknowledgements

This work is part of the SONAR project supported by the European Union / European Commission funding program Horizon 2020.

## References

[1] Landstorfer, M., C. Guhlke, and W. Dreyer. ‘Theory and Structure of the Metal-Electrolyte Interface Incorporating Adsorption and Solvation Effects’. Electrochimica Acta 201 (May 2016): 187–219. https://doi.org/10.1016/j.electacta.2016.03.013.

[2] Dreyer, Wolfgang, Clemens Guhlke, and Rüdiger Müller. ‘Bulk-Surface Electrothermodynamics and Applications to Electrochemistry’. Entropy 20, no. 12 (6 December 2018): 939. https://doi.org/10.3390/e20120939.
