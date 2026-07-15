# ---------------------------------------------------------------------------
# Partition functions.
#
# The internal partition function (Z) for a given element and ionization
# stage is computed by summing over all known energy levels:
#
#     Z(T) = Σ_i g_i · exp(−E_i / k_B T)
#
# where
#   g_i  = statistical weight (degeneracy) of level i
#   E_i  = energy of level i (cm⁻¹)
#   k_B T = thermal energy in cm⁻¹  (converted from eV: 1 eV = 8065.54 cm⁻¹)
# ---------------------------------------------------------------------------

"""
    partition_function(db, element, spectr_charge, temp)

Compute the internal partition function (electronic) for a given element
and ionization stage at temperature `temp`.

`temp` can be a `Quantity` in eV or K, or a bare number (assumed eV).

Returns Σ g·exp(−E / k_B T) over all NIST energy levels with valid
energy (≥ 0) and degeneracy (> 0).  Returns 0.0 if the element or
charge stage has no level data.
"""
function partition_function(db::LIBSDB, element::AbstractString, spectr_charge::Int, temp)
    temp_eV = to_internal_temp(temp)
    tbl = get(db.levels, element, nothing)
    tbl === nothing && return 0.0

    # Convert eV → cm⁻¹ for Boltzmann factor: 1 eV = 8065.54393734921 cm⁻¹
    temp_cm = temp_eV * 8065.54393734921

    z = 0.0
    for row in Tables.rows(tbl)
        # Only include levels with valid energy and degeneracy
        if row.spectr_charge == spectr_charge && row.energy >= 0 && row.g > 0
            z += row.g * exp(-row.energy / temp_cm)
        end
    end
    z
end
