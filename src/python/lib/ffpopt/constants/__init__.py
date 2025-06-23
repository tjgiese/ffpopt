#!/usr/bin/env python3

"""
Conversion factors and table-lookups

Brief summary of functions
--------------------------
GetAtomicNumber(ele) -> int
    Get atomic number from atomic symbol

GetAtomicSymbol(z) -> str
    Get atomic symbol from atomic number

GetAtomicMass(z) -> float
    Return the average naturally occuring mass in amu

GetStdIsotopeMass(z) -> float
    Return the standard isotope mass (light isotope) in amu

GetSecondIsotopeMass(z) -> float
    Return the second isotope mass (heavy isotope) in amu

GetCovalentRadius(z) -> float
    Return the covalent radius in Bohr

GetvdWRadius(z) -> float
    Return the vdW radius in Bohr

GetElectronegativity(z) -> float
    Return the electronegativity in atomic units (Hartree/electron)

GetHardness(z) -> float
    Return the hardness in atomic units (Hartree/electron)



"""

from . Conversions import PI
from . Conversions import PI_D_180
from . Conversions import PIO2
from . Conversions import TWO_PI
from . Conversions import FOUR_PI
from . Conversions import SQRT_PI
from . Conversions import SQRT2
from . Conversions import SQRT3
from . Conversions import EULER
from . Conversions import LN10
from . Conversions import RAD2DEG
from . Conversions import DEG2RAD
from . Conversions import AVOGADRO_CONSTANT
from . Conversions import SPEED_OF_LIGHT_SI
from . Conversions import BOLTZMANN_CONSTANT_SI
from . Conversions import PLANCK_CONSTANT_SI
from . Conversions import HBAR_SI
from . Conversions import ELECTRON_CHARGE_SI
from . Conversions import ELECTRON_REST_MASS_SI
from . Conversions import PROTON_REST_MASS_SI
from . Conversions import NEUTRON_REST_MASS_SI
from . Conversions import PERMITTIVITY_FREE_SPACE_SI
from . Conversions import FARADAY_CONSTANT_SI
from . Conversions import FINE_STRUCTURE_CONSTANT_SI
from . Conversions import RYDBERG_CONSTANT_SI
from . Conversions import AU_MASS_SI
from . Conversions import AU_CHARGE_SI
from . Conversions import AU_ANGULAR_MOMENTUM_SI
from . Conversions import AU_PERMITTIVITY_SI
from . Conversions import AU_DISTANCE_SI
from . Conversions import AU_ENERGY_SI
from . Conversions import AU_VELOCITY_SI
from . Conversions import AU_TIME_SI
from . Conversions import AU_ELECTROSTATIC_POTENTIAL_SI
from . Conversions import AU_FORCE_SI
from . Conversions import AU_MAGNETIC_DIPOLE_SI
from . Conversions import AU_PRESSURE_SI
from . Conversions import SPEED_OF_LIGHT_AU
from . Conversions import BOLTZMANN_CONSTANT_AU
from . Conversions import PLANCK_CONSTANT_AU
from . Conversions import PROTON_REST_MASS_AU
from . Conversions import NEUTRON_REST_MASS_AU
from . Conversions import RYDBERG_CONSTANT_AU
from . Conversions import AU_PER_PASCAL
from . Conversions import AU_PER_ATMOSPHERE
from . Conversions import AU_PER_ELECTRON_VOLT
from . Conversions import AU_PER_JOULE_PER_MOL
from . Conversions import AU_PER_KJOULE_PER_MOL
from . Conversions import AU_PER_KCAL_PER_MOL
from . Conversions import AU_PER_KCAL_PER_ANG_MOL
from . Conversions import AU_PER_INVERSE_CM
from . Conversions import AU_PER_ANGSTROM
from . Conversions import AU_PER_DEBYE
from . Conversions import AU_PER_GRAMS_PER_MOL
from . Conversions import AU_PER_CM3_PER_MOL
from . Conversions import AU_PER_ATOMIC_MASS_UNIT
from . Conversions import AU_PER_RYDBERG
from . Conversions import AU_PER_DYNE_PER_CM
from . Conversions import AU_PER_PKA_PER_KELVIN
from . Conversions import AU_PER_PKA_AT_300K
from . Conversions import AU_PER_PKA_AT_25C
from . Conversions import KT_AT_298K
from . Conversions import AU_PER_NM
from . Conversions import AU_PER_NM3
from . Conversions import AU_PER_ANG3
from . Conversions import AU_PER_GRAMS_PER_CM3
from . Conversions import AU_PER_CAL_PER_MOL
from . Conversions import AU_PER_BAR
from . Conversions import AU_PER_METER
from . Conversions import AU_PER_CM
from . Conversions import AU_PER_M3
from . Conversions import AU_PER_CM3
from . Conversions import AU_PER_FM2
from . PeriodicTable import GetAtomicNumber
from . PeriodicTable import GetAtomicSymbol
from . PeriodicTable import GetAtomicMass
from . PeriodicTable import GetStdIsotopeMass
from . PeriodicTable import GetSecondIsotopeMass
from . PeriodicTable import GetCovalentRadius
from . PeriodicTable import GetvdWRadius
from . PeriodicTable import GetElectronegativity
from . PeriodicTable import GetHardness


__all__ = ['PI',
           'PI_D_180',
           'PIO2',
           'TWO_PI',
           'FOUR_PI',
           'SQRT_PI',
           'SQRT2',
           'SQRT3',
           'EULER',
           'LN10',
           'RAD2DEG',
           'DEG2RAD',
           'AVOGADRO_CONSTANT',
           'SPEED_OF_LIGHT_SI',
           'BOLTZMANN_CONSTANT_SI',
           'PLANCK_CONSTANT_SI',
           'HBAR_SI',
           'ELECTRON_CHARGE_SI',
           'ELECTRON_REST_MASS_SI',
           'PROTON_REST_MASS_SI',
           'NEUTRON_REST_MASS_SI',
           'PERMITTIVITY_FREE_SPACE_SI',
           'FARADAY_CONSTANT_SI',
           'FINE_STRUCTURE_CONSTANT_SI',
           'RYDBERG_CONSTANT_SI',
           'AU_MASS_SI',
           'AU_CHARGE_SI',
           'AU_ANGULAR_MOMENTUM_SI',
           'AU_PERMITTIVITY_SI',
           'AU_DISTANCE_SI',
           'AU_ENERGY_SI',
           'AU_VELOCITY_SI',
           'AU_TIME_SI',
           'AU_ELECTROSTATIC_POTENTIAL_SI',
           'AU_FORCE_SI',
           'AU_MAGNETIC_DIPOLE_SI',
           'AU_PRESSURE_SI',
           'SPEED_OF_LIGHT_AU',
           'BOLTZMANN_CONSTANT_AU',
           'PLANCK_CONSTANT_AU',
           'PROTON_REST_MASS_AU',
           'NEUTRON_REST_MASS_AU',
           'RYDBERG_CONSTANT_AU',
           'AU_PER_PASCAL',
           'AU_PER_ATMOSPHERE',
           'AU_PER_ELECTRON_VOLT',
           'AU_PER_JOULE_PER_MOL',
           'AU_PER_KJOULE_PER_MOL',
           'AU_PER_KCAL_PER_MOL',
           'AU_PER_KCAL_PER_ANG_MOL',
           'AU_PER_INVERSE_CM',
           'AU_PER_ANGSTROM',
           'AU_PER_DEBYE',
           'AU_PER_GRAMS_PER_MOL',
           'AU_PER_CM3_PER_MOL',
           'AU_PER_ATOMIC_MASS_UNIT',
           'AU_PER_RYDBERG',
           'AU_PER_DYNE_PER_CM',
           'AU_PER_PKA_PER_KELVIN',
           'AU_PER_PKA_AT_300K',
           'AU_PER_PKA_AT_25C',
           'KT_AT_298K',
           'AU_PER_NM',
           'AU_PER_NM3',
           'AU_PER_ANG3',
           'AU_PER_GRAMS_PER_CM3',
           'AU_PER_CAL_PER_MOL',
           'AU_PER_BAR',
           'AU_PER_METER',
           'AU_PER_CM',
           'AU_PER_M3',
           'AU_PER_CM3',
           'AU_PER_FM2',
           'GetAtomicNumber',
           'GetAtomicSymbol',
           'GetAtomicMass',
           'GetStdIsotopeMass',
           'GetSecondIsotopeMass',
           'GetCovalentRadius',
           'GetvdWRadius',
           'GetElectronegativity',
           'GetHardness']
