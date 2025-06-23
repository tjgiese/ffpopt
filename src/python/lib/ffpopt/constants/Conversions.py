#!/usr/bin/env python3

def PI(): 
    '''Fundamental unitless constant'''
    return 3.141592653589793238462643383279502884197

def PI_D_180():
    '''Fundamental unitless constant'''
    return PI()/180.0


def PIO2():
    '''Fundamental unitless constant'''
    return PI() * 0.5


def TWO_PI():
    '''Fundamental unitless constant'''
    return 2.0 * PI()


def FOUR_PI():
    '''Fundamental unitless constant'''
    return 4.0 * PI()
    
    
def SQRT_PI():
    '''Fundamental unitless constant'''
    return 1.772453850905516027298167483341


def SQRT2():
    '''Fundamental unitless constant'''
    return 1.41421356237309504880168872420969807856967


def SQRT3():
    '''Fundamental unitless constant'''
    return 1.73205080756887729352744634150587236694281


def EULER():
    '''Fundamental unitless constant'''
    return 0.5772156649015328606065120900824024310422

def LN10():
    '''Fundamental unitless constant'''
    return 2.3025850929940459010936137929093092679977

def RAD2DEG():
    '''Converts a number in radians to a number in degrees'''
    return 180.0/PI()


def DEG2RAD():
    '''Converts a number from degrees to radians'''
    return PI()/180.0


def AVOGADRO_CONSTANT():
    '''Unitless'''
    return 6.02214179e+23


def SPEED_OF_LIGHT_SI():
    '''SI'''
    return 299792458.0


def BOLTZMANN_CONSTANT_SI():
    '''SI'''
    return 1.3806504e-23


def PLANCK_CONSTANT_SI():
    '''SI'''
    return 6.62606896e-34


def HBAR_SI():
    '''SI'''
    return PLANCK_CONSTANT_SI()/(2*PI())


def ELECTRON_CHARGE_SI():
    '''SI'''
    return 1.602176487e-19


def ELECTRON_REST_MASS_SI():
    '''SI'''
    return 9.10938215e-31


def PROTON_REST_MASS_SI():
    '''SI'''
    return 1.67262171e-27


def NEUTRON_REST_MASS_SI():
    ''' SI'''
    return 1.67492729e-27


def PERMITTIVITY_FREE_SPACE_SI():
  '''SI'''
  return 8.8541878176e-12



def FARADAY_CONSTANT_SI():
    '''SI  (derived)'''
    return AVOGADRO_CONSTANT()*ELECTRON_CHARGE_SI()



def FINE_STRUCTURE_CONSTANT_SI():
    ''' SI  (derived)'''
    return ELECTRON_CHARGE_SI()*ELECTRON_CHARGE_SI() / (2 * PERMITTIVITY_FREE_SPACE_SI() * PLANCK_CONSTANT_SI() * SPEED_OF_LIGHT_SI())



def RYDBERG_CONSTANT_SI():
    ''' SI  (derived)'''
    return FINE_STRUCTURE_CONSTANT_SI()*FINE_STRUCTURE_CONSTANT_SI() * ELECTRON_REST_MASS_SI() * SPEED_OF_LIGHT_SI() / (2 * PLANCK_CONSTANT_SI())




  #  ! REPLACED WITH CODATA FROM NIST
  #  !  FINE_STRUCTURE_CONSTANT = 7.2973525376E-3
  #  !  RYDBERG_CONSTANT 10973731.568527 ! 1/m
  #  !
  #  ! Atomic Units in terms of SI
  #  !
  #  ! Independent quantities
  #  !


def AU_MASS_SI():
    '''SI/AU'''
    return ELECTRON_REST_MASS_SI()



def AU_CHARGE_SI():
    '''SI/AU'''
    return ELECTRON_CHARGE_SI()


def AU_ANGULAR_MOMENTUM_SI():
    '''SI/AU'''
    return HBAR_SI()



def AU_PERMITTIVITY_SI():
    '''SI/AU (derived)'''
    return 4*PI()*PERMITTIVITY_FREE_SPACE_SI()




#  !
#  ! Dependent quantities
#  !


def AU_DISTANCE_SI():
    '''SI/AU (derived)'''
    return (AU_PERMITTIVITY_SI()*AU_ANGULAR_MOMENTUM_SI()*AU_ANGULAR_MOMENTUM_SI())/(AU_MASS_SI()*AU_CHARGE_SI()*AU_CHARGE_SI())



def AU_ENERGY_SI():
    '''SI/AU (derived)'''
    return (AU_CHARGE_SI()*AU_CHARGE_SI())/(AU_PERMITTIVITY_SI()*AU_DISTANCE_SI())



def AU_VELOCITY_SI():
    '''SI/AU (derived)'''
    return (AU_CHARGE_SI()*AU_CHARGE_SI())/(AU_PERMITTIVITY_SI()*AU_ANGULAR_MOMENTUM_SI())



def AU_TIME_SI():
    ''' SI/AU (derived)'''
    return (AU_DISTANCE_SI()/AU_VELOCITY_SI())



def AU_ELECTROSTATIC_POTENTIAL_SI():
    '''SI/AU (derived)'''
    return (AU_CHARGE_SI())/(AU_PERMITTIVITY_SI()*AU_DISTANCE_SI())



def AU_FORCE_SI():
    '''SI/AU (derived)'''
    return (AU_CHARGE_SI()*AU_CHARGE_SI())/(AU_PERMITTIVITY_SI()*AU_DISTANCE_SI()*AU_DISTANCE_SI())



def AU_MAGNETIC_DIPOLE_SI():
    '''SI/AU (derived)'''
    return (AU_CHARGE_SI()*AU_ANGULAR_MOMENTUM_SI())/AU_MASS_SI()



def AU_PRESSURE_SI():
    '''SI/AU (derived)'''
    return AU_FORCE_SI() / (AU_DISTANCE_SI()*AU_DISTANCE_SI())
  
  
  

#  !
#  ! Fundamental physical constants (AU)
#  !

def SPEED_OF_LIGHT_AU():
    ''' AU'''
    return SPEED_OF_LIGHT_SI() / AU_VELOCITY_SI()



def BOLTZMANN_CONSTANT_AU():
    '''AU'''
    return BOLTZMANN_CONSTANT_SI()/AU_ENERGY_SI()



def PLANCK_CONSTANT_AU():
    ''' AU'''
    return PLANCK_CONSTANT_SI()/(AU_TIME_SI()*AU_ENERGY_SI())



def PROTON_REST_MASS_AU():
    ''' AU'''
    return PROTON_REST_MASS_SI() / AU_MASS_SI()



def NEUTRON_REST_MASS_AU():
    ''' AU'''
    return NEUTRON_REST_MASS_SI() / AU_MASS_SI()


def RYDBERG_CONSTANT_AU():
    '''AU'''
    return RYDBERG_CONSTANT_SI() * AU_DISTANCE_SI()




#  !
#  ! Unit conversions (in terms of AU)
#  !


def AU_PER_PASCAL():
    '''AU/PASCAL (Conversion factor)'''
    return 1.0 / AU_PRESSURE_SI()



def AU_PER_ATMOSPHERE():
    '''AU/ATMOSPHERE  (Conversion factor)'''
    return AU_PER_PASCAL() * 101325.0



def AU_PER_ELECTRON_VOLT():
    '''AU/ELECTRON_VOLT  (Conversion factor)'''
    return AU_CHARGE_SI()/AU_ENERGY_SI()



def AU_PER_JOULE_PER_MOL():
    '''AU/JOULE_PER_MOL  (Conversion factor)'''
    return 1.0/(AU_ENERGY_SI()*AVOGADRO_CONSTANT())



def AU_PER_KJOULE_PER_MOL():
    '''AU/KJOULE_PER_MOL  (Conversion factor)'''
    return 1000.0 * AU_PER_JOULE_PER_MOL()



def AU_PER_KCAL_PER_MOL():
    '''AU/KCAL_PER_MOL  (Conversion factor)'''
    return 4.184 * AU_PER_KJOULE_PER_MOL()



def AU_PER_KCAL_PER_ANG_MOL():
    '''AU/KCAL_PER_ANG_MOL  (Conversion factor)'''
    return 4.184e13 / (AU_FORCE_SI() * AVOGADRO_CONSTANT())



def AU_PER_INVERSE_CM():
    '''AU/INVERSE_CM  (Conversion factor)'''
    return (100.*SPEED_OF_LIGHT_SI()*PLANCK_CONSTANT_SI())/AU_ENERGY_SI()



def AU_PER_ANGSTROM():
    '''AU/ANGSTROM  (Conversion factor)'''
    return 1.0e-10/AU_DISTANCE_SI()



def AU_PER_DEBYE():
    '''AU/DEBYE  (Conversion factor)'''
    return 3.335641e-30/(AU_CHARGE_SI()*AU_DISTANCE_SI())



def AU_PER_GRAMS_PER_MOL():
    '''AU/GRAMS_PER_MOL  (Conversion factor)'''
    return 1.0 / (1000.0 * AVOGADRO_CONSTANT()  * AU_MASS_SI())



def AU_PER_CM3_PER_MOL():
    '''AU/CM3_PER_MOL  (Conversion factor)'''
    return AVOGADRO_CONSTANT() * (100.0 * AU_DISTANCE_SI())*(100.0 * AU_DISTANCE_SI())*(100.0 * AU_DISTANCE_SI())



def AU_PER_ATOMIC_MASS_UNIT():
    '''AU/ATOMIC_MASS_UNIT  (Conversion factor)'''
    return 1.0/(AVOGADRO_CONSTANT() * 1000.0 * AU_MASS_SI()) #  ! au/amu


def AU_PER_RYDBERG():
    '''AU/RYDBERG  (Conversion factor)'''
    return RYDBERG_CONSTANT_SI() * AU_PER_INVERSE_CM() / 100.

def AU_PER_DYNE_PER_CM():
    '''AU/(DYNE/CM) (Conversion factor)'''
    # (100 cm/1 m) (m/au_distance) / ( 1e5 dyne/N) (N/au_force) ) = (cm/dyne)*(au_force/au_distance)
    return 0.001 * AU_DISTANCE_SI()/AU_FORCE_SI()

def AU_PER_PKA_PER_KELVIN():
    '''dG(T) (Hartree/mole/K) = k_b ln(10) pKa'''
    return BOLTZMANN_CONSTANT_AU()*LN10()

def AU_PER_PKA_AT_300K():
    '''dG(300K) (Hartree/mole) = k_b 300K ln(10) pKa'''
    return AU_PER_PKA_PER_KELVIN() * 300.

def AU_PER_PKA_AT_25C():
    '''dG(300K) (Hartree/mole) = k_b 300K ln(10) pKa'''
    return AU_PER_PKA_PER_KELVIN() * 298.15

def KT_AT_298K():
    '''AU of Boltzmann Constant * 298 K'''
    return BOLTZMANN_CONSTANT_AU() * 298.

def AU_PER_NM():
    '''AU/nanometer (conversion factor)'''
    return 10. * AU_PER_ANGSTROM()

def AU_PER_NM3():
    '''AU/nanometer^3 (conversion factor)'''
    return AU_PER_NM() * AU_PER_NM() * AU_PER_NM()

def AU_PER_ANG3():
    '''AU/Angstrom^3 (conversion factor)'''
    return AU_PER_ANGSTROM() * AU_PER_ANGSTROM() * AU_PER_ANGSTROM()

def AU_PER_GRAMS_PER_CM3():
    '''AU/(g/mol) (conversion factor)'''
    return AU_PER_GRAMS_PER_MOL() / AU_PER_CM3_PER_MOL()

def AU_PER_CAL_PER_MOL():
    '''AU/(cal/mol) (conversion factor)'''
    return AU_PER_KCAL_PER_MOL() / 1000.

def AU_PER_BAR():
    '''AU/bar (conversion factor)'''
    return  1.e+5 * AU_PER_PASCAL()
#    return 1.e+5 / (  AU_FORCE_SI() / ( AU_DISTANCE_SI() * AU_DISTANCE_SI() ) )

def AU_PER_METER():
    '''AU/meter (conversion factor)'''
    return 1. / AU_DISTANCE_SI()

def AU_PER_CM():
    '''AU/cm (conversion factor)'''
    return 1.e-3 / AU_DISTANCE_SI()

def AU_PER_M3():
    '''AU/meter^3 (conversion factor)'''
    return AU_PER_METER() * AU_PER_METER() * AU_PER_METER();

def AU_PER_CM3():
    '''AU/cm^3 (conversion factor)'''
    return AU_PER_CM() * AU_PER_CM() * AU_PER_CM()

def AU_PER_FM2():
    '''AU / (Farad * meter^2) (conversion factor)'''
    return AU_PER_M3() / ( FOUR_PI() * PERMITTIVITY_FREE_SPACE_SI() )




if __name__ == '__main__':

    print("%40s %23.14e" % ("PI()  = ",PI()  ))
    print("%40s %23.14e" % ("PI_D_180() = ",PI_D_180() ))
    print("%40s %23.14e" % ("PIO2() = ",PIO2() ))
    print("%40s %23.14e" % ("TWO_PI() = ",TWO_PI() ))
    print("%40s %23.14e" % ("FOUR_PI() = ",FOUR_PI() ))
    print("%40s %23.14e" % ("SQRT_PI() = ",SQRT_PI() ))
    print("%40s %23.14e" % ("SQRT2() = ",SQRT2() ))
    print("%40s %23.14e" % ("SQRT3() = ",SQRT3() ))
    print("%40s %23.14e" % ("EULER() = ",EULER() ))
    print("%40s %23.14e" % ("RAD2DEG() = ",RAD2DEG() ))
    print("%40s %23.14e" % ("DEG2RAD() = ",DEG2RAD() ))
    print("%40s %23.14e" % ("AVOGADRO_CONSTANT() = ",AVOGADRO_CONSTANT() ))
    
    print("\n\n")
    
    print("%40s %23.14e" % ("SPEED_OF_LIGHT_SI() = ",SPEED_OF_LIGHT_SI() ))
    print("%40s %23.14e" % ("BOLTZMANN_CONSTANT_SI() = ",BOLTZMANN_CONSTANT_SI() ))
    print("%40s %23.14e" % ("PLANCK_CONSTANT_SI() = ",PLANCK_CONSTANT_SI() ))
    print("%40s %23.14e" % ("HBAR_SI() = ",HBAR_SI() ))
    print("%40s %23.14e" % ("ELECTRON_CHARGE_SI() = ",ELECTRON_CHARGE_SI() ))
    print("%40s %23.14e" % ("ELECTRON_REST_MASS_SI() = ",ELECTRON_REST_MASS_SI() ))
    print("%40s %23.14e" % ("PROTON_REST_MASS_SI() = ",PROTON_REST_MASS_SI() ))
    print("%40s %23.14e" % ("NEUTRON_REST_MASS_SI() = ",NEUTRON_REST_MASS_SI() ))
    print("%40s %23.14e" % ("PERMITTIVITY_FREE_SPACE_SI() = ",PERMITTIVITY_FREE_SPACE_SI() ))
    print("%40s %23.14e" % ("FARADAY_CONSTANT_SI() = ",FARADAY_CONSTANT_SI() ))
    print("%40s %23.14e" % ("FINE_STRUCTURE_CONSTANT_SI() = ",FINE_STRUCTURE_CONSTANT_SI() ))
    print("%40s %23.14e" % ("RYDBERG_CONSTANT_SI() = ",RYDBERG_CONSTANT_SI() ))
    
    print("\n\n")


    print("%40s %23.14e" % ("AU_MASS_SI() = ",AU_MASS_SI() ))
    print("%40s %23.14e" % ("AU_CHARGE_SI() = ",AU_CHARGE_SI() ))
    print("%40s %23.14e" % ("AU_ANGULAR_MOMENTUM_SI() = ",AU_ANGULAR_MOMENTUM_SI() ))
    print("%40s %23.14e" % ("AU_PERMITTIVITY_SI() = ",AU_PERMITTIVITY_SI() ))
    print("%40s %23.14e" % ("AU_DISTANCE_SI() = ",AU_DISTANCE_SI() ))
    print("%40s %23.14e" % ("AU_ENERGY_SI() = ",AU_ENERGY_SI() ))
    print("%40s %23.14e" % ("AU_VELOCITY_SI() = ",AU_VELOCITY_SI() ))
    print("%40s %23.14e" % ("AU_TIME_SI() = ",AU_TIME_SI() ))
    print("%40s %23.14e" % ("AU_ELECTROSTATIC_POTENTIAL_SI() = ",AU_ELECTROSTATIC_POTENTIAL_SI() ))
    print("%40s %23.14e" % ("AU_FORCE_SI() = ",AU_FORCE_SI() ))
    print("%40s %23.14e" % ("AU_MAGNETIC_DIPOLE_SI() = ",AU_MAGNETIC_DIPOLE_SI() ))
    print("%40s %23.14e" % ("AU_PRESSURE_SI() = ",AU_PRESSURE_SI() ))
    
    print("\n\n")
    
    
    print("%40s %23.14e" % ("SPEED_OF_LIGHT_AU() = ",SPEED_OF_LIGHT_AU() ))
    print("%40s %23.14e" % ("BOLTZMANN_CONSTANT_AU() = ",BOLTZMANN_CONSTANT_AU() ))
    print("%40s %23.14e" % ("PLANCK_CONSTANT_AU() = ",PLANCK_CONSTANT_AU() ))
    print("%40s %23.14e" % ("PROTON_REST_MASS_AU() = ",PROTON_REST_MASS_AU() ))
    print("%40s %23.14e" % ("NEUTRON_REST_MASS_AU() = ",NEUTRON_REST_MASS_AU() ))
    print("%40s %23.14e" % ("RYDBERG_CONSTANT_AU() = ",RYDBERG_CONSTANT_AU() ))
    
    
    print("\n\n")
    
    
    print("%40s %23.14e" % ("AU_PER_PASCAL() = ",AU_PER_PASCAL() ))
    print("%40s %23.14e" % ("AU_PER_ATMOSPHERE() = ",AU_PER_ATMOSPHERE() ))
    print("%40s %23.14e" % ("AU_PER_BAR() = ",AU_PER_BAR() ))
    print("%40s %23.14e" % ("AU_PER_ELECTRON_VOLT() = ",AU_PER_ELECTRON_VOLT() ))
    print("%40s %23.14e" % ("AU_PER_JOULE_PER_MOL() = ",AU_PER_JOULE_PER_MOL() ))
    print("%40s %23.14e" % ("AU_PER_KJOULE_PER_MOL() = ",AU_PER_KJOULE_PER_MOL() ))
    print("%40s %23.14e" % ("AU_PER_CAL_PER_MOL() = ",AU_PER_CAL_PER_MOL() ))
    print("%40s %23.14e" % ("AU_PER_KCAL_PER_MOL() = ",AU_PER_KCAL_PER_MOL() ))
    print("%40s %23.14e" % ("AU_PER_KCAL_PER_ANG_MOL() = ",AU_PER_KCAL_PER_ANG_MOL() ))

    print("%40s %23.14e" % ("AU_PER_INVERSE_CM() = ",AU_PER_INVERSE_CM() ))
    print("%40s %23.14e" % ("AU_PER_ANGSTROM() = ",AU_PER_ANGSTROM() ))
    print("%40s %23.14e" % ("AU_PER_METER() = ",AU_PER_METER() ))
    print("%40s %23.14e" % ("AU_PER_CM() = ",AU_PER_CM() ))
    print("%40s %23.14e" % ("AU_PER_NM() = ",AU_PER_NM() ))
    print("%40s %23.14e" % ("AU_PER_M3() = ",AU_PER_M3() ))
    print("%40s %23.14e" % ("AU_PER_CM3() = ",AU_PER_CM3() ))
    print("%40s %23.14e" % ("AU_PER_ANG3() = ",AU_PER_ANG3() ))
    print("%40s %23.14e" % ("AU_PER_NM3() = ",AU_PER_NM3() ))
    print("%40s %23.14e" % ("AU_PER_FM2() = ",AU_PER_FM2() ))
    print("%40s %23.14e" % ("AU_PER_DEBYE() = ",AU_PER_DEBYE() ))
    print("%40s %23.14e" % ("AU_PER_GRAMS_PER_MOL() = ",AU_PER_GRAMS_PER_MOL() ))
    print("%40s %23.14e" % ("AU_PER_CM3_PER_MOL() = ",AU_PER_CM3_PER_MOL() ))
    print("%40s %23.14e" % ("AU_PER_GRAMS_PER_CM3() = ",AU_PER_GRAMS_PER_CM3() ))
    print("%40s %23.14e" % ("AU_PER_ATOMIC_MASS_UNIT() = ",AU_PER_ATOMIC_MASS_UNIT() ))
    print("%40s %23.14e" % ("AU_PER_RYDBERG() = ",AU_PER_RYDBERG() ))


