import numpy as np





 
#_______________________________________________________________________________________________________________________
 
# Function Starts Here
def SummaryFunction(
                    dataframe           =False,
                    wall                =0.55,
                    capillar_width      =0.05,  # mm
                    air_width           =0.1, 
                    fiber_diameter      =2, 
                    Z                   =1,     # Unit length
                    ):

############################ I    MATH FUNCTIONS ############################################################################################

    ##################
    # MATH FUNCTIONS #
    ##################

    def hollow_cylinder_volume(D, d, h):
        '''D = outer diameter, d = inner diameter, h = height'''
        return np.pi * h * (D**2 - d**2)/4

    def cylinder_volume(R, h):
        '''R = radius, h = height'''
        return np.pi * h * R**2

    def circle_area(R):
        '''R = radius'''
        return np.pi * R**2

    def hollow_circle_area(R, r):
        '''R = outer radius, r = inner radius'''
        return np.pi * R**2 - np.pi * r**2

    def compound_moliere_radius_gr(mass_fractions, moliere_radii):
        '''Returns the compound moliere radius in g/cm^2'''
        return 1 / np.sum(np.array(mass_fractions) / np.array(moliere_radii))

    def moliere_radius_cm(moliere_radius, density):
        '''Returns the moliere radius in cm'''
        return moliere_radius / density

    def compound_radiation_length_gr(mass_fractions, radiation_lengths):
        '''Returns the compound radiation length in g/cm^2'''
        return 1 / np.sum(np.array(mass_fractions) / np.array(radiation_lengths))

    def radiation_length_cm(radiation_length, density):
        '''Returns the radiation length in cm'''
        return radiation_length / density

############################ II   ATTRIBUTES SECTION #####################################################################################

    ##############
    # ATTRIBUTES #
    ##############
   
    # Absorber sizes mm
    X = 121.00
    Y = 121.00

    
    # Densities g/cm^3
    dens_fib = 1.05
    dens_air = 0.001225
    dens_cap = 7.9
    
    
    # Absorber Parameters
    #                    lead 
    #                      Pb 
    densities =         [11.34]       # g/cm^3
    mass_fraction =     [1]       # sum = 1
    rad_lengths =       [6.37]       # g/cm^2
    moliere_radii =     [18.18]       # g/cm^2


    capillar_dict = {
                    "columns":         ['C',      "Si",    "Mn",      "Ni",    "S",    "P",     "Cr",      "Cu",   "Fe" ],
                    "mass_fraction":   [0.0012,  0.008,   0.002,    0.0010,  0.0002, 0.00035,  0.0018,    0.003,  0.68725],
                    "rad_lengths":     [42.7,    21.82,   14.64,     12.68,   19.50,  21.21,    14.94,    12.86,   13.84 ],
                    "moliere_radii":   [11.03,   11.51,   13.74,     13.41,   11.53,  11.86,    13.47,    14.05,   13.53 ]
                    }
    
    # Radiation Lengths Fiber and Air, g/cm^2
    X0_fib = 43.79
    X0_air = 36.62
    
    # Moliere Radii of Fiber (Polystyrene and Air, g/cm^2
    Rm_fib = 9.97
    Rm_air = 8.83

    air_radius = air_width/2

    # Armour inner diam
    capillar_d = fiber_diameter + (air_radius * 2)

    # Armour outer diam
    capillar_D = capillar_d + capillar_width

############################ III  NUMBER OF FIBERS ###########################################################################################

    ####################
    # NUMBER OF FIBRES #
    ####################


    # Rounding precision according to wall decimal points
    pitch_precision = len(str(wall))

    # HOLE for Fiber, Air, Capillar
    hole_diameter = round(fiber_diameter + capillar_width + air_radius*2, pitch_precision)

    # Pitch of the given configuration (distance between two neighbouring fiber centers)
    pitch = round(hole_diameter + wall, pitch_precision)

    # Number of Fibres per axis
    fibres_per_axis = np.ceil(((X - fiber_diameter) / pitch))
 
    # Number of Fibres per Module (w/o staggering!)
    fibers_per_module = fibres_per_axis * fibres_per_axis # Only if we don't use staggering!
 
############################ IV   VOLUMES SECTION ######################################################################################
 
    ###################
    # VOLUMES SECTION #
    ###################
 
    # Full Module volume cm^3
    full_vol = (X * Y * Z)/1000

    # Volume of capillar cm^3
    cap_one_vol = hollow_cylinder_volume(capillar_D, capillar_d, Z)/1000

    # Volume of all Capillar
    cap_all_vol = cap_one_vol * fibers_per_module
 
    # Volume of all Fibres w/o staggering cm^3
    fiber_all_vol = cylinder_volume(fiber_diameter/2, Z) * fibers_per_module/1000
 
    # Volume of all Air+Fibres w/o staggering cm^3
    air_fiber_all_vol = cylinder_volume(((fiber_diameter / 2) + air_radius), Z) * fibers_per_module/1000
 
    # Volume of all Air cm^3
    air_all_vol = air_fiber_all_vol - fiber_all_vol

    # Volume of Hole (Fiber, Air, Capillar)
    hole_all_vol = fiber_all_vol + air_all_vol + cap_all_vol
 
    # Volume of Full Absorber (not Module!) cm^3
    abs_vol = full_vol - hole_all_vol
 
    # Volume per gram Absorber cm^3
    abs_per_gram_vol = np.sum(np.array(mass_fraction) / np.array(densities))

############################ V    DENSITY SECTION ###########################################################################################
 
    ###########
    # DENSITY #
    ###########
 
    # Absorber density g/cm^3
    dens_abs = 1 / abs_per_gram_vol
 
    # ASK EVGENII! Module density (not per gram, but full module values taken)
    dens_module = (abs_vol * dens_abs + fiber_all_vol * dens_fib + air_all_vol * dens_air + cap_one_vol * dens_cap) / full_vol
 
############################ VI   MASS SECTION ###########################################################################################
 
    ################
    # MASS SECTION #
    ################
 
    # Mass of all Fibres
    mass_fib = fiber_all_vol * dens_fib
 
    # Mass of Air
    mass_air = air_all_vol * dens_air

    # Mass of Capillar
    mass_cap = cap_all_vol * dens_cap
 
    # Mass of Absorber (not Module!)
    mass_abs = abs_vol / abs_per_gram_vol
 
    # Mass of Full Module
    mass_full = mass_fib + mass_air + mass_abs + mass_cap
 
############################ VII  MASS FRACTIONS ###########################################################################################
 
    ##################
    # MASS FRACTIONS #
    ##################
 
    # Mass fraction of Fibers
    mfrac_fib = mass_fib / mass_full
 
    # Mass fraction of Air
    mfrac_air = mass_air / mass_full

    # Mass fraction of Capillar
    mfrac_cap = mass_cap / mass_full
 
    # Mass fraction of Absorber
    mfrac_abs = mass_abs / mass_full
 
############################ VIII SAMPLING FRACTIONS ###########################################################################################
 
    ######################
    # SAMPLING FRACTIONS #
    ######################
 
    ######## 
    # MASS #
    ########
    
    # Mass Sampling fraction (Should make another version for vector implementation)
    samp_frac_af = mfrac_abs / mfrac_fib
 
    # Mass Sampling fraction (Should make another version for vector implementation)
    samp_frac_fa = mfrac_fib / mfrac_abs
    
    ##########
    # VOLUME #
    ##########
    
    # Volume Sampling fraction (Should make another version for vector implementation)
    vol_samp_frac_af = abs_vol / fiber_all_vol
    
    # Volume Sampling fraction (Should make another version for vector implementation)
    vol_samp_frac_fa = fiber_all_vol / abs_vol
     
############################ IX   RADIATION LENGTHS ###########################################################################################
 
    #####################
    # RADIATION LENGTHS #
    #####################
 
    # Absorber Radiation Length g/cm^2
    X0_abs = compound_radiation_length_gr(mass_fraction, rad_lengths)

    # Capillar Radiation Length g/cm^2
    X0_cap = compound_radiation_length_gr(capillar_dict.get("mass_fraction"), capillar_dict.get("rad_lengths"))
 
    # Module Radiation Length g/cm^2
    X0_module = 1 / ((mfrac_fib / X0_fib) + (mfrac_air / X0_air) + (mfrac_abs / X0_abs) + (mfrac_cap / X0_cap))
 
    # Absorber Radiation Length cm
    X0_abs_cm = X0_abs / dens_abs

    # Capillar Radiation Length cm
    X0_cap_cm = X0_cap / dens_cap
 
    # Module Radiation Length cm (not per gram, but full module values taken)
    X0_module_cm = X0_module / dens_module
 
############################ X    MOLIERE RADII ###########################################################################################
 
    #################
    # MOLIERE RADII #
    #################
 
    # Absorber Moliere Radius g/cm^2
    Rm_abs = compound_moliere_radius_gr(mass_fraction, moliere_radii)

    # Absorber Moliere Radius g/cm^2
    Rm_cap = compound_moliere_radius_gr(capillar_dict.get("mass_fraction"), capillar_dict.get("moliere_radii"))
 
    # Module Moliere Radius g/cm^2
    Rm_module = 1 / ((mfrac_fib / Rm_fib) + (mfrac_air / Rm_air) + (mfrac_abs / Rm_abs) + (mfrac_cap / Rm_cap))
 
    # Absorber Moliere Radius cm
    Rm_abs_cm = Rm_abs / dens_abs

    # Absorber Moliere Radius cm
    Rm_cap_cm = Rm_cap / dens_cap
 
    # Module Moliere Radius cm (not per gram, but full module values taken)
    Rm_module_cm = Rm_module / dens_module
 
############################ XI   MODULE NEW CONFIGURATION ###########################################################################################
 
    ############################
    # MODULE NEW CONFIGURATION #
    ############################
 
    # Module 7X0 cm
    X0_7 = X0_module_cm * 7
 
    # Module 25X0 cm
    X0_25 = X0_module_cm * 25
 
    # absorber_size_z mm
    absorber_size_z = np.ceil(X0_25)*10


    ###################
    # VOLUMES SECTION #
    ###################

    # Z in mm
    Z = absorber_size_z
 
    # Full Module volume cm^3
    full_vol = (X * Y * Z)/1000
    
    # Volume of ONE Fiber cm^3
    fiber_one_vol = cylinder_volume(fiber_diameter/2, Z)/1000

    # Volume of ALL Fibres w/o staggering cm^3
    fiber_all_vol = (fiber_one_vol* fibers_per_module)
    
    # Volume of ONE Air+Fibers cm^3
    air_fiber_one_vol = cylinder_volume(air_radius + (fiber_diameter/2), Z)/1000

    # Volume of ALL Air+Fibres w/o staggering cm^3
    air_fiber_all_vol = (air_fiber_one_vol * fibers_per_module)
 
    # Volume of ALL Air cm^3
    air_all_vol = air_fiber_all_vol - fiber_all_vol

    # # Second option for ALL Air Volume cm^3
    # air_round_second_all = (((np.pi * Z * (fiber_diameter + air_radius*2)**2 - fiber_diameter**2)/4)*fibers_per_module)/1000
    # print(f'air_round_second_all: {air_round_second_all:.2f}')

    # Volume of ONE Capillar cm^3
    cap_one_vol = hollow_cylinder_volume(capillar_D, capillar_d, Z)/1000

    # Volume of ALL Capillars cm^3
    cap_all_vol = cap_one_vol * fibers_per_module


    # Volume of ALL Holes cm^3
    hole_all_vol = fiber_all_vol + air_all_vol + cap_all_vol


    # Volume of Full Absorber (not Module!) cm^3
    abs_vol = full_vol - hole_all_vol
 
    # Volume per gram Absorber cm^3
    abs_per_gram_vol = np.sum(np.array(mass_fraction) / np.array(densities))


    ################
    # MASS SECTION #
    ################
 
    # Mass of all Fibres
    mass_fib = fiber_all_vol * dens_fib
 
    # Mass of Air
    mass_air = air_all_vol * dens_air

    # Mass of Capillar
    mass_cap = cap_all_vol * dens_cap
 
    # Mass of Absorber (not Module!)
    mass_abs = abs_vol / abs_per_gram_vol
 
    # Mass of Full Module
    mass_full = mass_fib + mass_air + mass_abs + mass_cap



    ##########################
    # VALUES FOR CONFIG FILE #
    ##########################

    # cell_crystal_size_z = | 0 | 1 |         Cell Size Z (absolute value)
    cell_crystal_size_z_0 = np.ceil(X0_7*10)
    cell_crystal_size_z_1 = np.ceil(X0_25)*10 - np.ceil(X0_7*10)
 
    # cell_pos_z = | 0 | 1 |                  Cell position Z (coordinates)
    cell_pos_z_0 = -(cell_crystal_size_z_1) / 2
    cell_pos_z_1 = cell_crystal_size_z_0 / 2
 
    # cell_separator_position = ...           Cell Separator position (between cells | 0 | 1 |)
    cell_separator_position = cell_pos_z_0 + cell_crystal_size_z_0/2
  
############################ SAVING TO A DATAFRAME ##################################################################################################################
    if dataframe == True:

        import pandas as pd
        
        dict_DF = {
            "X": [X, "mm"],
            "Y": [Y, "mm"],
            "Z": [Z, "mm"],
            "Pitch": [pitch, "mm"],
            "Hole diameter": [hole_diameter, "mm"],
            "Fibers per Axis": [fibres_per_axis, "un."],
            "Fibers per Module": [fibers_per_module, "un."],
            "Full Volume": [full_vol, "mm^3"],
            "Fibres Volume": [fiber_all_vol, "mm^3"],
            "Air Volume": [air_all_vol, "mm^3"],
            "Capillar Volume": [cap_all_vol, "mm^3"],
            "Hole Volume": [hole_all_vol, "mm^3"],
            "Absorber Volume": [abs_vol, "mm^3"],
            "Absorber Volume per Gram": [abs_per_gram_vol, "mm^3"],
            "Absorber Density": [dens_abs, "g/cm^3"],
            "Module Density": [dens_module, "g/cm^3"],
            "All Fibers Mass": [mass_fib, "gr"],
            "All Air Mass": [mass_air, "gr"],
            "All Capillar Mass": [mass_cap, "gr"],
            "Absorber Mass": [mass_abs, "gr"],
            "Module Mass": [mass_full, "gr"],
            "Fibers Mass Fractions": [mfrac_fib, "/1"],
            "Air Mass Fractions": [mfrac_air, "/1"],
            "Capillar Mass Fractions": [mfrac_cap, "/1"],
            "Absorber Mass Fractions": [mfrac_abs, "/1"],
            "Mass Sampling Fraction (absorber/fibres)": [samp_frac_af, "abs/fib [mas]"],
            "Mass Sampling Fraction (fibers/absorber)": [samp_frac_fa, "fib/abs [mas]"],
            "Volume Sampling Fraction (absorber/fibers)": [vol_samp_frac_af, "abs/fib [vol]"],
            "Volume Sampling Fraction (fibers/absorber)": [vol_samp_frac_fa, "fib/abs [vol]"],
            "Capillar Rad Length [gr/cm^2]": [X0_cap, "g/cm^2"],
            "Absorber Rad Length [gr/cm^2]": [X0_abs, "g/cm^2"],
            "Module Rad Length [gr/cm^2]": [X0_module, "g/cm^2"],
            "Capillar Rad Length [cm]": [X0_cap_cm, "cm"],
            "Absorber Rad Length [cm]": [X0_abs_cm, "cm"],
            "Module Rad Length [cm]": [X0_module_cm, "cm"],
            "Capillar Moliere Radius [gr/cm^2]": [Rm_cap, "g/cm^2"],
            "Absorber Moliere Radius [gr/cm^2]": [Rm_abs, "g/cm^2"],
            "Module Moliere Radius [gr/cm^2]": [Rm_module, "g/cm^2"],
            "Capillar Moliere Radius [cm]": [Rm_cap_cm, "cm"],
            "Absorber Moliere Radius [cm]": [Rm_abs_cm, "cm"],
            "Module Moliere Radius [cm]": [Rm_module_cm, "cm"],
            "7X0":[X0_7, "cm"],
            "25X0":[X0_25, "cm"]
            }

            
        return pd.DataFrame(dict_DF, index=["Results", "Units"])
