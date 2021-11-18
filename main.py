import numpy as np





 
#_______________________________________________________________________________________________________________________
 
# Function Starts Here
def SummaryFunction(
                    config              =False, 
                    print_results       =False, 
                    minimal_print       =False, 
                    return_values       =False, 
                    single_section      =False, 
                    automation          =False,
                    only_vis            =False,
                    DEF                 =False,
                    dataframe           =False,
                    wall                =0.55,
                    capillar_width      =0.05,  # mm
                    air_radius          =0.1, 
                    fiber_diameter      =2, 
                    seed                =42, 
                    absorber_material   =13,
                    Z                   =1,     # Unit length
                    energy              =1
                    ):

############################ I    MATH FUNCTIONS ############################################################################################

    ##################
    # MATH FUNCTIONS #
    ##################

    def hollow_cylinder_volume(D, d, h):
        return np.pi * h * (D**2 - d**2)/4

    def cylinder_volume(R, h):
        return np.pi * h * R**2

    def circle_area(R):
        return np.pi * R**2

    def hollow_circle_area(R, r):
        return np.pi * R**2 - np.pi * r**2

    def compound_moliere_radius_gr(mass_fractions, moliere_radii):
        return 1 / np.sum(np.array(mass_fractions) / np.array(moliere_radii))

    def moliere_radius_cm(moliere_radius, density):
        return moliere_radius / density

    def compound_radiation_length_gr(mass_fractions, radiation_lengths):
        return 1 / np.sum(np.array(mass_fractions) / np.array(radiation_lengths))

    def radiation_length_cm(radiation_length, density):
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
    #                    lead   antimony  tin
    #                      Pb      Sb     Sn
    densities =         [11.34,  6.697, 5.75]       # g/cm^3
    mass_fraction =     [0.84,   0.12,  0.04]       # sum = 1
    rad_lengths =       [6.37,   8.73,  8.82]       # g/cm^2
    moliere_radii =     [18.18, 15.81, 15.77]       # g/cm^2


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

############################ XII  PRINTING SECTION ######################################################################################################
 
    ####################
    # PRINTING SECTION #
    ####################
 
    if print_results == True:

        print("\nCONFIGURATION\n")
        print("X size of a module:_______| {} mm".format(X))
        print("Y size of a module:_______| {} mm".format(Y))
        print("Z size of a module:_______| {} mm".format(Z))
        print()
        print("Pitch:____________________| {} mm".format(pitch))
        print("Hole diameter:____________| {} mm".format(hole_diameter))
        # print(pitch_precision)
    
        print("\nNUMBER OF FIBERS\n")
        print('Fibers per Axis:__________|   {}'.format(fibres_per_axis))
        print('Fibers per Module_________| {}'.format(fibers_per_module))
    
        print("\nVOLUMES\n")
        print('Full Volume:______________| {} cm^3'.format(full_vol))
        print('Fibres Volume:____________| {:.2f} cm^3'.format(fiber_all_vol))
        print('Air + Fibres Volume:______| {:.2f} cm^3'.format(air_fiber_all_vol))
        print('Air Volume:_______________|  {:.2f} cm^3'.format(air_all_vol))
        print('Capillar Volume:__________|  {:.2f} cm^3'.format(cap_all_vol))
        print('Hole Volume:______________|  {:.2f} cm^3'.format(hole_all_vol))
        print('Absorber Volume:__________| {:.2f} cm^3'.format(abs_vol))
        print('Per gram Absorber Volume:_|    {:.2f} cm^3'.format(abs_per_gram_vol))
    
        print("\nDENSITIES\n")
        print('Absorber Density:_________| {:.2f} g/cm^3'.format(dens_abs))
        print('Module Density:___________|  {:.2f} g/cm^3'.format(dens_module))
    
        print("\nMASSES\n")
        print('All Fibres Mass:__________|  {:.2f} gr'.format(mass_fib))
        print('All Air Mass:_____________|     {:.2f} gr'.format(mass_air))
        print('All Capillar Mass:________|  {:.2f} gr'.format(mass_cap))
        print('Absorber Mass:____________| {:.2f} gr'.format(mass_abs))
        print('Module Mass:______________| {:.2f} gr'.format(mass_full))
    
        print("\nMASS FRACTIONS\n")
        print('Fibres Mass Fraction:_____| {:.2f}'.format(mfrac_fib))
        print('Air Mass Fraction:________| {:.10f}'.format(mfrac_air))
        print('Capillar Mass Fraction:___| {:.10f}'.format(mfrac_cap))
        print('Absorber Mass Fraction:___| {:.2f}'.format(mfrac_abs))
    
        print("\nMASS SAMPLING FRACTIONS\n")
        print('Mass Sampling Fraction (absorber/fibres): {:.2f}'.format(samp_frac_af))
        print('Mass Sampling Fraction (fibres/absorber): {:.2f}'.format(samp_frac_fa))
        
        print("\nVOLUME SAMPLING FRACTIONS\n")
        print('Volume Sampling Fraction (absorber/fibres): {:.2f}'.format(vol_samp_frac_af))
        print('Volume Sampling Fraction (fibres/absorber): {:.2f}'.format(vol_samp_frac_fa))
    
        print("\nRADIATION LENGTHS\n")
        print('Capillar Rad Length:______| {:.2f} g/cm^2'.format(X0_cap))
        print('Absorber Rad Length:______| {:.2f} g/cm^2'.format(X0_abs))
        print('Module Rad Length:________| {:.2f} g/cm^2'.format(X0_module))
        print('Capillar Rad Length:______| {:.2f} cm'.format(X0_cap_cm))
        print('Absorber Rad Length:______| {:.2f} cm'.format(X0_abs_cm))
        print('Module Rad Length:________| {:.2f} cm'.format(X0_module_cm))
    
        print("\nMOLIERE RADII\n")
        print('Capillar Moliere Radius:__| {:.2f} g/cm^2'.format(Rm_cap))
        print('Absorber Moliere Radius:__| {:.2f} g/cm^2'.format(Rm_abs))
        print('Module Moliere Radius:____| {:.2f} g/cm^2'.format(Rm_module))
        print('Capillar Moliere Radius:__| {:.2f} cm'.format(Rm_cap_cm))
        print('Absorber Moliere Radius:__| {:.2f} cm'.format(Rm_abs_cm))
        print('Module Moliere Radius:____| {:.2f} cm'.format(Rm_module_cm))
    
    
        print("\nPROPOSED SECTION LENGTHS\n")
        print('7X0:   {:.2f} cm'.format(X0_7))
        print('25X0:  {:.2f} cm'.format(X0_25))
    
############################ XIII CAN COMMENT BELOW: IT IS ALL FOR .CFG FILES ############################################################################################
    
        # ###########################################################
        # ####### CAN COMMENT BELOW: IT IS ALL FOR .CFG FILES #######
        # ###########################################################
    if config == True:


        ##################
        # DOUBLE SECTION #
        ##################


        if single_section == False:

    
            print("\nVALUES GO TO CONFIG\n")
            print('absorber_size_x = {}'.format(np.float(X)))
            print('absorber_size_y = {}'.format(np.float(Y)))
            print('absorber_size_z = {}'.format(absorber_size_z))
            print()
            print('cell_crystal_size_x = | {} | {} |'.format(fiber_diameter, fiber_diameter))
            print('cell_crystal_size_y = | {} | {} |'.format(fiber_diameter, fiber_diameter))
            print('cell_crystal_size_z = | {} | {} |'.format(cell_crystal_size_z_0,
                                                            cell_crystal_size_z_1))
            print()
            print('cell_crystal_pitch_x = | {} | {} |'.format(pitch, pitch))
            print('cell_crystal_pitch_y = | {} | {} |'.format(pitch, pitch))
            print()
            print('cell_pos_z = | {} | {} |'.format(cell_pos_z_0,
                                                    cell_pos_z_1))
            print('cell_separator_position = {}'.format(cell_separator_position))
            print()
            # print('cell_staggering_size = | {} | {} |'.format(pitch/2, pitch/2))
    
    
            #############################
            # WRITING TO A FILE SECTION #
            #############################
            
        
        
            lines = open('cfg_files/double_section_canvas.cfg').read().splitlines()
            # Seed
            lines[23] = 'seed = {}'.format(seed)
        
            # ABSORBER
            lines[305] = 'absorber_size_x   = {}'.format(X)
            lines[306] = 'absorber_size_y   = {}'.format(Y)
            lines[307] = 'absorber_size_z   = {}'.format(absorber_size_z)
        
            # ABSORBER MATERIAL
            lines[312] = 'absorber_material = {}'.format(absorber_material)
        
            # CELLS
            lines[342] = 'cell_pos_z  = |{}|{}|'.format(cell_pos_z_0, cell_pos_z_1)
            lines[343] = 'cell_x_elements  = |{}|{}|'.format(fibres_per_axis, fibres_per_axis)
            lines[344] = 'cell_y_elements  = |{}|{}|'.format(fibres_per_axis, fibres_per_axis)
            lines[347] = 'cell_crystal_size_x  = |{}|{}|'.format(fiber_diameter, fiber_diameter)
            lines[348] = 'cell_crystal_size_y  = |{}|{}|'.format(fiber_diameter, fiber_diameter)
            lines[349] = 'cell_crystal_size_z  = |{}|{}|'.format(cell_crystal_size_z_0, cell_crystal_size_z_1)
            lines[350] = 'cell_crystal_pitch_x  = |{}|{}|'.format(pitch, pitch)
            lines[351] = 'cell_crystal_pitch_y  = |{}|{}|'.format(pitch, pitch)
            lines[356] = 'cell_separator_position = {}'.format(cell_separator_position)
        
            # STAGGERING
            # lines[378] = 'cell_staggering_size   = |{}|{}|'.format(pitch / 2, pitch / 2)
        
            # OUTPUT FILE NAME
            output_name = 'Pb-Sb-Sn_diam{}_pitch{}_DOUBLE_SECTION_seed42.cfg'.format(fiber_diameter, pitch)
            open("my_configs/" + output_name, 'w').write('\n'.join(lines))
        




        elif single_section == True:


        ##################
        # SINGLE SECTION #
        ##################


    
            print("\nVALUES GO TO CONFIG\n")
            print('absorber_size_x = {}'.format(np.float(X)))
            print('absorber_size_y = {}'.format(np.float(Y)))
            print('absorber_size_z = {}'.format(absorber_size_z))
            print()
            print('cell_crystal_size_x = | {} |'.format(fiber_diameter))
            print('cell_crystal_size_y = | {} |'.format(fiber_diameter))
            print('cell_crystal_size_z = | {} |'.format(absorber_size_z))
            print()
            print('cell_crystal_pitch_x = | {} |'.format(pitch))
            print('cell_crystal_pitch_y = | {} |'.format(pitch))
            print()
            print('cell_pos_z = | {} |'.format(0))
            # print('cell_separator_position = {}'.format(cell_separator_position))
            print()
            # print('cell_staggering_size = | {} | {} |'.format(pitch/2, pitch/2))
    
    
            #############################
            # WRITING TO A FILE SECTION #
            #############################
            
        
        
            lines = open('my_configs/single_section_canvas.cfg').read().splitlines()
            # Seed
            lines[23] = 'seed = {}'.format(seed)
        
            # ABSORBER
            lines[305] = 'absorber_size_x   = {}'.format(X)
            lines[306] = 'absorber_size_y   = {}'.format(Y)
            lines[307] = 'absorber_size_z   = {}'.format(absorber_size_z)
        
            # ABSORBER MATERIAL
            lines[312] = 'absorber_material = {}'.format(absorber_material)
        
            # CELLS
            lines[342] = 'cell_pos_z  = |{}|'.format(0)
            lines[343] = 'cell_x_elements  = |{}|'.format(fibres_per_axis)
            lines[344] = 'cell_y_elements  = |{}|'.format(fibres_per_axis)
            lines[347] = 'cell_crystal_size_x  = |{}|'.format(fiber_diameter)
            lines[348] = 'cell_crystal_size_y  = |{}|'.format(fiber_diameter)
            lines[349] = 'cell_crystal_size_z  = |{}|'.format(absorber_size_z)
            lines[350] = 'cell_crystal_pitch_x  = |{}|'.format(pitch)
            lines[351] = 'cell_crystal_pitch_y  = |{}|'.format(pitch)
            # lines[356] = 'cell_separator_position = {}'.format(cell_separator_position)
        
            # STAGGERING
            # lines[378] = 'cell_staggering_size   = |{}|{}|'.format(pitch / 2, pitch / 2)
        
            # OUTPUT FILE NAME
            output_name = 'Pb-Sb-Sn_diam{}_pitch{}_SINGLE_SECTION_seed42.cfg'.format(fiber_diameter, pitch)
            open("my_configs/" + output_name, 'w').write('\n'.join(lines))


    
    # print(summary_dict)
    if minimal_print==True:

        print("\n")
        print('Fiber Diameter:_______________| {:.2f} mm'.format(fiber_diameter))
        print('Pitch:________________________| {:.2f} mm'.format(pitch))
        print('Absorber Rad Length:______| {:.2f} cm'.format(X0_abs_cm))
        print('Module Rad Length:________| {:.2f} cm'.format(X0_module_cm))
        print('Absorber Moliere Radius:__| {:.2f} cm'.format(Rm_abs_cm))
        print('Module Moliere Radius:____| {:.2f} cm'.format(Rm_module_cm))

    if return_values == True:
        return [fiber_diameter, pitch, wall, hole_diameter, air_radius, X0_module_cm, 
                Rm_module_cm, fibres_per_axis, fibers_per_module, samp_frac_af, 
                fiber_all_vol, full_vol, mass_fib, mass_air, mass_full, Z]

############################ XIV  AUTOMATION SECTION ##########################################################################################################

    ##############
    # AUTOMATION #
    ##############

    if automation == True:
        import ROOT
        import root_numpy as rn
        
        # Importing neccessary module - Python talking to Shell
        import subprocess

        # Variable subscription
        diam = fiber_diameter    # mm
        pitch = pitch           # mm

        # Control flow in for folder designation
        if single_section == False:
            section = "double_section"

        elif single_section == True:
            section = "single_section"



        mkdir_folder = f"mkdir root_outputs/{section}/{int(diam*10)}_{int(pitch*100)}"

        out_final_folder = f"{int(diam*10)}_{int(pitch*100)}"

        build_program = "build/FibresCalo"

        cfg_name = f"my_configs/Pb-Sb-Sn_diam{diam}_pitch{pitch}_{section.upper()}_seed42.cfg"

        out_path = f"root_outputs/{section}/{out_final_folder}"

        out_name = f"Pb_Sb_Polysterene_round_3x3deg_{diam*10}mm_{pitch*100}pitch_{energy}_GeV"

        gps_name = f"{energy}_GeV_my_gps.mac"




        subprocess.Popen(f"{mkdir_folder}", shell=True)

        

        if only_vis == True:
            subprocess.Popen(f"{build_program} {cfg_name}",  shell=True)

        elif only_vis == False:
            # subprocess.Popen(f"{build_program} {cfg_name}",  shell=True)
            subprocess.Popen(f"{build_program} {cfg_name} {out_path}/{out_name} my_gps/{gps_name}", shell=True)
        

        if DEF == True:
            myFile = f"{out_path}/{out_name}.root"

            results = np.array([])
    
            # File Object. Reading root file as variable. File Object
            rootFile = ROOT.TFile.Open(myFile)
            
            # Getting access to the "shower" branch
            rootTree = rootFile.Get("shower")
            
            # Turning rootTree into Numpy Array using root_numpy function
            numpyArray = rn.tree2array(rootTree, branches=["totalEnDep", "event"], selection=f'isInCrystal==1')
        


            
            for i in range(1000):
                
                results = np.append(results, numpyArray[numpyArray['event'] == i]['totalEnDep'].sum()/1000)
                
        #       


            return results

############################ SAVING TO A DATAFRAME ##################################################################################################################
    if dataframe == True:
        import pandas as pd
        


        indeces = np.arange(1,3)
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


    # print("\nPROPOSED SECTION LENGTHS\n")
    # print('7X0:   {:.2f} cm'.format(X0_7))
    # print('25X0:  {:.2f} cm'.format(X0_25))