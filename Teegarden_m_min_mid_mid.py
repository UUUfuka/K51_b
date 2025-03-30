import posidonius
import numpy as np
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename
    #filename = posidonius.constants.BASE_DIR+"target/case7.json"

    #initial_time = 4.5e6*365.25 # time [days] where simulation starts
    initial_time = 8.e9*365.25 # time [days] where simulation starts 
    time_step = 0.08/1.e1 # days
    #time_step = 0.05 # days
    #time_limit   = 4*time_step # days
    time_limit   = 365.25 * 1.0e5 # days
    historic_snapshot_period = 365.25/1.e1 # days
    recovery_snapshot_period = 100.*historic_snapshot_period # days
    consider_tides = True
    consider_rotational_flattening = True
    #consider_general_relativity = False
    consider_general_relativity = "Kidder1995" # Assumes one central massive body
    #consider_general_relativity = "Anderson1975" # Assumes one central massive body
    #consider_general_relativity = "Newhall1983" # Considers all bodies
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_general_relativity)

    star_mass = 0.097 # 
    star_radius_factor = 0.120 
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_dissipation_factor = 2.006*3.845764e4  # trappist1
    star_dissipation_factor_scale = 0.01
    star_radius_of_gyration_2 = 2.00e-1 # 
    star_love_number = 0.307 # 
    fluid_love_number = star_love_number
    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    star_rotation_period = 96.2*24. # hours #The CARMENES search for exoplanets around M dwarfs
    #Mapping stellar activity indicators across the M dwarf domain
    #star_rotation_period = 19.0*24 # hours
    star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    star_spin = posidonius.Axes(0., 0., star_angular_frequency)
    #star_evolution_type = posidonius.GalletBolmont2017(star_mass) # mass = 0.30 .. 1.40
    #star_evolution_type = posidonius.BolmontMathis2016(star_mass) # mass = 0.40 .. 1.40
    #star_evolution_type = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
    #star_evolution_type = posidonius.Leconte2011(star_mass) # mass = 0.01 .. 0.08
    #star_evolution_type = posidonius.Baraffe1998(star_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #star_evolution_type = posidonius.LeconteChabrier2013() # Jupiter
    star_evolution_type = posidonius.NonEvolving()
    universe.add_particle(star_mass, star_radius, star_dissipation_factor, star_dissipation_factor_scale, star_radius_of_gyration_2, star_love_number, fluid_love_number, star_position, star_velocity, star_spin, star_evolution_type)

  
    # Radiuses in R_EARTH
    planet_radiuses = (1.083537445281188, 1.0531407937532873 ,0.9808538041740507 ) #Exoplex (~~~)
    # Masses in M_SUN
    planet_masses = (3.483172777185608e-06,3.152871910383525e-06,2.4622428252518955e-06)
    # Semi-major axis in AU
    planet_a = (0.025966039916238987,0.0455941832496055,0.07918787875458583)
    # Eccentricity
    planet_e = ( 0.03, 0.04, 0.07)
    # Inclination in degrees
    planet_i = (90., 90.,90.)
    # Mean anomaly in degrees
    planet_l = (56.8,115.9,200.4)
    # argument of pericentre in degrees
    planet_p = (338,301,345) 

    for r, m, a, e, i, l, p in zip(planet_radiuses, planet_masses, planet_a, planet_e, planet_i, planet_l, planet_p):
        planet_mass = m # Solar masses (3.0e-6 solar masses = 1 earth mass)
        planet_radius_factor = r # R Earth
        planet_radius = planet_radius_factor * posidonius.constants.R_EARTH
        # Terrestrial:
        k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
        planet_dissipation_factor = 2. * posidonius.constants.G * k2pdelta/(3. * np.power(planet_radius, 5))
        planet_dissipation_factor_scale = 1.0
        #planet_dissipation_factor_scale = 0.1
        planet_radius_of_gyration_2 = 0.3308
        planet_love_number = 0.299
        planet_fluid_love_number = 0.9532

        #////////// Specify initial position and velocity for a stable orbit
        #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
        #a = 0.01111;                             # semi-major axis (in AU)
        e = e;                              # eccentricity
        i = i * posidonius.constants.DEG2RAD;                      # inclination (degrees)
        p = p * posidonius.constants.DEG2RAD;                                # argument of pericentre (degrees)
        n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
        l = l * posidonius.constants.DEG2RAD;           # mean anomaly (degrees)
        p = (p + n);                 # Convert to longitude of perihelion !!
        q = a * (1.0 - e);                     # perihelion distance
        gm = posidonius.constants.G*(planet_mass+star_mass);
        x, y, z, vx, vy, vz = posidonius.calculate_cartesian_coordinates(gm, q, e, i, p, n, l);
        planet_position = posidonius.Axes(x, y, z)
        planet_velocity = posidonius.Axes(vx, vy, vz)

        #////// Initialization of planetary spin
        planet_obliquity = 1.0e-4 # rad
        # Pseudo-synchronization period
        planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(a, e, star_mass, planet_mass) # days
        planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1
        planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star_mass*planet_mass, planet_position, planet_velocity)
        planet_inclination = planet_keplerian_orbital_elements[3]
        planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)

        planet_evolution_type = posidonius.NonEvolving()
        universe.add_particle(planet_mass, planet_radius, planet_dissipation_factor, planet_dissipation_factor_scale, planet_radius_of_gyration_2, planet_love_number, planet_fluid_love_number, planet_position, planet_velocity, planet_spin, planet_evolution_type)

    ############################################################################

    whfast_alternative_coordinates="DemocraticHeliocentric"
    #whfast_alternative_coordinates="WHDS"
    #whfast_alternative_coordinates="Jacobi"
    universe.write(filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)
    #universe.write(filename, integrator="IAS15")
    #universe.write(filename, integrator="LeapFrog")
