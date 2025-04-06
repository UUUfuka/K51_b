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
    initial_time = 3.0e8*365.25 # time [days] where simulation starts
    time_step = 0.08*1.5e1 # days
    #time_step = 0.05 # days
    #time_limit   = 4*time_step # days
    time_limit   = 365.25 * 1.0e8 # days
    historic_snapshot_period = 365.25*15. # days 
    recovery_snapshot_period = 100.*historic_snapshot_period # days
    #consider_tides = True
    #consider_rotational_flattening = True
    #consider_general_relativity = False
    #consider_general_relativity = "Kidder1995" # Assumes one central massive body
    #consider_general_relativity = "Anderson1975" # Assumes one central massive body
    #consider_general_relativity = "Newhall1983" # Considers all bodies
    #universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_general_relativity)
    consider_effects = posidonius.ConsiderEffects({
         "tides": True,
         "rotational_flattening": True,
         "general_relativity": True,
         "disk": False,
         "wind": False,
         "evolution": False,
     })
    #universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_general_relativity)
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

    star_mass = 0.985 # The Featureless Transmission Spectra of Two Super-puff Planets #
    star_radius_factor = 0.881 #
    star_radius = star_radius_factor * posidonius.constants.R_SUN #
    star_radius_of_gyration_2 = 2.00e-1 # 
    star_position = posidonius.Axes(0., 0., 0.)#
    star_velocity = posidonius.Axes(0., 0., 0.)#
    star_rotation_period = 8.222*24. #The Featureless Transmission Spectra of Two Super-puff Planets
    star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    star_spin = posidonius.Axes(0., 0., star_angular_frequency)
    star=posidonius.Particle(star_mass,star_radius,star_radius_of_gyration_2,star_position,star_velocity,star_spin)

    star_tides_parameters = {"dissipation_factor_scale": 0.01,
                             "dissipation_factor": 2.006*3.845764e4,
                             "love_number": 0.307,}
    star_tides = posidonius.effects.tides.CentralBody(star_tides_parameters)

    star_fluid_love_number = 0.307 
    star_rotational_flattening_parameters = {"love_number": star_fluid_love_number }
    star_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(star_rotational_flattening_parameters)

    star_general_relativity = posidonius.effects.general_relativity.CentralBody("Kidder1995")

    #star_evolution_type = posidonius.GalletBolmont2017(star_mass) # mass = 0.30 .. 1.40
    #star_evolution_type = posidonius.BolmontMathis2016(star_mass) # mass = 0.40 .. 1.40
    #star_evolution_type = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
    #star_evolution_type = posidonius.Leconte2011(star_mass) # mass = 0.01 .. 0.08
    #star_evolution_type = posidonius.Baraffe1998(star_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #star_evolution_type = posidonius.LeconteChabrier2013() # Jupiter


    star_wind = posidonius.effects.wind.Disabled()
    star_disk = posidonius.effects.disk.Disabled()
    star_evolution = posidonius.NonEvolving()

    star.set_tides(star_tides)
    star.set_rotational_flattening(star_rotational_flattening)
    star.set_general_relativity(star_general_relativity)
    star.set_wind(star_wind)
    star.set_disk(star_disk)
    star.set_evolution(star_evolution)
    universe.add_particle(star)
    # universe.Particle(star_mass, star_radius, star_dissipation_factor, star_dissipation_factor_scale, star_radius_of_gyration_2, star_love_number, fluid_love_number, star_position, star_velocity, star_spin, star_evolution_type)

  #FIRST outside 2:1
  
    # Radiuses in R_EARTH
    planet_radiuses = (6.83, 6.4, 9.32, 6. ) 
    # Masses in M_SUN
    #planet_masses = (6.7, 6.09, 6.6, 5.4 ) 
    planet_masses = (2.0123378637300316e-05, 1.829125013450133e-05, 1.9823029702415233e-05, 1.6218842483794285e-05) 
    # Semi-major axis in AU
    planet_a = (0.24691452906203384, 0.3773679903833344, 0.500178611205457, 0.801917005809947)
    # Eccentricity
    planet_e = ( 0.0162, 0.0093, 0.0061,0.020)
    # Inclination in degrees
    planet_i = (0.,0.,0.,0.)
    # Mean anomaly in degrees
    planet_l = (0.0, 145.24391693342372, 213.68515767464623, 219.15542825142626)  #nominal from mean longitude #?(model dolgota voshod yzel = 0) normirovan na 2*pi 
    # argument of pericentre in degrees
    planet_p = (184.3987053549955, 354.28940686250036, 20.136303428248134, 211.72360295768402)

    for r, m, a, e, i, l, p in zip(planet_radiuses, planet_masses, planet_a, planet_e, planet_i, planet_l, planet_p):
        planet_mass = m # Solar masses (3.0e-6 solar masses = 1 earth mass)
        planet_radius_factor = r # R Earth
        planet_radius = planet_radius_factor * posidonius.constants.R_EARTH
        
        # Terrestrial:
        planet_radius_of_gyration_2 = 0.3308
        
        k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
        planet_tides_parameters = {
            "dissipation_factor_scale": 1.0,
            "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
            "love_number": 0.299,
        }
        planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_parameters)
        # planet_dissipation_factor = 2. * posidonius.constants.G * k2pdelta/(3. * np.power(planet_radius, 5))
        #posidonius.constants.K2 = posidonius.constants.G
        # planet_dissipation_factor_scale = 1.0
        #planet_dissipation_factor_scale = 0.1
        # planet_love_number = 0.299

        
        planet_rotational_flattening_parameters = {"love_number": 0.9532}
        planet_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(planet_rotational_flattening_parameters)
        # planet_fluid_love_number = 0.9532

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
        # gm = posidonius.constants.G*(planet_mass+star_mass);
        planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        # x, y, z, vx, vy, vz = posidonius._calculate_cartesian_coordinates(gm, q, e, i, p, n, l);
        # planet_position = posidonius.Axes(x, y, z)
        # planet_velocity = posidonius.Axes(vx, vy, vz)

        #////// Initialization of planetary spin
        planet_obliquity = 1.0e-4 # rad
        # Pseudo-synchronization period
        planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(a, e, star_mass, planet_mass) # days
        planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1
        planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity,masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        planet_inclination = planet_keplerian_orbital_elements[3]
        planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)

        planet_general_relativity = posidonius.effects.general_relativity.OrbitingBody()
        planet_wind = posidonius.effects.wind.Disabled()
        planet_disk = posidonius.effects.disk.Disabled()
        planet_evolution_type = posidonius.NonEvolving()
        planet_evolution = posidonius.NonEvolving()
        planet = posidonius.Particle(planet_mass, planet_radius, planet_radius_of_gyration_2, planet_position, planet_velocity, planet_spin)
        planet.set_tides(planet_tides)
        planet.set_rotational_flattening(planet_rotational_flattening)
        planet.set_general_relativity(planet_general_relativity)
        planet.set_wind(planet_wind)
        planet.set_disk(planet_disk)
        planet.set_evolution(planet_evolution)
        universe.add_particle(planet)


        # universe.add_particle(planet_mass, planet_radius, planet_dissipation_factor, planet_dissipation_factor_scale, planet_radius_of_gyration_2, planet_love_number, planet_fluid_love_number, planet_position, planet_velocity, planet_spin, planet_evolution_type)

    ############################################################################

    whfast_alternative_coordinates="DemocraticHeliocentric"
    #whfast_alternative_coordinates="WHDS"
    #whfast_alternative_coordinates="Jacobi"
    universe.write(filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)
    #universe.write(filename, integrator="IAS15")
    #universe.write(filename, integrator="LeapFrog")
