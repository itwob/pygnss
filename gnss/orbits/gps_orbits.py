from numpy import ones, asarray, fabs, sqrt, sin, cos, pi, arccos, arctan2
import numpy
npmax = numpy.max

from gnss.time import gpsseconds


def compute_gps_satellite_position_from_yuma(alm, t, mu = 3.986005e14, Omega_e_dot = 7.2921151467e-5, compute_velocity=False, rollover=0):
    """Computes and returns the GPS satellite ECEF positions given a Yuma
    almanac and a set of times (seconds into GPS week).
    
    parameters:
    alm - the Yuma almanac (an ordered tuple containing Yuma almanac information)
        The Yuma almanac should contain: health, e, t_oa, i_0, Omega_dot, sqrt_a, Omega_week, omega, M_0, A_f0, A_f1, week
    t - the desired times for which to compute satellite locations
    
    optional parameters:
    mu - Earth gravitational constant (m^3/s^2)
    Omega_dot - WGS Earth rotation rate (rad/s)
    compute_velocity - flag indicating whether or not to compute satellite velocity as well (default False)
    """
    health, e, t_oa, i_0, Omega_dot, sqrt_a, Omega_week, omega, M_0, A_f0, A_f1, week = alm
    a = sqrt_a**2.
    delta_t = t - gpsseconds(week, t_oa, rollover)
    n_0 = sqrt(mu / a**3)
    n = n_0  # mean motion
    M = M_0 + n * delta_t  # mean anomaly
    E = 0.  # eccentric anomaly
    err = 1.
    tol = 1e-8
    while err > tol:  # iteratively solve for eccentric anomaly
        _temp = M + e * sin(E)
        err = npmax(fabs(_temp - E))
        E = _temp
    nu = arctan2(sqrt(1. - e**2) * sin(E), cos(E) - e)  # true anomaly
    E = arccos((e + cos(nu)) / (1. + e * cos(nu)))  # eccentric anomaly
    Phi = nu + omega  # argument of latitude
    u = Phi  # argument of latitude
    r = a * (1. - e * cos(E))  # orbital radius
    i = i_0  # inclination
    Omega = Omega_week + (Omega_dot - Omega_e_dot) * delta_t - Omega_e_dot * t_oa  # longitude of ascensing node (corrected)
    x_orb, y_orb = r * cos(u), r * sin(u)  # compute x,y in orbital plane
    x_ecef = x_orb * cos(Omega) - y_orb * sin(Omega) * cos(i)  # transform from orbital system to ECEF system
    y_ecef = x_orb * sin(Omega) + y_orb * cos(Omega) * cos(i)
    z_ecef = y_orb * sin(i) * ones(asarray(Omega).size)
    if compute_velocity:  # compute x,y velocity in orbital plane
        v_x_orb = n * a * sin(E) / (1. - e * cos(E))
        v_y_orb = -n * a * sqrt(1. - e**2) * cos(E) / (1. - e * cos(E))
        v_x_ecef = v_x_orb * cos(Omega) - v_y_orb * sin(Omega) * cos(i)  # transform from orbital system to ECEF system
        v_y_ecef = v_x_orb * sin(Omega) + v_y_orb * cos(Omega) * cos(i)
        v_z_ecef = v_y_orb * sin(i) * ones(asarray(Omega).size)
        return asarray([[x_ecef, y_ecef, z_ecef], [v_x_ecef, v_y_ecef, v_z_ecef]])
    else:
        return asarray([x_ecef, y_ecef, z_ecef])



def compute_gps_satellite_position_from_ephemeris(eph, t, mu=3.986005e14, Omega_e_dot=7.2921151467e-5, compute_velocity=False, rollover=0):
    """Computes and returns the GPS satellite ECEF positions given ephemeris
    parameters and a set of times.
    
    parameters:
    eph - the ephemeris set (an ordered tuple containing ephemeris information)
        `eph` should contain: e, t_oa, i_0, Omega_dot, sqrt_a, Omega_week, omega, M_0, A_f0, A_f1, week  TODODODO
    t - the desired times for which to compute satellite locations
    
    optional parameters:
    mu - Earth gravitational constant (m^3/s^2)
    Omega_dot - WGS Earth rotation rate (rad/s)
    compute_velocity - flag indicating whether or not to compute satellite velocity as well (default False)
    """
    e, t_oe, i_0, Omega_dot, a, Omega_week, omega, M_0, week, delta_n, i_dot, c_us, c_rs, c_is, c_uc, c_rc, c_ic = eph
    delta_t = t - gpsseconds(week, t_oe, rollover)
    n_0 = sqrt(mu / a**3)
    n = n_0 + delta_n  # corrected mean motion
    M = M_0 + n * delta_t  # mean anomaly
    E = 0.  # eccentric anomaly
    err = 1.
    tol = 1e-8
    while err > tol:  # iteratively solve for eccentric anomaly
        _temp = M + e * sin(E)
        err = npmax(fabs(_temp - E))
        E = _temp
    nu = arctan2(sqrt(1. - e**2) * sin(E), cos(E) - e)  # true anomaly
    E = arccos((e + cos(nu)) / (1. + e * cos(nu)))  # eccentric anomaly
    Phi = nu + omega  # argument of latitude
    delta_u = c_us * sin(2 * Phi) + c_uc * cos(2 * Phi)  # secord-order harmonic perturbation corrections
    delta_r = c_rs * sin(2 * Phi) + c_rc * cos(2 * Phi)
    delta_i = c_is * sin(2 * Phi) + c_ic * cos(2 * Phi)
    u = Phi + delta_u  # argument of latitude (corrected)
    r = a * (1. - e * cos(E)) + delta_r  # orbital radius (corrected)
    i = i_0 + delta_i + i_dot * delta_t  # inclination (corrected)
    Omega = Omega_week + (Omega_dot - Omega_e_dot) * delta_t - Omega_e_dot * t_oe  # longitude of ascensing node (corrected)
    x_orb, y_orb = r * cos(u), r * sin(u)  # compute x,y in orbital plane
    x_ecef = x_orb * cos(Omega) - y_orb * sin(Omega) * cos(i)  # transform from orbital system to ECEF system
    y_ecef = x_orb * sin(Omega) + y_orb * cos(Omega) * cos(i)
    z_ecef = y_orb * sin(i) * ones(asarray(Omega).size)
    if compute_velocity:  # compute x,y velocity in orbital plane
        v_x_orb = n * a * sin(E) / (1. - e * cos(E))
        v_y_orb = -n * a * sqrt(1. - e**2) * cos(E) / (1. - e * cos(E))
        v_x_ecef = v_x_orb * cos(Omega) - v_y_orb * sin(Omega) * cos(i)  # transform from orbital system to ECEF system
        v_y_ecef = v_x_orb * sin(Omega) + v_y_orb * cos(Omega) * cos(i)
        v_z_ecef = v_y_orb * sin(i) * ones(asarray(Omega).size)
        return asarray([[x_ecef, y_ecef, z_ecef], [v_x_ecef, v_y_ecef, v_z_ecef]])
    else:
        return asarray([x_ecef, y_ecef, z_ecef])


def compute_gps_orbital_parameters_from_ephemeris(eph, t, mu=3.986005e14, Omega_e_dot=7.2921151467e-5, rollover=0):
    """Computes and returns the following orbital parameters that define GPS satellite position/velocity:
    argument of latitude
    orbital radius
    orbital inclination
    longitude of ascensing node
    corrected mean motion
    eccentric anomaly
    
    parameters:
    eph - the ephemeris set (an ordered tuple containing ephemeris information)
        `eph` should contain: e, t_oa, i_0, Omega_dot, sqrt_a, Omega_week, omega, M_0, A_f0, A_f1, week  TODODODO
    t - the desired times for which to compute satellite locations
    
    optional parameters:
    mu - Earth gravitational constant (m^3/s^2)
    Omega_dot - WGS Earth rotation rate (rad/s)
    rollover - rollover to include when computing GPS time from week number (default 0)
    """
    e, t_oe, i_0, a, omega_dot, omega_0, omega, M_0, week, delta_n, i_dot, c_us, c_rs, c_is, c_uc, c_rc, c_ic = \
            eph.e, eph.t_oe, eph.i_0, eph.a, eph.omega_dot, eph.omega_0, eph.omega, eph.m_0, eph.week, eph.delta_n, eph.i_dot, \
            eph.c_us, eph.c_rs, eph.c_is, eph.c_uc, eph.c_rc, eph.c_ic
    delta_t = t - gpsseconds(week, t_oe, rollover)
    n_0 = sqrt(mu / a**3)
    n = n_0 + delta_n  # corrected mean motion
    M = M_0 + n * delta_t  # mean anomaly
    E = 0.  # eccentric anomaly
    err = 1.
    tol = 1e-8
    while err > tol:  # iteratively solve for eccentric anomaly
        _temp = M + e * sin(E)
        err = npmax(fabs(_temp - E))
        E = _temp
    nu = arctan2(sqrt(1. - e**2) * sin(E), cos(E) - e)  # true anomaly
    E = arccos((e + cos(nu)) / (1. + e * cos(nu)))  # eccentric anomaly
    Phi = nu + omega  # argument of latitude
    delta_u = c_us * sin(2 * Phi) + c_uc * cos(2 * Phi)  # secord-order harmonic perturbation corrections
    delta_r = c_rs * sin(2 * Phi) + c_rc * cos(2 * Phi)
    delta_i = c_is * sin(2 * Phi) + c_ic * cos(2 * Phi)
    u = Phi + delta_u  # argument of latitude (corrected)
    r = a * (1. - e * cos(E)) + delta_r  # orbital radius (corrected)
    i = i_0 + delta_i + i_dot * delta_t  # inclination (corrected)
    Omega = omega_0 + (omega_dot - Omega_e_dot) * delta_t - Omega_e_dot * t_oe  # longitude of ascensing node (corrected)
    return u, r, i, Omega, n, E

def compute_gps_satellite_position_from_orbital_parameters(u, r, i, Omega):
    '''Computes and returns the satellite position given a set of orbital parameters
    that define the satellite position
    
    parameters:
    u - argument of latitude
    r - orbital radius
    i - orbital inclination
    Omega - longitude of ascending node
    '''
    x_orb, y_orb = r * cos(u), r * sin(u)  # compute x,y in orbital plane
    x_ecef = x_orb * cos(Omega) - y_orb * sin(Omega) * cos(i)  # transform from orbital system to ECEF system
    y_ecef = x_orb * sin(Omega) + y_orb * cos(Omega) * cos(i)
    z_ecef = y_orb * sin(i) * ones(asarray(Omega).size)
    return asarray([x_ecef, y_ecef, z_ecef])

def compute_gps_satellite_velocity_from_orbital_parameters(e, i, Omega, n, E):
    '''Computes and returns the satellite velocity given a set of orbital parameters
    that define the satellite velocity.
    
    parameters:
    e - orbital eccentricity
    i - orbital inclination
    Omega - longitude of ascending node
    n - mean motion
    E - eccentric anomaly
    '''
    v_x_orb = n * a * sin(E) / (1. - e * cos(E))
    v_y_orb = -n * a * sqrt(1. - e**2) * cos(E) / (1. - e * cos(E))
    v_x_ecef = v_x_orb * cos(Omega) - v_y_orb * sin(Omega) * cos(i)  # transform from orbital system to ECEF system
    v_y_ecef = v_x_orb * sin(Omega) + v_y_orb * cos(Omega) * cos(i)
    v_z_ecef = v_y_orb * sin(i) * ones(asarray(Omega).size)
    return asarray([v_x_ecef, v_y_ecef, v_z_ecef])


