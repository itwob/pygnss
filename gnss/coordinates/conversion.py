"""
coordinate_utils.py

@author Brian Breitsch
@email brianbreitsch@gmail.com
"""

from numpy import sin, cos, tan, array, zeros, sqrt, arctan2, rad2deg, deg2rad, column_stack, absolute, arcsin, pi, asarray
from datetime import datetime, timedelta
from pytz import UTC
import numpy

def make_3_tuple_array(arr):
    """Reshapes ndarray so that it has dimensions (N,3)
    """
    if type(arr) == list:
        arr = array(arr)
    assert(arr.ndim <= 2 and arr.size >= 3)
    if arr.shape[0] == 3:
        arr = arr.reshape((1,3)) if arr.ndim == 1 else arr.T
    assert(arr.shape[1] == 3)
    return arr


def ecef2geo(x_ref, N=None):
    """Converts ECEF coordinates to geodetic coordinates, 

    Parameters
    ----------
    x_ref : an ndarray of N ECEF coordinates with shape (N,3).

    Returns
    -------
    output : (N,3) ndarray
        geodetic coordinates
    `N` : (optional) the radius of curvature, passed in as output parameter

    Notes
    -----
    >>> from numpy import array, deg2rad
    >>> geo = array([27.174167, 78.042222, 0])
    >>> ecef = geo2ecef(deg2rad(geo))
    >>> new_geo = ecef2geo(ecef)
    array([[             nan],
    	   [  7.08019709e+01],
           [ -6.37805436e+06]])
    >>> # [1176.45, 5554.887, 2895.397] << should be this 
    >>> ecef2geo(array([
    	[27.174167, 78.042222, 0],
    	[39.5075, -84.746667, 0]])).reshaps((3,2))
    array([[             nan,              nan],
           [  7.08019709e+01,  -6.50058423e+01],
           [ -6.37805436e+06,  -6.37804350e+06]])
    [1176.45, 5554.887, 2895.397]
    [451.176, -4906.978, 4035.946]
    """
    x_ref = make_3_tuple_array(x_ref)
    # we = 7292115*1**-11 # Earth angular velocity (rad/s)
    # c = 299792458    # Speed of light in vacuum (m/s)
    rf = 298.257223563 # Reciprocal flattening (1/f)
    a = 6378137.       # Earth semi-major axis (m)
    b = a - a / rf     # Earth semi-minor axis derived from f = (a - b) / a
    if x_ref.shape == (3,):
        x_ref = x_ref.reshape((1,3))
    x = x_ref[:,0]; y = x_ref[:,1]; z = x_ref[:,2];

    # We must iteratively derive N
    lat = arctan2(z, sqrt(x**2 + y**2))
    h = z / sin(lat)
    d_h = 1.; d_lat = 1.

    if numpy.any(N) is None:
        N = zeros(x_ref.shape[0])

    while (d_h > 1e-10) and (d_lat > 1e-10):
    
        N[:] = a**2 / (sqrt(a**2 * cos(lat)**2 + b**2 * sin(lat)**2))
        N1 = N * (b / a)**2
    
        temp_h = sqrt(x**2 + y**2) / cos(lat) - N
        temp_lat = arctan2(z / (N1 + h), sqrt(x**2 + y**2) / (N + h))
        d_h = numpy.max(absolute(h - temp_h))
        d_lat = numpy.max(absolute(lat - temp_lat))

        h = temp_h
        lat = temp_lat

    lon = arctan2(y,x)

    lat = rad2deg(lat)
    lon = rad2deg(lon)

    geo = column_stack((lat, lon, h))
    return geo.squeeze()


def ecef2enu(x_ref, x_obj):
    """Converts satellite ECEF coordinates to user-relative ENU coordinates.

    Parameters
    ----------
    x_ref : ndarray of shape (3,)
        observer coordinate
    x_obj : ndarray of shape(N,3)
        object coordinates

    Returns
    -------
    output : ndarray of shape(N,3)
        The east-north-up coordinates

    Notes
    -----
    """
    x_ref = make_3_tuple_array(x_ref)
    x_obj = make_3_tuple_array(x_obj)
    # get the lat and lon of the user position
    geo = ecef2geo(x_ref)
    if geo.shape != (3,):
        raise Exception('ecef2enu can only handle one reference position')
    lat = deg2rad(geo[0]); lon = deg2rad(geo[1])

    # create the rotation matrix
    Rl = array([[-sin(lon),               cos(lon),              0],
          [-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat)],
          [cos(lat) * cos(lon),  cos(lat) * sin(lon), sin(lat)]])
    # NOTE: This matrix is not fixed to do multiple user locations yet...

    dx = x_obj - x_ref
    return Rl.dot(dx.T).T.squeeze()


def ecef2sky(x_ref, x_obj):
    """Converts user and satellite ecef coordinates to azimuth and elevation
    from user on Earth.

    Parameters
    ----------
    x_ref : ndarray of shape (3,)
        observer coordinate
    xs : ndarray of shape(N,3)
        object coordinates

    Returns
    -------
    output : ndarray of shape(N,2)
        The objects' sky coordinatescoordinates
        (azimuth, elevation in degrees)

    Notes
    -----
    """

    enu = ecef2enu(x_ref, x_obj)
    enu = make_3_tuple_array(enu)
    e = enu[:, 0]; n = enu[:, 1]; u = enu[:, 2]
    az = arctan2(e, n)
    el = arcsin(u / sqrt(e**2 + n**2 + u**2))

    return column_stack((rad2deg(az), rad2deg(el))).squeeze()


def geo2ecef(geo):
    """Converts geodetic coordinates to ECEF coordinates

    Parameters
    ----------
    geo : ndarray of shape (N,3)
        geodetic coordinates

    Returns
    -------
    output : ndarray of shape(N,3)
        ECEF coordinates

    Notes
    -----
    """
    geo = make_3_tuple_array(geo)
    a = 6378137. # Earth semi-major axis (m)
    rf = 298.257223563 # Reciprocal flattening (1/f)
    b = a * (rf - 1) / rf # Earth semi-minor axis derived from f = (a - b) / a
    lat = deg2rad(geo[:, 0]); lon = deg2rad(geo[:, 1]); h = geo[:, 2]

    N = a**2 / sqrt(a**2 * cos(lat)**2 + b**2 * sin(lat)**2)
    N1 = N * (b / a)**2

    x = (N + h) * cos(lat) * cos(lon)
    y = (N + h) * cos(lat) * sin(lon)
    z = (N1 + h) * sin(lat)

    x_ref = column_stack((x, y, z))
    return x_ref.squeeze()


def geo2sky(geo_ref, geo_obj):
    """Converts object geodetic coordinates to azimuth and elevation from
    reference geodetic coordinates on Earth.

    Parameters
    ----------
    geo_ref : ndarray of shape (3,)
        geodetic (lat, lon, alt) coordinates of observer
    geo_obj : ndarray of shape (N,3)
        geodetic (lat, lon, alt) coordinates of object

    Returns
    -------
    output : ndarray of shape(N,2)
        sky coordinates (azimuth, elevation in degrees)

    Notes
    -----
    """
    geo_ref = make_3_tuple_array(geo_ref)
    geo_obj = make_3_tuple_array(geo_obj)
    x_ref = geo2ecef(geo_ref)
    x_obj = geo2ecef(geo_obj)

    sky = ecef2sky(x_ref, x_obj)
    return sky


def eci2ecef(eci, time):
    """Converts eci coordinates to ecef coordinates given
    the time or times for said coordinates

    Parameters
    ----------
    eci : ndarray of shape (N,3)
    time : time array in GPS seconds of shape (N,)

    Returns
    -------
    output : ndarray of shape (N,3)
        ecef coordinates

    Notes
    Time of 2000 January 1, 12 UTC is (in GPS seconds) 630763213.0
    -----
    See: http://physics.stackexchange.com/questions/98466/radians-to-rotate-earth-to-match-eci-lat-lon-with-ecef-lat-lon

    TODO: only accepts one time right now
    """
    #delta_days = (time - datetime64(datetime(2000, 1, 1, 12, tzinfo=UTC))) / timedelta64(1, 'D')
    delta_days = (time - 630763213.0) / (3600 * 24)
    gmst = 18.697374558 + 24.065709824419 * delta_days
    theta = 2 * pi / 24 * gmst
    rot = asarray([[ cos(gmst), sin(gmst), 0],
                   [-sin(gmst), cos(gmst), 0],
                   [ 0,         0,         1]])
    return rot.dot(eci.T).T


def ecef2eci(ecef, time):
    """Converts ecef coordinates to eci coordinates given
    the time or times for said coordinates

    Parameters
    ----------
    ecef : ndarray of shape (N,3)
    time : time array in GPS seconds of shape (N,)

    Returns
    -------
    output : ndarray of shape (N,3)
        eci coordinates

    Notes
    -----
    See: http://physics.stackexchange.com/questions/98466/radians-to-rotate-earth-to-match-eci-lat-lon-with-ecef-lat-lon
    Time of 2000 January 1, 12 UTC is (in GPS seconds) 630763213.0

    TODO: only accepts one time right now
    """
    delta_days = (time - 630763213.0) / (3600 * 24)
    gmst = 18.697374558 + 24.065709824419 * delta_days
    theta = 2 * pi / 24 * gmst
    rot = asarray([[cos(gmst), -sin(gmst), 0],
                   [sin(gmst),  cos(gmst), 0],
                   [ 0,         0,         1]])
    return rot.dot(ecef.T).T
