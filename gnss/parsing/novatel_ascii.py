#!/usr/bin/env python

"""
Reads data from NovAtel ASCII formatted text file into a buffer.
See here for file format information: http://www.novatel.com/assets/Documents/Manuals/om-20000047.pdf
"""

__author__ = "Brian Breitsch"
__copyright__ = "Copyright 2013"
__credits__ = ["Brian Breitsch"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Brian Breitsch"
__email__ = "brianbreitsch@gmail.com"
__status__ = "Infant"

from collections import namedtuple
from numpy import uint32, float64, nan, zeros
from gnss.time import gpsseconds

def count_novatel_ascii_logs(_file):
    """
    Determines the number of log entries in the NovAtel ASCII formatted file at
    filepath.
    Returns a dictionary with message strings as keys and count as values.
    """
    # the following messages are those expected to be in the file
    messages = {
                'RANGEA':0,
                'BESTPOSA':0,
                'GPSEPHEMA':0,
                'IONUTCA':0,
               }
    with open(_file, 'r') as f:
        for line in f:
            message = line[1:line.find(',')]
            messages[message] += 1
    return messages


class ChannelObservation():

    def __init__(self, adr=None, psr=None, cnr=None):
        self.adr = adr
        self.psr = psr
        self.cnr = cnr


class Observation():

    def __init__(self, time=None, l1=None, l2=None):
        self.time = time
        self.l1 = l1 if l1 else ChannelObservation()
        self.l2 = l2 if l2 else ChannelObservation()


class BestPos():

    def __init__(self, time=None, geo=None):
        self.time = time
        self.geo = geo


class Ephemeris():

    def __init__(self, prn=None, tow=None, health=None, iode1=None, iode2=None, week=None, zweek=None, t_oe=None, a=None, delta_n=None, m_0=None,
            e=None, omega=None, c_uc=None, c_us=None, c_rc=None, c_rs=None, c_ic=None, c_is=None, i_0=None, i_dot=None, omega_0=None, omega_dot=None, 
            iodc=None, t_oc=None, t_gd=None, a0=None, a1=None, a2=None, n=None):
        self.prn = prn
        self.tow = tow
        self.health = health
        self.iode1 = iode1
        self.iode2 = iode2
        self.week = week
        self.zweek = zweek
        self.t_oe = t_oe
        self.a = a
        self.delta_n = delta_n
        self.m_0 = m_0
        self.e = e
        self.omega = omega
        self.c_uc = c_uc
        self.c_us = c_us
        self.c_rc = c_rc
        self.c_rs = c_rs
        self.c_ic = c_ic
        self.c_is = c_is
        self.i_0 = i_0
        self.i_dot = i_dot
        self.omega_0 = omega_0
        self.omega_dot = omega_dot
        self.iodc = iodc
        self.t_oc = t_oc
        self.t_gd = t_gd
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.n = n
        

def read_novatel_ascii(_file, svids=range(1, 33)):
    """
    Reads gps log data from file at filepath. File should be NovAtel ASCII
    formatted text file. Returns namedtuple structs containing all file data.
    `obs, eph, bestpos`
    """
    message_sizes = count_novatel_ascii_logs(_file)
    POS_SIZE = message_sizes['BESTPOSA']
    EPH_SIZE = message_sizes['GPSEPHEMA']
    RNG_SIZE = message_sizes['RANGEA']
    NUM_FIELDS = 10

    obs = {}
    for svid in svids:
        obs[svid] = Observation()
        obs[svid].time = zeros((RNG_SIZE,))
        obs[svid].l1 = ChannelObservation()
        obs[svid].l2 = ChannelObservation()
        obs[svid].l1.adr = zeros((RNG_SIZE,))
        obs[svid].l1.dopp= zeros((RNG_SIZE,))
        obs[svid].l1.psr = zeros((RNG_SIZE,))
        obs[svid].l1.cnr = zeros((RNG_SIZE,))
        obs[svid].l2.adr = zeros((RNG_SIZE,))
        obs[svid].l2.dopp= zeros((RNG_SIZE,))
        obs[svid].l2.psr = zeros((RNG_SIZE,))

        obs[svid].time.fill(nan)
        obs[svid].l1.adr.fill(nan)
        obs[svid].l1.dopp.fill(nan)
        obs[svid].l1.psr.fill(nan)
        obs[svid].l1.cnr.fill(nan)
        obs[svid].l2.adr.fill(nan)
        obs[svid].l2.dopp.fill(nan)
        obs[svid].l2.psr.fill(nan)

    eph = {}
    for svid in svids:
        eph[svid] = Ephemeris()

    bestpos = BestPos()
    bestpos.time = zeros((POS_SIZE,))
    bestpos.geo = zeros((POS_SIZE, 3))

    pos_ind = 0
    rng_ind = 0
   
    with open(_file, 'r') as f:
        for line in f:
            line = line.strip().split(';')
            header, data = line[0].split(','), line[1].split('*')  # crc after *
            data, crc = data[0].split(','), data[1]
            week_no, seconds = int(header[5]), float(header[6])
            gpstime = gpsseconds(week_no, seconds)  # no rollover for NovAtel
            if header[0] == '#BESTPOSA':
                bestpos.time[pos_ind] = gpstime
                bestpos.geo[pos_ind, :] = [float(data[2]), float(data[3]), float(data[4])]
                pos_ind += 1
            elif header[0] == '#GPSEPHEMA':
                prn             = uint32( data[0])  # Satellite PRN
                eph[prn].prn    = prn
                eph[prn].tow    = float64(data[1])  # time of week
                eph[prn].health = uint32( data[2])  # satellite health
                eph[prn].iode1  = uint32( data[3])  # issue of ephemeris data
                eph[prn].iode2  = uint32( data[4])  # issue of ephemeris data 2
                eph[prn].week   = uint32( data[5])  # GPS week number
                eph[prn].zweek  = uint32( data[6])  # GPS zweek number
                eph[prn].t_oe   = float64(data[7])  # time of applicability for ephemeris (s)
                eph[prn].a      = float64(data[8])  # semi-major axis (m)
                eph[prn].delta_n = float64(data[9]) # mean-motion diff (rad/s)
                eph[prn].m_0    = float64(data[10]) # mean anomoly of ref time (rad)
                eph[prn].e      = float64(data[11]) # eccentricity
                eph[prn].omega  = float64(data[12]) # argument of perigee
                eph[prn].c_uc   = float64(data[13]) # argument of latitude (amplitude of cosine, radians)
                eph[prn].c_us   = float64(data[14]) # argument of latitude (amplitude of sine, radians)
                eph[prn].c_rc   = float64(data[15]) # orbit radius (amplitude of cosine, meters)
                eph[prn].c_rs   = float64(data[16]) # orbit radius (amplitude of sine, meters)
                eph[prn].c_ic   = float64(data[17]) # inclination (amplitude of cosine, meters)
                eph[prn].c_is   = float64(data[18]) # inclination (amplitude of sine, meters)
                eph[prn].i_0    = float64(data[19]) # inclination angle at reference time, (rad)
                eph[prn].i_dot  = float64(data[20]) # rate of inclination angle, (rad/s)
                eph[prn].omega_0 = float64(data[21]) # right ascension at week (rad)
                eph[prn].omega_dot  = float64(data[22]) # rate of right ascension, (rad/s)
                eph[prn].iodc   = uint32( data[23]) # issue of data clock
                eph[prn].t_oc   = float64(data[24]) # SV clock correction term, seconds
                eph[prn].t_gd   = float64(data[25]) # estimated group delay difference, seconds
                eph[prn].a0     = float64(data[26]) # clock aging parameter, seconds (s)
                eph[prn].a1     = float64(data[27]) # clock aging parameter, (s/s)
                eph[prn].a2     = float64(data[28]) # clock aging parameter, (s/s/s)
                eph[prn].n      = float64(data[30]) # corrected mean motion, (rad/s)
            elif header[0] == '#RANGEA':
                num_records = int(data[0])
                data = data[1:]
                for i in range(0, num_records * NUM_FIELDS, NUM_FIELDS):
                    # check channel tracking status (NovAtel Manuel p 237)
                    prn = int(data[i])
                    obs[prn].time[rng_ind] = gpstime
                    tracking_status = int(data[i + 9], 16)
                    if not tracking_status & (0x1 << 21):  # L1
                        obs[prn].l1.psr[rng_ind] = float(data[i + 2])
                        obs[prn].l1.adr[rng_ind] = float(data[i + 4])
                        obs[prn].l1.dopp[rng_ind]= float(data[i + 6])
                        obs[prn].l1.cnr[rng_ind] = float(data[i + 7])
                    else:  # L2
                        obs[prn].l2.psr[rng_ind] = float(data[i + 2])
                        obs[prn].l2.dopp[rng_ind]= float(data[i + 6])
                        obs[prn].l2.adr[rng_ind] = float(data[i + 4])
                rng_ind += 1
    return obs, eph, bestpos
