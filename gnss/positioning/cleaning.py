



class 







    '''Returns cleaned pseudorange measurement given GPS L1 and L2 measurements of code range and carrier phase (cycles).'''

    # estimate propagation time
    tau = l1_psr / c
    # estimate transmission time
    sv_time = rx_time - tau  # => + 72m error w/o
    
    # compute satellite position
    u, r, i, Omega, n, E = compute_gps_orbital_parameters_from_ephemeris(eph, sv_time)
    sv_ecef[svid] = compute_gps_satellite_position_from_orbital_parameters(u, r, i, Omega)
    
    # compute satellite clock drift error and relativistic corrections
    sv_time_correction_epoch = eph.zweek * 604800 + eph.t_oc
    elapsed_time = sv_time - sv_time_correction_epoch
    sv_time_diff = eph.a0 + eph.a1 * elapsed_time + eph.a2 * elapsed_time**2 \
                    - 4.442807633e-10 * eph.e * sqrt(eph.a) * sin(E)

    # compute satellite sky coordinates to get obliquity
    sv_sky = ecef2sky(rx_ecef_ref, sv_ecef[svid])
    el = radians(sv_sky[1])
    
    # correct ECEF for Earth rotation
    theta = omega_e_dot * tau
    rot_omega_e = asarray([[cos(theta), sin(theta), 0],
                            [-sin(theta), cos(theta), 0],
                            [0, 0, 1]])
    sv_ecef[svid] = rot_omega_e.dot(sv_ecef[svid].T).T  # <-- 300m satellite pos error w/o => +20m rx pos error

def clean_pseudorange_gps_l1_l2(l1_psr, l2_psr, l1_adr, l2_adr, tropo_model, iono_model):
    # compute partially cleaned carrier phase measurements (for iono error--better than psr, clock stuff doesn't matter)
    adr_l1_c = l1_adr + c * sv_time_diff / lambda_l1
    adr_l2_c = l2_adr + c * sv_time_diff / lambda_l2
    
    # compute tec and ionospheric/tropospheric delays
    tropo_delay = tropo_model.delay(el)
    iono_delay = IonosphereDelayModel.compute_delay(f_carr_l1, f_carr_l2, obs.l1.psr[page], obs.l2.psr[page])

    # TODO mask erroneous delay error

    # compute cleaned pseudorange; sv_time_diff => +160km pos error, tropo_error => 15m error, iono -2m error
    clean_psr = l1_psr + c * sv_time_diff - tropo_delay # - iono_delay
