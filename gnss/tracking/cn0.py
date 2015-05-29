from numpy import sqrt, mean, var, log10



class CN0Algorithm(object):
    '''
    Encapsulates algorithm to compute carrier-to-noise
    ratio.
    '''

    def __init__(self):
        pass

    def compute(self, i_corr_data, q_corr_data, integration_time):
        '''
        Computes carrier-to-noise ratio given I and Q correlation samples
        and the integration time and accumulation times.
        `integration_time` -- the time of coherent correlation used to generate 
            `i_corr_data` and `q_corr_data` samples.
        '''
        z = i_corr_data**2 + q_corr_data**2
        z_ave = mean(z)
        z_var = var(z)
        ave_power = sqrt(z_ave**2 - z_var)
        var_iq = 0.5 * (z_ave - ave_power)
        return 10 * log10(1. / integration_time * ave_power / (2. * var_iq))


