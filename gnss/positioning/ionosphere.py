

class IonosphereDelayModel(object):
    '''
    Ionosphere delay model.
    '''
    c = 299792458.  # Speed of light in vacuum (m/s)
    
    def __init__(self, carr_freq_1, carr_freq_2):
        '''
        `carr_freq_1` and `carr_freq_2` -- carrier phases of two signals to be used for model delay computation
        '''
        self.carr_freq_1 = carr_freq_1
        self.carr_freq_2 = carr_freq_2

    def compute_tec_from_psr(self, carr_freq_1, carr_freq_2, psr_1, psr_2):
        '''Computes total electron content given
        `carr_freq_1` and `carr_freq_2` -- the two carrier frequencies
        `psr_1` and `psr_2` -- the two pseudorange measurements'''
        return carr_freq_1**2 * carr_freq_2**2 / (40.3 * (carr_freq_1**2 - carr_freq_2**2)) * (psr_2 - psr_1)

    def compute_relative_tec_from_phase(self, carr_freq_1, carr_freq_2, phi_1, phi_2):
        '''Computes relative total electron content given
        `carr_freq_1` and `carr_freq_2` -- the two carrier frequencies
        `phi_1` and `phi_2` -- the two phase measurements'''
        c = IonosphereDelayModel.c  # speed of light
        lambda_1 = c / carr_freq_1
        lambda_2 = c / carr_freq_2
        return carr_freq_1**2 * carr_freq_2**2 / (40.3 * (carr_freq_1**2 - carr_freq_2**2)) * (lambda_1 * phi_1 - lambda_2 * phi_2)

    def compute_delay(self, carr_freq_1, carr_freq_2, psr_1, psr_2):
        '''Computes ionospheric delay for signal 1 given
        `psr_1` and `psr_2` -- pseudorange measurements corresponding to signals at frequencies
        `carr_freq_1` and `carr_freq_2` -- Hz
        '''
        tec = self.compute_tec_from_psr(carr_freq_1, carr_freq_2, psr_1, psr_2)
        return 40.3 * tec / carr_freq_1**2

    def compute_delay_from_phase(self, carr_freq_1, carr_freq_2, phi_1, phi_2):
        '''Computes ionospheric delay for signal 1 given
        `phi_1` and `phi_2` -- carrier measurements corresponding to signals at frequencies
        `carr_freq_1` and `carr_freq_2` -- Hz
        '''
        tec = self.compute_relative_tec_from_phase(carr_freq_1, carr_freq_2, phi_1, phi_2)
        return 40.3 * tec / carr_freq_1**2
        
    def delay(self, psr_1, psr_2):
        '''Computes delay from two pseudorange values using `self.compute_delay`'''
        return self.compute_delay(self.carr_freq_1, self.carr_freq_2, psr_1, psr_2)
