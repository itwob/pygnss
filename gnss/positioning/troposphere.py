from numpy import exp, sin, cos, tan, radians


class TroposphereDelayModel(object):

    dry_boundary_altitude = 43e3  # km
    wet_boundary_altitude = 12e3  # km

    @classmethod
    def h20_partial_pressure_from_relative_humidity(cls, relative_humidity, temperature):
        '''Compute the partial pressure of water vapor from the relative humidity'''
        return 6.108 * relative_humidity * exp((17.15 * temperature - 4684) / (temperature - 38.45))

    def compute_dry_refractivity(pressure, temperature, altitude):
        '''Compute dry refractivity given
        `pressure` -- total pressure (mbar)
        `temperature` -- in Kelvin
        `altitude` -- meters (or same scale as `dry_boundary_altitude`)'''
        surface_dry_refractivity = 77.64 * pressure / temperature
        return surface_dry_refractivity * (1. - altitude / dry_boundary_altitude)**4

    def compute_wet_refractivity(temperature, h20_partial_pressure, altitude):
        '''Compute wet refractivity given
        `temperature` -- in Kelvin
        `h20_partial_pressure` -- water vapor partial pressure (mbar)
        `altitude` -- meters (or same scale as `wet_boundary_altitude`)'''
        surface_wet_refractivity = 3.73e5 * h20_partial_pressure / temperature**2
        return surface_wet_refractivity * (1. - altitude / dry_boundary_altitude)**4

    def compute_delay(elevation, dry_zenith_delay, wet_zenith_delay):
        '''Computes tropospheric delay given satellite elevation and
        dry/wet zenith delay estimates'''
        # the following are empricial formulas for obliquity factor for dry/wet atmosphere contribution
        obliq_dry = 1. / (sin(elevation) + 1.43e-3 / (tan(elevation) + 44.5e-3))
        obliq_wet = 1. / (sin(elevation) + 0.35e-3 / (tan(elevation) + 17.0e-3))
        return obliq_dry * dry_zenith_delay + obliq_wet * wet_zenith_delay

    def __init__(self, pressure, h20_partial_pressure, temperature, altitude=0):
        '''
        `pressure` -- total pressure (mbar)
        `h20_partial_pressure` -- water vapor partial pressure (mbar)
        `temperature` -- in Kelvin!!!
        `altitude` -- altitude (m) of measurements? I think? TODO
        '''
        self.pressure = pressure
        self.h20_partial_pressure = h20_partial_pressure
        self.temperature = temperature
        self.altitude = altitude

    @property
    def dry_refractivity(self):
        '''Computes dry refractivity
        See: `compute_dry_refractivity`'''
        return TroposphereDelayModel.compute_dry_refractivity(self.pressure, self.temperature, self.altitude)

    @property
    def wet_refractivity(self):
        '''Computes wet refractivity
        See: `compute_wet_refractivity`'''
        return TroposphereDelayModel.compute_wet_refractivity(self.temperature, self.h20_partial_pressure, self.altitude)


class SaastomoinenModel(TroposphereDelayModel):
    
    def compute_dry_zenith_delay(latitude, altitude):
        '''Compute dry component of tropospheric zenith delay given
        `latitude` -- in degrees
        `altitude` -- in meters'''
        return 2.277e-3 * (1. + 2.6e-3 * cos(2 * radians(latitude)) + 0.28e-6 * altitude)
    
    def compute_wet_zenith_delay(temperature, h20_partial_pressure):
        '''Compute wet component of tropospheric zenith delay given
        `temperature` -- degrees Kelvin
        `h20_partial_pressure -- water vapor partial pressure (mbar)'''
        return 2.2277e-3 * (1255 / temperature + 0.05) * h20_partial_pressure

    def __init__(self, pressure, h20_partial_pressure, temperature, latitude, altitude=0):
        super().__init__(pressure, h20_partial_pressure, temperature, altitude)
        self.latitude = latitude

    @property
    def dry_zenith_delay(self):
        '''Computes dry zenith delay
        See: `compute_dry_zenith_delay`'''
        return SaastomoinenModel.compute_dry_zenith_delay(self.latitude, self.altitude)

    @property
    def wet_zenith_delay(self):
        '''Computes wet zenith delay
        See: `compute_wet_zenith_delay`'''
        return SaastomoinenModel.compute_wet_zenith_delay(self.temperature, self.h20_partial_pressure)

    def delay(self, elevation):
        '''Computes total tropospheric delay given satellite elevation.
        See `compute_delay`'''
        return TroposphereDelayModel.compute_delay(elevation, self.dry_zenith_delay, self.wet_zenith_delay)


class HopfieldModel(TroposphereDelayModel):

    def compute_dry_zenith_delay(pressure, temperature, dry_boundary_altitude=TroposphereDelayModel.dry_boundary_altitude):
        '''Compute dry component of tropospheric zenith delay given
        `pressure` -- total pressure in millibars (i.e. mbar)
        `temperature` -- in Kelvin
        `dry_boundary_altitude` -- (optional) boundary altitude of dry atmosphere'''
        return 77.6e-6 * pressure / temperature * dry_boundary_altitude / 5.
    
    def compute_wet_zenith_delay(temperature, h20_partial_pressure, wet_boundary_altitude=TroposphereDelayModel.wet_boundary_altitude):
        '''Compute wet component of tropospheric zenith delay given
        `temperature` -- in Kelvin
        `h20_partial_pressure` -- partial pressure of water vapor (mbar)
        `wet_boundary_altitude` -- (optional) boundary of wet atmosphere'''
        return 0.373 * h20_partial_pressure / temperature**2 * wet_boundary_altitude / 5.

    @property
    def dry_zenith_delay(self):
        '''Computes dry zenith delay
        See: `compute_dry_zenith_delay`'''
        return HopfieldModel.compute_dry_zenith_delay(self.pressure, self.temperature, self.dry_boundary_altitude)
    
    @property
    def wet_zenith_delay(self):
        '''Computes wet zenith delay
        See: `compute_wet_zenith_delay`'''
        return HopfieldModel.compute_wet_zenith_delay(self.temperature, self.h20_partial_pressure, self.wet_boundary_altitude)

    def delay(self, elevation):
        '''Computes total tropospheric delay given satellite elevation.
        See `compute_delay`'''
        return TroposphereDelayModel.compute_delay(elevation, self.dry_zenith_delay, self.wet_zenith_delay)




