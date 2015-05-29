from collections import namedtuple
from numpy import ones, zeros, cos, sin, radians, asarray, pi
import numpy
npsum = numpy.sum
from bokeh.plotting import figure, cursession


SkyCoordinates = namedtuple('SkyCoordinates', 'azimuth elevation')


class SkyPlot:
    
    def __init__(self, title='sky plot'):
        r = ones((1,))
        theta = zeros((1,))
        self.plot = figure(title=title, plot_width=500, plot_height=500, x_range=[-90, 90], y_range=[-90, 90])
        self.plot.grid.grid_line_color = None
        self.plot.axis.axis_line_color = None
        self.plot.axis.major_tick_line_color = None
        x = [0, 0, 1, 1]
        y = [0, 1, 1, 0]
        self.plot.scatter([x], [y], size=12, alpha=0.7, name='points')
        svid_markers = []
        # todo plot markers

        # data is a dictionary of az/el coordinates for a given svid
        self.data = {}
        
    
    def update(self, svid, azimuth, elevation):
        self.data[svid] = SkyCoordinates(azimuth, elevation)
        
        azimuth = asarray([sky.azimuth for sky in self.data.values()]).flatten()
        elevation = asarray([sky.elevation for sky in self.data.values()]).flatten()

        r = 90 - elevation
        theta = pi / 2 - radians(azimuth)
        mask = r < 90
        print(azimuth.shape, elevation.shape, npsum(mask))
        r = r[mask]
        theta = theta[mask]
        x = r * cos(theta)
        y = r * sin(theta)
        print(x.shape, y.shape)
        renderer = self.plot.select(dict(name='points'))[0]
        renderer.data_source.data['x'] = x
        renderer.data_source.data['y'] = y
        cursession().store_objects(renderer.data_source)
