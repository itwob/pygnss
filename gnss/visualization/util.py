from numpy import ones, zeros, cos, sin, real, imag, histogram2d, meshgrid, any, iscomplex
import numpy
npmax = numpy.max
from bokeh.plotting import figure, cursession

class PolarPlot:
    
    def __init__(self, title='polar plot', r=ones((1,)), theta=zeros((1,))):
        self.plot = figure(title=title, plot_width=500, plot_height=500)
        x = r * cos(theta)
        y = r * sin(theta)
        self.plot.scatter(x, y, size=12, alpha=0.7, name='scatter')
    
    def update(self, r, theta):
        renderer = self.plot.select(dict(name='scatter'))[0]
        x = r * cos(theta)
        y = r * sin(theta)
        renderer.data_source.data['x'] = x
        renderer.data_source.data['y'] = y
        cursession().store_objects(renderer.data_source)


def sample_histogram(samples):
    '''Returns a bokeh plot object histogram of samples.
    If samples are complex, the plot is a 2D scatter plot
    with marker sizes indicating magnitude.
    If samples are real, the plot is a 1D histogram.
    NOTE: currently assumes complex are 8-bit complex samples
    '''
    if any(iscomplex(samples)):
        r = real(samples)
        i = imag(samples)
        magn, x_bins, y_bins = histogram2d(r, i, 16, range=((-8,8), (-8,8)))
        x, y = meshgrid(x_bins, y_bins)
        sizes = magn / npmax(magn) * .7
        x = x[:-1, :-1].flatten()
        y = y[:-1, :-1].flatten()
        sizes = sizes.flatten()
        f = figure(title='complex samples scatter histogram')
        f.scatter(x, y, radius=sizes)
        return f
    else:
        raise Exception('not implemented')

