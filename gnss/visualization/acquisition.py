
import numpy
npmax = numpy.max
from numpy import ones, uint8, uint32, arange, angle, unwrap, pi
from bokeh.plotting import hplot, figure, show, cursession
from gnss.acquisition import CoarseAcquirer, CoarseAcquirerLowMem, FineAcquirer

def plot_coarse_acquisition_results(acq):
    if isinstance(acq, CoarseAcquirerLowMem):
        dopp_min, dopp_max = acq.dopp_bins[0], acq.dopp_bins[-1]
        p = figure(title='correlation', x_range=[dopp_min, dopp_max], 
                x_axis_label='doppler (Hz)', y_axis_label='correlation strength',
                plot_width=500, plot_height=500)
        p.line(acq.dopp_bins, acq.corr_max)
        return p
    elif isinstance(acq, CoarseAcquirer):
        c1, c2 = acq.plot_code_window
        r, c = acq.plot_corr.shape
        img = ones((r, c), dtype=uint32)
        view = img.view(dtype=uint8).reshape((r, c, 4))
        view[:, :, 0] = (acq.plot_corr / npmax(acq.plot_corr) * 255).astype(uint8)
        view[:, :, 1] = (acq.plot_corr / npmax(acq.plot_corr) * 128).astype(uint8)
        view[:, :, 3] = (acq.plot_corr / npmax(acq.plot_corr) * 255).astype(uint8)
        dopp_min, dopp_max = acq.dopp_bins[0], acq.dopp_bins[-1]
        p = figure(title='correlation', x_range=[c1, c2], y_range=[dopp_min, dopp_max], 
                x_axis_label='code phase (samples)', y_axis_label='doppler (Hz)', plot_width=500, plot_height=500)
        p.image_rgba(image=[img], x=[c1], y=[dopp_min], dw=[c2 - c1], dh=[dopp_max - dopp_min])
        return p


def plot_fine_acquisition_results(acq):
    if isinstance(acq, FineAcquirer):
        p = figure(title='fine phase', x_axis_label='time (s)', y_axis_label='phase (rad)', plot_width=500, plot_height=500)
        t = arange(len(acq.phases)) * acq.block_length
        p.line(t, unwrap(angle(acq.phases)))
        return p


class FinePhasePlot(object):
    '''Makes a line plot showing fine acquisition block phases'''

    def __init__(self):
        self.plot = figure(title='fine phase', x_axis_label='time (s)', y_axis_label='phase (rad)', plot_width=500, plot_height=500)
        self.plot.line([0, 1], [0, pi], size=12, alpha=0.7, name='line')

    def update(self, acq):
        renderer = self.plot.select(dict(name='line'))[0]
        angles = unwrap(angle(acq.phases))
        t = arange(len(angles)) * acq.block_length

        print(t[0], t[-1], angles[0], angles[-1])
        renderer.data_source.data['x'] = t
        renderer.data_source.data['y'] = angles
        renderer.data_source._dirty = True
        cursession().store_objects(renderer.data_source)


class CorrelationLinePlot(object):
    '''Makes a line plot showing the max correlation for different Doppler bins'''

    def __init__(self):
        self.plot = figure(title='correlation', x_axis_label='doppler (Hz)', y_axis_label='correlation strength', plot_width=500, plot_height=500)
        self.plot.line([0, 1], [0, 1], name='line')
        self.renderer = self.plot.select(dict(name='line'))[0]

    def update(self, acq):
        self.renderer.data_source.data['x'] = acq.dopp_bins
        self.renderer.data_source.data['y'] = acq.corr_max
        cursession().store_objects(self.renderer.data_source)


class CorrelationGridPlot(object):
    '''Makes a heat map showing correlation grid for CoarseAcquirer object'''

    def __init__(self):
        self.plot = figure(title='correlation', x_axis_label='code phase (samples)', y_axis_label='doppler (Hz)', plot_width=500, plot_height=500)
        fake_corr = asarray([[0, 1], [0, 1]])
        c1, c2 = 0, 1
        dopp_min, dopp_max = 0, 1
        img = self.make_image(fake)
        self.plot.image_rgba(image=[img], x=[c1], y=[c2], dw=[c2 - c1], dh=[dopp_max - dopp_min], name='img')
        self.renderer = self.plot.select(dict(name='img'))[0]

    def make_image(self, corr):
        r, c = corr.shape
        img = ones((r, c), dtype=uint32)
        view = img.view(dtype=uint8).reshape((r, c, 4))
        view[:, :, 0] = (corr / npmax(corr) * 255).astype(uint8)
        view[:, :, 1] = (corr / npmax(corr) * 128).astype(uint8)
        view[:, :, 3] = (corr / npmax(corr) * 255).astype(uint8)
        return img

    def update(self, acq):
        c1, c2 = acq.plot_code_window
        dopp_min, dopp_max = acq.dopp_bins[0], acq.dopp_bins[-1]
        img = self.make_image(acq.plot_corr)
        self.renderer.data_source.data['image'] = img
        self.renderer.data_source.data['x'] = [c1]
        self.renderer.data_source.data['y'] = [c2]
        self.renderer.data_source.data['dw'] = [c2 - c1]
        self.renderer.data_source.data['dh'] = [dopp_max - dopp_min]
        cursession().store_objects(self.renderer.data_source)


