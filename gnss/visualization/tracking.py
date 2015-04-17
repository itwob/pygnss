from numpy import sqrt, ceil, arange
import matplotlib.pyplot as pyplot
from bokeh.plotting import hplot, figure, show, cursession
from gnss.util import RingBuffer

def plot_outputs(store, library='bokeh'):
    keys = store.outputs.keys()
    n = len(keys)
    fig = None
    rows = cols = int(ceil(sqrt(n)))
    if library is 'matplotlib':
        fig = fig if fig else pyplot.figure()
        for i, key in enumerate(keys, 1):
            ax = fig.add_subplot(rows, cols, i)
            ax.plot(store.buffers[key])
            ax.set_title(key)
        return fig
    elif library is 'bokeh':
        plots = []
        for key in store.outputs.keys():
            plot = figure(title=key, plot_width=250, plot_height=250, tools="pan,wheel_zoom,box_zoom,reset,save")
            plot.line(arange(store.outputs[key]['size']), store.buffers[key], size=12, alpha=0.7)
            plots.append(plot)
        plot = hplot(*plots, name="tracking outputs")
        return plot
    return None


class TrackingPlot(object):
    
    def __init__(self, size):
        self.size = size
        
        self.time = RingBuffer(size)
        self.i_corr = RingBuffer(size)
        self.q_corr = RingBuffer(size)
        self.phase_error = RingBuffer(size)
        self.delay_error = RingBuffer(size)
        
        plot = figure(title='i corr', plot_width=250, plot_height=250, tools="pan,wheel_zoom,box_zoom,reset,save")
        plot.line(self.time.get(), self.i_corr.get(), size=12, alpha=0.7, name='i_corr')
        self.i_corr_plot = plot
        
        plot = figure(title='q corr', plot_width=250, plot_height=250, tools="pan,wheel_zoom,box_zoom,reset,save")
        plot.line(self.time.get(), self.q_corr.get(), size=12, alpha=0.7, name='q_corr')
        self.q_corr_plot = plot
        
        plot = figure(title='phase error', plot_width=250, plot_height=250, tools="pan,wheel_zoom,box_zoom,reset,save")
        plot.line(self.time.get(), self.phase_error.get(), size=12, alpha=0.7, name='phase_error')
        self.phase_error_plot = plot
        
        plot = figure(title='delay error', plot_width=250, plot_height=250, tools="pan,wheel_zoom,box_zoom,reset,save")
        plot.line(self.time.get(), self.delay_error.get(), size=12, alpha=0.7, name='delay_error')
        self.delay_error_plot = plot
        
        children = [self.i_corr_plot, self.q_corr_plot, self.phase_error_plot, self.delay_error_plot]
        self.plot = hplot(*children, name="tracking outputs")
        
        self.i_corr_renderer = self.i_corr_plot.select(dict(name='i_corr'))[0]
        self.q_corr_renderer = self.q_corr_plot.select(dict(name='q_corr'))[0]
        self.phase_error_renderer = self.phase_error_plot.select(dict(name='phase_error'))[0]
        self.delay_error_renderer = self.delay_error_plot.select(dict(name='delay_error'))[0]
    
    def show(self):
        show(self.plot)
        self.i_corr_renderer = self.i_corr_plot.select(dict(name='i_corr'))[0]
        self.q_corr_renderer = self.q_corr_plot.select(dict(name='q_corr'))[0]
        self.phase_error_renderer = self.phase_error_plot.select(dict(name='phase_error'))[0]
        self.delay_error_renderer = self.delay_error_plot.select(dict(name='delay_error'))[0]
    
    def push_data(self, time, i_corr, q_corr, phase_error, delay_error):
        self.time.extend(time)
        self.i_corr.extend(i_corr)
        self.q_corr.extend(q_corr)
        self.phase_error.extend(phase_error)
        self.delay_error.extend(delay_error)
    
    def update(self):
        self.i_corr_renderer.data_source.data['x'] = self.time.get()
        self.i_corr_renderer.data_source.data['y'] = self.i_corr.get()
        self.q_corr_renderer.data_source.data['x'] = self.time.get()
        self.q_corr_renderer.data_source.data['y'] = self.q_corr.get()
        self.phase_error_renderer.data_source.data['x'] = self.time.get()
        self.phase_error_renderer.data_source.data['y'] = self.phase_error.get()
        self.delay_error_renderer.data_source.data['x'] = self.time.get()
        self.delay_error_renderer.data_source.data['y'] = self.delay_error.get()
        self.delay_error_renderer.data_source._dirty = True
        cursession().store_objects(self.i_corr_renderer.data_source,
                                   self.q_corr_renderer.data_source,
                                   self.phase_error_renderer.data_source,
                                   self.delay_error_renderer.data_source)
#         cursession().store_objects([i_corr_ds])
