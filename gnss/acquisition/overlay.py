

def sync_overlay_code(signal, source, time, chip0, f_dopp, f_overlay, block_size):
    '''Finds the true chip at `time` for a signal that has an overlay
    code with chipping frequency `f_overlay`.'''
    
    for chip in chips:
        correlator.correlate(signal, source, block_size, time, chip, f_dopp, theta)
