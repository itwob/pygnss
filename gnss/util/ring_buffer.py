from numpy import zeros, arange, asarray

class RingBuffer(object):
    '''Modulo indexing implementation of  Ring Buffer'''
    
    def __init__(self, size, initial_value=None):
        self.buffer = zeros((size,))
        self.size = size
        self.index = 0
        if initial_value:
            self.buffer[:] = initial_value
    
    def __len__(self):
        return self.size
    
    def get(self, size=None):
        if size:
            indices = ((self.index - size + arange(size)) % self.size).astype(int)
            return self.buffer[indices]
        indices = ((self.index + arange(self.size)) % self.size).astype(int)
        return self.buffer[indices]
    
    def extend(self, vals):
        vals = asarray([vals])
        indices = ((self.index + arange(len(vals))) % self.size).astype(int)
        self.buffer[indices] = vals
        self.index = indices[-1] + 1

    def __iadd__(self, vals):
        self.extend(vals)
