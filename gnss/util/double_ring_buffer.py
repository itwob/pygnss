class RingDoubleBuffer(object):
    '''Double buffer implementation of Ring Buffer'''
    
    def __init__(self, size):
        self.buffer = zeros((2 * size,))
        self.size = size
        self.index = 0
    
    def __len__(self):
        return self.size
    
    def __getitem__(self, key):
        '''Adds `self.index` to the slice start and stop. Stop index should not
        go out of bounds of buffer. (Roll backs take place in `__setitem__`'''
        if isinstance(key, int):
            return self.buffer[key]
        key = slice(key.start if key.start else 0, key.stop if key.stop else self.size, key.step if key.step else None)
        return self.buffer[self.index + key.start: self.index + key.stop: key.step]
    
    def __setitem__(self, key, vals):
        '''Set items in buffer directly'''
        if isinstance(key, int):
            self.buffer[self.index] = vals
        else:
            self.buffer[self.index + key.start: self.index + key.stop: key.step] = val
    
    def extend(self, vals):
        if self.index + len(vals) >= len(self.buffer):
            self.buffer = roll(self.buffer, -self.index)  #-(len(self.buffer) - self.index))
            self.index = 0
        self.buffer[self.index:self.index + len(vals)] = vals
        self.index += len(vals)

    def __iadd__(self, vals):
        self.extend(vals)
