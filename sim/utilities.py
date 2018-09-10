import numpy as np

class Plottable(object):
    '''
    A simple base class that allows to store frames in a frame processing
    process so that they can be plotted at the end. This is useful for debbuging.
    '''

    def store(self, name, data):
        ''' Store some data for later plotting '''

        if 'data' not in self.__dict__:
            self.data = {}

        if name in self.data:
            self.data[name].append(data)
        else:
            self.data[name] = [data]


    def plot(self):
        ''' Plot debugging data '''

        if 'data' in self.__dict__:

            import matplotlib.pyplot as plt

            for name, a in self.data.items():
                plt.figure()
                plt.title(name)
                plt.xlabel('Time')
                plt.ylabel('Frequency')
                plt.imshow(np.array(a).T, origin='lower', aspect='auto')
                plt.colorbar()
            plt.show()


class LeakyIntegration(object):
    '''
    A class that simplifies leaky integration of any values
    '''

    def __init__(self, ff, func=None, init=None):
        '''
        Parameters
        ----------
        ff: float
            forgetting factor
        func: func, optional
            A function to apply before the integration
        '''

        self.sum = init
        self.ff = ff

        if func is None:
            self.func = lambda x : x
        else:
            self.func = func

    def update(self, X):
        '''
        Update the integration with new data

        Parameters
        ----------
        X: array_like
            The new data
        '''
        if self.sum is None:
            self.sum = self.func(X)
        else:
            #self.sum = self.ff * self.sum + (1 - self.ff) * self.func(X)
            self.sum = self.ff * self.sum + self.func(X)

    def get(self):
        return self.sum


class ShortTimeAverage(object):
    '''
    A class that simplifies doing short time averaging

    Parameters
    ----------
    L: int
        length of the short time average
    f: func
        function of the input data to average
    '''

    def __init__(self, L, func):

        self.buf = []
        self.sum = None
        self.val = None

        self.f = func
        self.L = L

    def __len__(self):
        return len(self.buf)

    def update(self, X):
        '''
        Update the short time average with new data

        Parameters
        ----------
        X: array_like
            The new data
        '''

        self.buf.insert(0, self.f(X))

        # initialize the output variable
        if self.sum is None:
            self.sum = np.zeros_like(self.buf[0])

        # remove oldest element
        if len(self.buf) > self.L:
            self.sum -= self.buf[-1]
            self.buf.pop()

        # add newest element
        self.sum += self.buf[0]

        # average over the buffer length
        self.val = self.sum / len(self.buf)

    def get(self):
        ''' Return current short time average value '''
        return self.val

    def last(self):
        ''' Return the most recent element in the buffer '''
        return self.buf[0]
