from GMesh import GMesh
import numpy

class Stats:
    """Container for statistics fields

    shape - shape of these arrays
    min   - minimum value
    max   - maximum value
    mean  - mean value
    """
    def __init__(self, shape, mean=None, min=None, max=None):
        self.shape = shape
        self.min = numpy.zeros(shape)
        self.max = numpy.zeros(shape)
        self.mean = numpy.zeros(shape)
        if mean is not None: self.set_equal(mean)
        if min is not None: self.set_equal(min)
        if max is not None: self.set_equal(max)
    def __repr__(self):
        return '<Stats shape:(%i,%i)>'%(self.shape[0], self.shape[1])
    def __copy__(self):
        return Stats(self.shape, mean=self.mean, min=self.min, max=self.max)
    def copy(self):
        """Returns new instance with copied values"""
        return self.__copy__()
    def dump(self):
        print('min:')
        print(self.min)
        print('mean:')
        print(self.mean)
        print('max:')
        print(self.max)
    def set_equal(self, values):
        assert values.shape == self.shape, 'Data has the wrong shape!'
        self.mean = values.copy()
        self.min = values.copy()
        self.max = values.copy()
    def set(self, min, max, mean):
        assert min.shape == self.shape, 'Min data has the wrong shape!'
        assert max.shape == self.shape, 'Max data has the wrong shape!'
        assert mean.shape == self.shape, 'Mean data has the wrong shape!'
        self.mean = mean.copy()
        self.min = min.copy()
        self.max = max.copy()
    def mean4(self):
        """Return 2d/4-point mean"""
        return 0.25*( (self.mean[::2,::2]+self.mean[1::2,1::2]) + (self.mean[::2,1::2]+self.mean[1::2,::2]) )
    def min4(self):
        """Return 2d/4-point minimum"""
        return numpy.minimum( numpy.minimum( self.min[::2,::2], self.min[1::2,1::2]),
                              numpy.minimum( self.min[::2,1::2], self.min[1::2,::2]) )
    def max4(self):
        """Return 2d/4-point maximum"""
        return numpy.maximum( numpy.maximum( self.max[::2,::2], self.max[1::2,1::2]),
                              numpy.maximum( self.max[::2,1::2], self.max[1::2,::2]) )
    def mean2u(self):
        """Return 2d/2-point mean on u-edges"""
        return 0.5*( self.mean[::2,::2] + self.mean[1::2,::2] )
    def min2u(self):
        """Return 2d/2-point minimum on u-edges"""
        return numpy.minimum( self.min[::2,::2], self.min[1::2,::2] )
    def max2u(self):
        """Return 2d/2-point maximum on u-edges"""
        return numpy.maximum( self.max[::2,::2], self.max[1::2,::2] )
    def mean2v(self):
        """Return 2d/2-point mean on v-edges"""
        return 0.5*( self.mean[::2,::2] + self.mean[::2,1::2] )
    def min2v(self):
        """Return 2d/2-point minimum on v-edges"""
        return numpy.minimum( self.min[::2,::2], self.min[::2,1::2] )
    def max2v(self):
        """Return 2d/2-point maximum on v-edges"""
        return numpy.maximum( self.max[::2,::2], self.max[::2,1::2] )

class ThinWalls(GMesh):
    """Container for thin wall topographic data and mesh.

    Additional members:
    c_simple - elevation statistics of cell, shape (nj,ni)
    u_simple - elevation statistics of western edge of cell, shape (nj,ni+1)
    v_simple - elevation statistics of southern edge of cell, shape (nj+1,nj)
    shapeu  - shape of zu_simple_mean, ie. =(nj,ni+1)
    shapev  - shape of zv_simple_mean, ie. =(nj+1,ni)

    Extends the GMesh class.
    """

    def __init__(self, *args, **kwargs):
        """Constructor for ThinWalls."""
        GMesh.__init__(self, *args, **kwargs)
        self.shapeu = (self.nj, self.ni+1)
        self.shapev = (self.nj+1, self.ni)
        self.c_simple = Stats(self.shape)
        self.u_simple = Stats(self.shapeu)
        self.v_simple = Stats(self.shapev)
    def refine(self):
        """Returns new ThinWalls instance with twice the resolution."""
        M = super().refineby2()
        return ThinWalls(lon=M.lon, lat=M.lat)
    def dump(self):
        """Dump Mesh to tty."""
        super().dump()
        self.c_simple.dump()
        self.u_simple.dump()
        self.v_simple.dump()
        self.c_effective.dump()
        self.u_effective.dump()
        self.v_effective.dump()
    def set_cell_mean(self, data):
        """Set elevation of cell center."""
        assert data.shape==self.shape, 'data argument has wrong shape'
        self.c_simple.set_equal(data)
    def set_edge_mean(self, datau, datav):
        """Set elevation of cell edges u,v."""
        assert datau.shape==self.shapeu, 'datau argument has wrong shape'
        assert datav.shape==self.shapev, 'datav argument has wrong shape'
        self.u_simple.set_equal(datau)
        self.u_simple.set_equal(datau)
        self.v_simple.set_equal(datav)
    def init_effective_values(self):
        """Initialize effective values by setting equal to simple values."""
        self.c_effective = self.c_simple.copy()
        self.u_effective = self.u_simple.copy()
        self.v_effective = self.v_simple.copy()
    def set_edge_to_step(self):
        """Set elevation of cell edges to step topography."""
        tmp = numpy.zeros(self.shapeu)
        tmp[:,1:-1] = numpy.maximum( self.c_simple.mean[:,:-1], self.c_simple.mean[:,1:] )
        tmp[:,0] = self.c_simple.mean[:,0]
        tmp[:,-1] = self.c_simple.mean[:,-1]
        self.u_simple.set_equal( tmp )
        tmp = numpy.zeros(self.shapev)
        tmp[1:-1,:] = numpy.maximum( self.c_simple.mean[:-1,:], self.c_simple.mean[1:,:] )
        tmp[0,:] = self.c_simple.mean[0,:]
        tmp[-1,:] = self.c_simple.mean[-1,:]
        self.v_simple.set_equal( tmp )
    def coarsen(self):
        M = ThinWalls(lon=self.lon[::2,::2],lat=self.lat[::2,::2])
        M.c_simple.mean = self.c_simple.mean4()
        M.c_simple.min = self.c_simple.min4()
        M.c_simple.max = self.c_simple.max4()
        M.u_simple.mean =self.u_simple.mean2u()
        M.u_simple.min =self.u_simple.min2u()
        M.u_simple.max =self.u_simple.max2u()
        M.v_simple.mean = self.v_simple.mean2v()
        M.v_simple.min = self.v_simple.min2v()
        M.v_simple.max = self.v_simple.max2v()
        return M
    def plot(self, axis, thickness=0.2, metric='mean', measure='simple', *args, **kwargs):
        """Plots ThinWalls data."""
        def copy_coord(xy):
            XY = numpy.zeros( (2*self.nj+2,2*self.ni+2) )
            dr = xy[1:,1:] - xy[:-1,:-1]
            dl = xy[:-1,1:] - xy[1:,:-1]
            XY[::2,::2] = xy
            XY[2::2,2::2] = XY[2::2,2::2] - dr*thickness/2
            XY[1::2,::2] = xy
            XY[1:-1:2,2::2] = XY[1:-1:2,2::2] - dl*thickness/2
            XY[::2,1::2] = xy
            XY[2::2,1:-1:2] = XY[2::2,1:-1:2] + dl*thickness/2
            XY[1::2,1::2] = xy
            XY[1:-1:2,1:-1:2] = XY[1:-1:2,1:-1:2] + dr*thickness/2 
            return XY
        lon = copy_coord(self.lon)
        lat = copy_coord(self.lat)
        def pcol_elev(c,u,v):
            tmp = numpy.ma.zeros( (2*self.nj+1,2*self.ni+1) )
            tmp[::2,::2] = numpy.ma.masked # Mask corner values
            tmp[1::2,1::2] = c
            tmp[1::2,::2] = u
            tmp[::2,1::2] = v
            return axis.pcolormesh(lon, lat, tmp, *args, **kwargs)
        if measure is 'simple':
            c,u,v = self.c_simple, self.u_simple, self.v_simple
        elif measure is 'effective':
            c,u,v = self.c_effective, self.u_effective, self.v_effective
        else: raise Exception('Unknown "measure"')
        if metric is 'mean': return pcol_elev( c.mean, u.mean, v.mean )
        elif metric is 'min': return pcol_elev( c.min, u.min, v.min )
        elif metric is 'max': return pcol_elev( c.max, u.max, v.max )
        else: raise Exception('Unknown "metric"')
    def plot_grid(self, axis, *args, **kwargs):
        """Plots ThinWalls mesh."""
        super().plot(axis, *args, **kwargs)
