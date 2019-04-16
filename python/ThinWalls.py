from GMesh import GMesh
import numpy

class Stats:
    """Container for statistics fields

    shape - shape of these arrays
    min   - minimum value
    max   - maximum value
    mean  - mean value
    """
    def __init__(self, shape, values=None):
        self.shape = shape
        self.min = numpy.zeros(shape)
        self.max = numpy.zeros(shape)
        self.mean = numpy.zeros(shape)
        if values is not None: self.set_equal(values)
    def __repr__(self):
        return '<Stats shape:(%i,%i)>'%(self.shape[0], self.shape[1])
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
    def set_cell_mean(self, data):
        """Set elevation of cell center."""
        assert data.shape==self.shape, 'data argument has wrong shape'
        self.c_simple.set_equal(data)
    def set_edge_mean(self, datau, datav):
        """Set elevation of cell edges u,v."""
        assert datau.shape==self.shapeu, 'datau argument has wrong shape'
        assert datav.shape==self.shapev, 'datav argument has wrong shape'
        self.u_simple.set_equal(datau)
        self.v_simple.set_equal(datav)
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
        M.c_simple.mean = 0.25*( (self.c_simple.mean[::2,::2]+self.c_simple.mean[1::2,1::2])
                         +(self.c_simple.mean[::2,1::2]+self.c_simple.mean[1::2,::2]) )
        M.c_simple.min = numpy.minimum(
                    numpy.minimum( self.c_simple.min[::2,::2], self.c_simple.min[1::2,1::2]),
                    numpy.minimum( self.c_simple.min[::2,1::2], self.c_simple.min[1::2,::2]) )
        M.c_simple.max = numpy.maximum(
                    numpy.maximum( self.c_simple.max[::2,::2], self.c_simple.max[1::2,1::2]),
                    numpy.minimum( self.c_simple.max[::2,1::2], self.c_simple.max[1::2,::2]) )
        M.u_simple.mean = 0.5*( self.u_simple.mean[::2,::2] + self.u_simple.mean[1::2,::2] )
        M.u_simple.min = numpy.minimum( self.u_simple.min[::2,::2], self.u_simple.min[1::2,::2] )
        M.u_simple.max = numpy.maximum( self.u_simple.max[::2,::2], self.u_simple.max[1::2,::2] )
        M.v_simple.mean = 0.5*( self.v_simple.mean[::2,::2] + self.v_simple.mean[::2,1::2] )
        M.v_simple.min = numpy.minimum( self.v_simple.min[::2,::2], self.v_simple.min[::2,1::2] )
        M.v_simple.max = numpy.maximum( self.v_simple.max[::2,::2], self.v_simple.max[::2,1::2] )
        return M
    def plot(self, axis, thickness=0.2, metric='mean', *args, **kwargs):
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
        if metric is 'mean':
            return pcol_elev( self.c_simple.mean, self.u_simple.mean, self.v_simple.mean )
        elif metric is 'min':
            return pcol_elev( self.c_simple.min, self.u_simple.min, self.v_simple.min )
        elif metric is 'max':
            return pcol_elev( self.c_simple.max, self.u_simple.max, self.v_simple.max )
        else: raise Exception('Unknown "metric"')
    def plot_grid(self, axis, *args, **kwargs):
        """Plots ThinWalls mesh."""
        super().plot(axis, *args, **kwargs)
