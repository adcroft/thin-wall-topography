from GMesh import GMesh
import numpy

class ThinWalls(GMesh):
    """Container for thin wall topographic data and mesh.

    Additional members:
    zc_mean - mean elevation of cell, shape (nj,ni)
    zu_mean - mean elevation of western edge of cell, shape (nj,ni+1)
    zv_mean - mean elevation of southern edge of cell, shape (nj+1,nj)
    shapeu  - shape of zu_mean, ie. =(nj,ni+1)
    shapev  - shape of zv_mean, ie. =(nj+1,ni)

    Extends the GMesh class.
    """

    def __init__(self, *args, **kwargs):
        """Constructor for ThinWalls."""
        GMesh.__init__(self, *args, **kwargs)
        self.shapeu = (self.nj, self.ni+1)
        self.shapev = (self.nj+1, self.ni)
        self.zc_mean = None
        self.zu_mean = None
        self.zv_mean = None
    def refine(self):
        """Returns new ThinWalls instance with twice the resolution."""
        M = super().refineby2()
        self.shapeu = (self.nj, self.ni+1)
        self.shapev = (self.nj+1, self.ni)
        M.zc_mean = None
        M.zu_mean = None
        M.zv_mean = None
        return M
    def dump(self):
        """Dump Mesh to tty."""
        super().dump()
        print('zc_mean =',self.zc_mean)
        print('zu_mean =',self.zu_mean)
        print('zv_mean =',self.zv_mean)
    def set_cell_mean(self, data):
        """Set elevation of cell center."""
        assert data.shape==self.shape, 'data argument has wrong shape'
        self.zc_mean = data.copy()
    def set_edge_mean(self, datau, datav):
        """Set elevation of cell edges u,v."""
        assert datau.shape==self.shapeu, 'datau argument has wrong shape'
        assert datav.shape==self.shapev, 'datav argument has wrong shape'
        self.zu_mean = datau.copy()
        self.zv_mean = datav.copy()
    def set_edge_mean_to_step(self):
        """Set elevation of cell edges to step topography."""
        self.zu_mean = numpy.zeros(self.shapeu)
        self.zu_mean[:,1:-1] = numpy.maximum( self.zc_mean[:,:-1], self.zc_mean[:,1:] )
        self.zu_mean[:,0] = self.zc_mean[:,0]
        self.zu_mean[:,-1] = self.zc_mean[:,-1]
        #self.zu_mean[:,0] = numpy.maximum( self.zc_mean[:,0], self.zc_mean[:,-1] )
        #self.zu_mean[:,-1] = numpy.maximum( self.zc_mean[:,0], self.zc_mean[:,-1] )
        self.zv_mean = numpy.zeros(self.shapev)
        self.zv_mean[1:-1,:] = numpy.maximum( self.zc_mean[:-1,:], self.zc_mean[1:,:] )
        self.zv_mean[0,:] = self.zc_mean[0,:]
        self.zv_mean[-1,:] = self.zc_mean[-1,:]
        #self.zv_mean[0,:] = numpy.maximum( self.zc_mean[0,:], self.zc_mean[-1,:] )
        #self.zv_mean[-1,:] = numpy.maximum( self.zc_mean[0,:], self.zc_mean[-1,:] )
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
            return pcol_elev( self.zc_mean, self.zu_mean, self.zv_mean )
        else: raise Exception('Unknown "metric"')
