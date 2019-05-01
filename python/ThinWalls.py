from GMesh import GMesh
import numpy

class Stats:
    """Container for statistics fields

    shape - shape of these arrays
    low   - minimum value
    hgh   - maximum value
    ave   - mean value
    """
    def __init__(self, shape, mean=None, min=None, max=None):
        self.shape = shape
        self.low = numpy.zeros(shape)
        self.hgh = numpy.zeros(shape)
        self.ave = numpy.zeros(shape)
        if mean is not None: self.set_equal(mean)
        if min is not None: self.set_equal(min)
        if max is not None: self.set_equal(max)
    def __repr__(self):
        return '<Stats shape:(%i,%i)>'%(self.shape[0], self.shape[1])
    def __copy__(self):
        return Stats(self.shape, mean=self.ave, min=self.low, max=self.hgh)
    def copy(self):
        """Returns new instance with copied values"""
        return self.__copy__()
    def dump(self):
        print('min:')
        print(self.low)
        print('mean:')
        print(self.ave)
        print('max:')
        print(self.hgh)
    def set_equal(self, values):
        assert values.shape == self.shape, 'Data has the wrong shape!'
        self.ave = values.copy()
        self.low = values.copy()
        self.hgh = values.copy()
    def set(self, min, max, mean):
        assert min.shape == self.shape, 'Min data has the wrong shape!'
        assert max.shape == self.shape, 'Max data has the wrong shape!'
        assert mean.shape == self.shape, 'Mean data has the wrong shape!'
        self.ave = mean.copy()
        self.low = min.copy()
        self.hgh = max.copy()
    def mean4(self):
        """Return 2d/4-point mean"""
        return 0.25*( (self.ave[::2,::2]+self.ave[1::2,1::2]) + (self.ave[::2,1::2]+self.ave[1::2,::2]) )
    def min4(self):
        """Return 2d/4-point minimum"""
        return numpy.minimum( numpy.minimum( self.low[::2,::2], self.low[1::2,1::2]),
                              numpy.minimum( self.low[::2,1::2], self.low[1::2,::2]) )
    def max4(self):
        """Return 2d/4-point maximum"""
        return numpy.maximum( numpy.maximum( self.hgh[::2,::2], self.hgh[1::2,1::2]),
                              numpy.maximum( self.hgh[::2,1::2], self.hgh[1::2,::2]) )
    def mean2u(self):
        """Return 2d/2-point mean on u-edges"""
        return 0.5*( self.ave[::2,::2] + self.ave[1::2,::2] )
    def min2u(self):
        """Return 2d/2-point minimum on u-edges"""
        return numpy.minimum( self.low[::2,::2], self.low[1::2,::2] )
    def max2u(self):
        """Return 2d/2-point maximum on u-edges"""
        return numpy.maximum( self.hgh[::2,::2], self.hgh[1::2,::2] )
    def mean2v(self):
        """Return 2d/2-point mean on v-edges"""
        return 0.5*( self.ave[::2,::2] + self.ave[::2,1::2] )
    def min2v(self):
        """Return 2d/2-point minimum on v-edges"""
        return numpy.minimum( self.low[::2,::2], self.low[::2,1::2] )
    def max2v(self):
        """Return 2d/2-point maximum on v-edges"""
        return numpy.maximum( self.hgh[::2,::2], self.hgh[::2,1::2] )
    def flip(self, axis):
        """Flip the data along the given axis"""
        self.low = numpy.flip(self.low, axis=axis)
        self.ave = numpy.flip(self.ave, axis=axis)
        self.hgh = numpy.flip(self.hgh, axis=axis)
    def transpose(self):
        """Transpose data swapping i-j indexes"""
        self.low = self.low.T
        self.ave = self.ave.T
        self.hgh = self.hgh.T
        #self.shape = self.low.shape

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
        self.c_effective = Stats(self.shape)
        self.u_effective = Stats(self.shapeu)
        self.v_effective = Stats(self.shapev)
    def __copy__(self):
        copy = ThinWalls(shape=self.shape, lon=self.lon, lat=self.lat)
        copy.c_simple = self.c_simple.copy()
        copy.u_simple = self.u_simple.copy()
        copy.v_simple = self.v_simple.copy()
        copy.c_effective = self.c_effective.copy()
        copy.u_effective = self.u_effective.copy()
        copy.v_effective = self.v_effective.copy()
    def copy(self):
        """Returns new instance with copied values"""
        return self.__copy__()
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
        tmp[:,1:-1] = numpy.maximum( self.c_simple.ave[:,:-1], self.c_simple.ave[:,1:] )
        tmp[:,0] = self.c_simple.ave[:,0]
        tmp[:,-1] = self.c_simple.ave[:,-1]
        self.u_simple.set_equal( tmp )
        tmp = numpy.zeros(self.shapev)
        tmp[1:-1,:] = numpy.maximum( self.c_simple.ave[:-1,:], self.c_simple.ave[1:,:] )
        tmp[0,:] = self.c_simple.ave[0,:]
        tmp[-1,:] = self.c_simple.ave[-1,:]
        self.v_simple.set_equal( tmp )
    def push_corners(self, update_interior_mean_max=True):
        """Folds out tallest corners. Acts only on "effective" values.

        A convex corner within a coarse grid cell can be made into a
        concave corner without changing connectivity across the major
        parts of the cell. The cross-corner connection for the minor
        part of the cell is eliminated."""

        self.push_corners_sw(update_interior_mean_max=update_interior_mean_max)
        # Alias
        C, U, V = self.c_effective, self.u_effective, self.v_effective
        # Flip in j direction
        C.flip(axis=0)
        U.flip(axis=0)
        V.flip(axis=0)
        self.push_corners_sw(update_interior_mean_max=update_interior_mean_max) # Push NW
        # Flip in i direction
        C.flip(axis=1)
        U.flip(axis=1)
        V.flip(axis=1)
        self.push_corners_sw(update_interior_mean_max=update_interior_mean_max) # Push NE
        # Flip in j direction
        C.flip(axis=0)
        U.flip(axis=0)
        V.flip(axis=0)
        self.push_corners_sw(update_interior_mean_max=update_interior_mean_max) # Push SE
        # Flip in i direction
        C.flip(axis=1)
        U.flip(axis=1)
        V.flip(axis=1)
    def push_corners_sw(self, update_interior_mean_max=True):
        """Folds out SW corner is it is the highest ridge. Acts only on "effective" values.

        A convex corner within a coarse grid cell can be made into a
        concave corner without changing connectivity across the major
        parts of the cell. The cross-corner connection for the minor
        part of the cell is eliminated."""
        # Alias
        C,U,V = self.c_effective,self.u_effective,self.v_effective
        # Inner SW corner
        crnr_min = numpy.minimum( U.low[::2,1::2], V.low[1::2,::2] )    # Min or "sill" for SW corner
        crnr_mean = 0.5*( U.ave[::2,1::2] + V.ave[1::2,::2] )         # Mean for SW corner
        crnr_max = numpy.maximum( U.hgh[::2,1::2], V.hgh[1::2,::2] )    # Max for SW corner
        # Values for the coarse cell outside of the SW corner
        opp_ridge = numpy.maximum( U.low[1::2,1::2], V.low[1::2,1::2] ) # Ridge for NE corner
        opp_cmean = ( ( C.ave[::2,1::2] + C.ave[1::2,::2] ) + C.ave[1::2,1::2] )/3 # Mean of outer cells
        j,i = numpy.nonzero( crnr_min>opp_ridge )  # Find where the SW corner has the highest sill
        if len(i)>0:
            J,I = 2*j,2*i
            # Replace inner minimum values with ridge value
            # - set inner SW corner sill to peak of the NW ridge to avoid introducing a new deep diagonal
            #   connection across the interior of the coarse cell
            U.low[J,I+1] = opp_ridge[j,i]
            V.low[J+1,I] = opp_ridge[j,i]
            # ????? No replace inner mean and max ???? Not used?
            # Override outer SW edge values with SW corner inner values
            U.low[J,I] = numpy.maximum( U.low[J,I], crnr_min[j,i] )
            V.low[J,I] = numpy.maximum( V.low[J,I], crnr_min[j,i] )
            U.ave[J,I] = numpy.maximum( U.ave[J,I], crnr_mean[j,i] )
            V.ave[J,I] = numpy.maximum( V.ave[J,I], crnr_mean[j,i] )
            U.hgh[J,I] = numpy.maximum( U.hgh[J,I], crnr_max[j,i] )
            V.hgh[J,I] = numpy.maximum( V.hgh[J,I], crnr_max[j,i] )
            # Override SW cell values with outer values from coarse cell
            C.low[J,I] = opp_ridge[j,i] # This will be taller than other minimums but is it lower than opp_cmean ????
            if update_interior_mean_max:
                C.ave[J,I] = numpy.maximum( C.ave[J,I], opp_cmean[j,i] ) # Avoids changing the mean of the remaining coarse cell
                C.hgh[J,I] = numpy.maximum( C.hgh[J,I], opp_ridge[j,i] )   # Will be taller than cell means?
                #opp_ridge = 0.5*( U.ave[1::2,1::2] + V.ave[1::2,1::2] ) # Ridge for NE corner
                U.ave[J,I+1] = opp_ridge[j,i]
                V.ave[J+1,I] = opp_ridge[j,i]
                #opp_ridge = numpy.maximum( U.hgh[1::2,1::2], V.hgh[1::2,1::2] ) # Ridge for NE corner
                U.hgh[J,I+1] = opp_ridge[j,i]
                V.hgh[J+1,I] = opp_ridge[j,i]
    def lower_tallest_buttress(self):
        """Lower tallest barrier to remove buttress"""
        # Alias lowest
        C,U,V = self.c_effective.low,self.u_effective.low,self.v_effective.low
        # Find where the S ridge is higher than other 3
        oppo3 = numpy.maximum( U[1::2,1::2], numpy.maximum( V[1::2,::2], V[1::2,1::2] ) )
        j,i = numpy.nonzero( U[::2,1::2]>oppo3 )
        U[2*j,2*i+1] = oppo3[j,i]
        # Find where the N ridge is higher than other 3
        oppo3 = numpy.maximum( U[::2,1::2], numpy.maximum( V[1::2,::2], V[1::2,1::2] ) )
        j,i = numpy.nonzero( U[1::2,1::2]>oppo3 )
        U[2*j+1,2*i+1] = oppo3[j,i]
        # Find where the W ridge is higher than other 3
        oppo3 = numpy.maximum( V[1::2,1::2], numpy.maximum( U[::2,1::2], U[1::2,1::2] ) )
        j,i = numpy.nonzero( V[1::2,::2]>oppo3 )
        V[2*j+1,2*i] = oppo3[j,i]
        # Find where the E ridge is higher than other 3
        oppo3 = numpy.maximum( V[1::2,::2], numpy.maximum( U[::2,1::2], U[1::2,1::2] ) )
        j,i = numpy.nonzero( V[1::2,1::2]>oppo3 )
        V[2*j+1,2*i+1] = oppo3[j,i]
        # Alias for averages
        C,U,V = self.c_effective.ave,self.u_effective.ave,self.v_effective.ave
        # Find where the S ridge is higher than other 3
        oppo3 = numpy.maximum( U[1::2,1::2], numpy.maximum( V[1::2,::2], V[1::2,1::2] ) )
        j,i = numpy.nonzero( U[::2,1::2]>oppo3 )
        U[2*j,2*i+1] = oppo3[j,i]
        # Find where the N ridge is higher than other 3
        oppo3 = numpy.maximum( U[::2,1::2], numpy.maximum( V[1::2,::2], V[1::2,1::2] ) )
        j,i = numpy.nonzero( U[1::2,1::2]>oppo3 )
        U[2*j+1,2*i+1] = oppo3[j,i]
        # Find where the W ridge is higher than other 3
        oppo3 = numpy.maximum( V[1::2,1::2], numpy.maximum( U[::2,1::2], U[1::2,1::2] ) )
        j,i = numpy.nonzero( V[1::2,::2]>oppo3 )
        V[2*j+1,2*i] = oppo3[j,i]
        # Find where the E ridge is higher than other 3
        oppo3 = numpy.maximum( V[1::2,::2], numpy.maximum( U[::2,1::2], U[1::2,1::2] ) )
        j,i = numpy.nonzero( V[1::2,1::2]>oppo3 )
        V[2*j+1,2*i+1] = oppo3[j,i]
    def fold_out_central_ridges(self):
        """Folded out interior ridges to the sides of the coarse cell"""
        self.fold_out_central_ridge_s()
        # Alias
        C, U, V = self.c_effective, self.u_effective, self.v_effective
        # Flip in j direction so j=S, i=E
        C.flip(axis=0)
        U.flip(axis=0)
        V.flip(axis=0)
        self.fold_out_central_ridge_s()
        # Transpose so j=E, i=S
        C.transpose()
        U.transpose()
        V.transpose()
        self.u_effective,self.v_effective = self.v_effective,self.u_effective
        C, U, V = self.c_effective, self.u_effective, self.v_effective
        self.fold_out_central_ridge_s()
        # Flip in j direction so j=W, i=S
        C.flip(axis=0)
        U.flip(axis=0)
        V.flip(axis=0)
        self.fold_out_central_ridge_s()
        # Undo transformations
        C.transpose()
        U.transpose()
        V.transpose()
        self.u_effective,self.v_effective = self.v_effective,self.u_effective
        C, U, V = self.c_effective, self.u_effective, self.v_effective
        C.flip(axis=0)
        U.flip(axis=0)
        V.flip(axis=0)
        C.flip(axis=1)
        U.flip(axis=1)
        V.flip(axis=1)
    def fold_out_central_ridge_s(self):
        """An interior east-west ridge is folded out to the southern outer edges if it
        is the tallest central ridge and the south is the taller half to expand to."""
        # Alias
        C,U,V = self.c_effective,self.u_effective,self.v_effective
        ew_ridge_low = numpy.minimum( V.low[1::2,::2], V.low[1::2,1::2] )
        #ew_ridge_hgh = numpy.maximum( V.hgh[1::2,::2], V.hgh[1::2,1::2] )
        #ew_ridge_ave = 0.5*( V.low[1::2,::2] + V.low[1::2,1::2] )
        ns_ridge_low_min = numpy.minimum( U.low[::2,1::2], U.low[1::2,1::2] )
        ns_ridge_low_max = numpy.maximum( U.low[::2,1::2], U.low[1::2,1::2] )
        # Coarse cell index j,i
        j,i = numpy.nonzero(
              ( ( ew_ridge_low>ns_ridge_low_min) & (ew_ridge_low>=ns_ridge_low_max ) ) # E-W ridge is the taller ridge
              & (
                  ( U.low[::2,1::2] > U.low[1::2,1::2] ) # Southern buttress is taller than north
                  | (
                      ( U.low[::2,1::2] >= U.low[1::2,1::2] ) # Southern buttress is equal to the north
                      & (
                          ( C.low[::2,::2]+C.low[::2,1::2] > C.low[1::2,::2]+C.low[1::2,1::2] ) | # Southern cells are higher than north on average
                          ( V.low[:-1:2,::2]+V.low[:-1:2,1::2] > V.low[2::2,::2]+V.low[2::2,1::2] ) # Southern edges are higher than north on average
                ) ) ) )
        J,I = 2*j,2*i
        # Outer edges of southern half
        U.low[J,I] = numpy.maximum( U.low[J,I], ew_ridge_low[j,i] )
        V.low[J,I] = numpy.maximum( V.low[J,I], ew_ridge_low[j,i] )
        V.low[J,I+1] = numpy.maximum( V.low[J,I+1], ew_ridge_low[j,i] )
        U.low[J,I+2] = numpy.maximum( U.low[J,I+2], ew_ridge_low[j,i] )
        # Replace E-W ridge
        V.low[J+1,I] = ns_ridge_low_min[j,i]
        V.low[J+1,I+1] = ns_ridge_low_min[j,i]
        # Southern cells
        C.low[J,I] = ns_ridge_low_min[j,i]
        C.low[J,I+1] = ns_ridge_low_min[j,i]
        U.low[J,I+1] = ns_ridge_low_min[j,i]
    def coarsen(self):
        M = ThinWalls(lon=self.lon[::2,::2],lat=self.lat[::2,::2])
        M.c_simple.ave = self.c_simple.mean4()
        M.c_simple.low = self.c_simple.min4()
        M.c_simple.hgh = self.c_simple.max4()
        M.u_simple.ave =self.u_simple.mean2u()
        M.u_simple.low =self.u_simple.min2u()
        M.u_simple.hgh =self.u_simple.max2u()
        M.v_simple.ave = self.v_simple.mean2v()
        M.v_simple.low = self.v_simple.min2v()
        M.v_simple.hgh = self.v_simple.max2v()
        M.c_effective.ave = self.c_effective.mean4()
        M.c_effective.low = self.c_effective.min4()
        M.c_effective.hgh = self.c_effective.max4()
        M.u_effective.ave =self.u_effective.mean2u()
        M.u_effective.low =self.u_effective.min2u()
        M.u_effective.hgh =self.u_effective.max2u()
        M.v_effective.ave = self.v_effective.mean2v()
        M.v_effective.low = self.v_effective.min2v()
        M.v_effective.hgh = self.v_effective.max2v()
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
        if metric is 'mean': return pcol_elev( c.ave, u.ave, v.ave )
        elif metric is 'min': return pcol_elev( c.low, u.low, v.low )
        elif metric is 'max': return pcol_elev( c.hgh, u.hgh, v.hgh )
        else: raise Exception('Unknown "metric"')
    def plot_grid(self, axis, *args, **kwargs):
        """Plots ThinWalls mesh."""
        super().plot(axis, *args, **kwargs)
