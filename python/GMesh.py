#!/usr/bin/env python

import numpy as np

def is_mesh_uniform(lon,lat):
    """Returns True if the input grid (lon,lat) is uniform and False otherwise"""
    def compare(array):
        eps = np.finfo( array.dtype ).eps # Precision of datatype
        delta = np.abs( array[1:] - array[:-1] ) # Difference along first axis
        error = np.abs( array )
        error = np.maximum( error[1:], error[:-1] ) # Error in difference
        derror = np.abs( delta - delta.flatten()[0] ) # Tolerance to which comparison can be made
        return np.all( derror < ( error + error.flatten()[0] ) )
    assert len(lon.shape) == len(lat.shape), "Arguments lon and lat must have the same rank"
    if len(lon.shape)==2: # 2D arralat
        assert lon.shape == lat.shape, "Arguments lon and lat must have the same shape"
    if len(lon.shape)>2 or len(lat.shape)>2:
        raise Exception("Arguments must be either both be 1D or both be 2D arralat")
    return compare(lat) and compare(lon.T)

def refine_loop(trg_lon_grid,trg_lat_grid, src_lon_grid,src_lat_grid, max_stages=5):
    """This function refines the target grid until all points in the source grid are sampled."""
    """It returns the list of the refined grids."""
    GMesh_list = []
    GMesh_list.append(GMesh(x=trg_lon_grid,y=trg_lat_grid))
    i=0
    hits = GMesh_list[i].source_hits(src_lon_grid,src_lat_grid)
    while(not np.all(hits) and i <= max_stages):
        print("Missed some! Must Refine! Stage ", i+1, "grid shape", GMesh_list[i].x.shape)
        GMesh_list.append(GMesh_list[i].refineby2())
        i=i+1
        hits = GMesh_list[i].source_hits(src_lon_grid,src_lat_grid)

    if(i > max_stages):
        print("Warning: Maximum number of allowed refinements reached without all source points hit.")
    else:
        print("Hit all! Done refining after ",i, " steps!")

    return GMesh_list, hits

class GMesh:
    """Describes 2D meshes for ESMs.

    Meshes have shape=(nj,ni) cells with (nj+1,ni+1) vertices with coordinates (x,y).

    When constructing, either provide 1d or 2d coordinates (x,y), or assume a
    uniform spherical grid with 'shape' cells covering the whole sphere with
    longitudes starting at x0.

    Attributes:

    shape - (nj,ni)
    ni    - number of cells in x-direction (last)
    nj    - number of cells in y-direction (first)
    x     - longitude of mesh (cell corners), shape (nj+1,ni=1)
    y     - latitude of mesh (cell corners), shape (nj+1,ni=1)
    area  - area of cells, shape (nj,ni)
    """

    def __init__(self, shape=None, x=None, y=None, area=None, x0=-180., rfl=0):
        """Constructor for Mesh:
        shape - shape of cell array, (nj,ni)
        ni    - number of cells in x-direction (last index)
        nj    - number of cells in y-direction (first index)
        x     - longitude of mesh (cell corners) (1d or 2d)
        y     - latitude of mesh (cell corners) (1d or 2d)
        area  - area of cells (2d)
        x0    - used when generating a spherical grid in absence of (x,y)
        rfl   - refining level of this mesh
        """
        if (shape is None) and (x is None) and (y is None): raise Exception('Either shape must be specified or both x and y')
        if (x is None) and (y is not None): raise Exception('Either shape must be specified or both x and y')
        if (x is not None) and (y is None): raise Exception('Either shape must be specified or both x and y')
        # Determine shape
        if shape is not None:
            (nj,ni) = shape
        else: # Determine shape from x and y
            if (x is None) or (y is None): raise Exception('Either shape must be specified or both x and y')
            if len(x.shape)==1: ni = x.shape[0]-1
            elif len(x.shape)==2: ni = x.shape[1]-1
            else: raise Exception('x must be 1D or 2D.')
            if len(y.shape)==1 or len(y.shape)==2: nj = y.shape[0]-1
            else: raise Exception('y must be 1D or 2D.')
        self.ni = ni
        self.nj = nj
        self.shape = (nj,ni)
        # Check shape of arrays and construct 2d coordinates
        if x is not None and y is not None:
            if len(x.shape)==1:
                if len(y.shape)>1: raise Exception('x and y must either be both 1d or both 2d')
                if x.shape[0] != ni+1: raise Exception('x has the wrong length')
            if len(y.shape)==1:
                if len(x.shape)>1: raise Exception('x and y must either be both 1d or both 2d')
                if y.shape[0] != nj+1: raise Exception('y has the wrong length')
            if len(x.shape)==2 and len(y.shape)==2:
                if x.shape != y.shape: raise Exception('x and y are 2d and must be the same size')
                if x.shape != (nj+1,ni+1): raise Exception('x has the wrong size')
                self.x = x
                self.y = y
            else:
                self.x, self.y = np.meshgrid(x,y)
        else: # Construct coordinates
            y1d = np.linspace(-90.,90.,nj+1)
            x1d = np.linspace(x0,x0+360.,ni+1)
            self.x, self.y = np.meshgrid(x1d,y1d)
        if area is not None:
            if area.shape != (nj,ni): raise Exception('area has the wrong shape or size')
            self.area = area
        else:
            self.area = None

        self.rfl = rfl #refining level

    def __repr__(self):
        return '<GMesh nj:%i ni:%i shape:(%i,%i)>'%(self.nj,self.ni,self.shape[0],self.shape[1])
    def __getitem__(self, key):
        return getattr(self, key)

    def dump(self):
        print(self)
        print('x.rfl   =',self.rfl)
        print('x.shape =',self.x.shape)
        print('y.shape =',self.y.shape)

    def refineby2(self):
        """Returns new Mesh instance with twice the resolution"""
        x = np.zeros( (2*self.nj+1, 2*self.ni+1) )
        y = np.zeros( (2*self.nj+1, 2*self.ni+1) )
        #area = numpy.zeros( (2*self.nj, 2*self.ni) )
        x[::2,::2] = self.x
        x[::2,1::2] = 0.5 * ( self.x[:,:-1] + self.x[:,1:] )
        x[1::2,::2] = 0.5 * ( self.x[:-1,:] + self.x[1:,:] )
        x[1::2,1::2] = 0.25 * ( ( self.x[:-1,:-1] + self.x[1:,1:] ) + ( self.x[:-1,1:] + self.x[1:,:-1] ) )
        y[::2,::2] = self.y
        y[::2,1::2] = 0.5 * ( self.y[:,:-1] + self.y[:,1:] )
        y[1::2,::2] = 0.5 * ( self.y[:-1,:] + self.y[1:,:] )
        y[1::2,1::2] = 0.25 * ( ( self.y[:-1,:-1] + self.y[1:,1:] ) + ( self.y[:-1,1:] + self.y[1:,:-1] ) )
        return GMesh(x=x, y=y, rfl=self.rfl+1)

    def coarsenby2(self, coarser_mesh):
        """Set the height for lower level Mesh by coarsening"""
        if(self.rfl == 0):
            raise Exception('Coarsest grid, no more coarsening possible!')

        coarser_mesh.height = 0.25*( self.height[::2,::2]
                                   + self.height[1::2,1::2]
                                   + self.height[1:-1:2,0:-1:2]
                                   + self.height[0:-1:2,1:-1:2])

    def find_nn_uniform_source(self,xs,ys):
        """Returns the i&j arrays for the indexes of the nearest neighbor point to each mesh point"""
        #Here we assume that the source mesh {(xs,ys)} is a uniform lat-lon mesh!
        #In this case the index of the closest source point can be easily found by arithmetic.
        assert is_mesh_uniform(xs,ys), 'Grid is not uniform, this method will not work properly'
        delxs = xs[0,1] - xs[0,0]
        delys = ys[1,0] - ys[0,0]
#        nn_i = np.rint((self.x-xs[0,0])/delxs) #Nearest integer (the even one if equidistant)
#        nn_j = np.rint((self.y-ys[0,0])/delys)
        nn_i = np.floor(0.5+(self.x-xs[0,0])/delxs) #Nearest integer (the upper one if equidistant)
        nn_j = np.floor(0.5+(self.y-ys[0,0])/delys)
        #These must be bounded by the extents of the arrays
        upper_i=xs.shape[1]-1
        upper_j=ys.shape[0]-1
        nn_i = np.where(nn_i>upper_i, upper_i, nn_i)
        nn_j = np.where(nn_j>upper_j, upper_j, nn_j)
        nn_i = np.where(nn_i<0, 0, nn_i)
        nn_j = np.where(nn_j<0, 0, nn_j)
        return nn_i.astype(int),nn_j.astype(int)

    def source_hits(self, xs, ys):
        """Returns the number of times each source data point is sampled by this mesh"""
        #This depends on the sampling method
        #Here we assume a Nearest Neighbor sampling.
        #For each GMesh point (x,y):
        #   find the nearest point on the source mesh {(xs,ys)}
        #   increment the number of hits for that source point
        #
        if xs.shape != ys.shape: raise Exception('xs and ys must be the same shape')
        nns_i,nns_j = self.find_nn_uniform_source(xs,ys)
        hits = np.zeros(xs.shape)
        hits[nns_j[:,:],nns_i[:,:]] = 1
#Niki: Deal with the degenerate cases where source points are well outside the target domain
#      and are never going to be hit.
        return hits

    def project_source_data_onto_target_mesh(self,xs,ys,zs):
        """Returns the array on target mesh with values equal to the nearest-neighbor source point data"""
        if xs.shape != ys.shape: raise Exception('xs and ys must be the same shape')
        nns_i,nns_j = self.find_nn_uniform_source(xs,ys)
        self.height = np.zeros(self.x.shape)
        self.height[:,:] = zs[nns_j[:,:],nns_i[:,:]]
        return
