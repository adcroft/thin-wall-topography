#!/usr/bin/env python

import numpy as np
import time

def is_coord_uniform(coord, tol=1.e-5):
    """Returns True if the coordinate "coord" is uniform along the first axis, and False otherwise

    tol is the allowed fractional variation in spacing, i.e. ( variation in delta ) / delta < tol"""
    eps = np.finfo( coord.dtype ).eps # Precision of datatype
    # abscoord = np.abs( coord ) # Magnitude of coordinate values
    # abscoord = np.maximum( abscoord[1:], abscoord[:-1] ) # Largest magnitude of coordinate used at centers
    # roundoff = eps * abscoord # This is the roundoff error in calculating "coord[1:] - coord[:-1]"
    delta = np.abs( coord[1:] - coord[:-1] ) # Spacing along first axis
    roundoff = tol * delta[0] # Assuming delta is approximately uniform, use first value to estimate allowed variance
    derror = np.abs( delta - delta.flatten()[0] ) # delta should be uniform so delta - delta[0] should be zero
    return np.all( derror <= roundoff )

def is_mesh_uniform(lon,lat):
    """Returns True if the input grid (lon,lat) is uniform and False otherwise"""
    assert len(lon.shape) == len(lat.shape), "Arguments lon and lat must have the same rank"
    if len(lon.shape)==2: # 2D array
        assert lon.shape == lat.shape, "Arguments lon and lat must have the same shape"
    assert len(lon.shape)<3 and len(lat.shape)<3, "Arguments must be either both be 1D or both be 2D arralat"
    return is_coord_uniform(lat) and is_coord_uniform(lon.T)

def reg_cell_index_of_x( coord, x, debug=True ):
    assert len(coord.shape) == 1, "coords must be 1D"
    assert is_coord_uniform( coord ), "coord must be uniform"
    k = np.argmin( np.abs( coord[:-1] ) ) # Index of smallest magnitude coordinate (avoid allowing last cell)
    if debug: print('cell with smallest magnitude coord =',k)
    dc = coord[1:] - coord[:-1] # Most accurate estimate of delta (using smallest magnitude coordinate values)
    rdc = 1. / dc[k]
    if debug: print('delta coord =',dc[k],'1/dc =',rdc)
    assert dc[k]>0, "coord must increase with increasing index"
    c_left, c_right = coord[0] - 0.5 * dc[k], coord[-1] + 0.5 * dc[k]
    if debug: print('left,right =',c_left,c_right)
    ind_c = rdc * ( x - c_left )
    if debug: print('non-dim (x-x0)/dx =',ind_c)
    ind_c = ind_c.astype(int)
    periodic_lon = np.abs( ( c_right - c_left ) - 360 ) < 1.e-5 * dc[k] # Detect if this is global longitude in degrees)
    if debug: print('periodic =',periodic_lon)
    if periodic_lon:
        ind_c = np.mod( ind_c, len(coord) )
    else:
        ind_c = np.maximum( 0, np.minimum( len(coord)-1, ind_c ) )
    if debug: print('ind_c=',ind_c)
    return ind_c

class GMesh:
    """Describes 2D meshes for ESMs.

    Meshes have shape=(nj,ni) cells with (nj+1,ni+1) vertices with coordinates (lon,lat).

    When constructing, either provide 1d or 2d coordinates (lon,lat), or assume a
    uniform spherical grid with 'shape' cells covering the whole sphere with
    longitudes starting at lon0.

    Attributes:

    shape - (nj,ni)
    ni    - number of cells in i-direction (last)
    nj    - number of cells in j-direction (first)
    lon   - longitude of mesh (cell corners), shape (nj+1,ni=1)
    lat   - latitude of mesh (cell corners), shape (nj+1,ni=1)
    area  - area of cells, shape (nj,ni)
    """

    def __init__(self, shape=None, lon=None, lat=None, area=None, lon0=-180., from_cell_center=False, rfl=0):
        """Constructor for Mesh:
        shape - shape of cell array, (nj,ni)
        ni    - number of cells in i-direction (last index)
        nj    - number of cells in j-direction (first index)
        lon   - longitude of mesh (cell corners) (1d or 2d)
        lat   - latitude of mesh (cell corners) (1d or 2d)
        area  - area of cells (2d)
        lon0  - used when generating a spherical grid in absence of (lon,lat)
        rfl   - refining level of this mesh
        """
        if (shape is None) and (lon is None) and (lat is None): raise Exception('Either shape must be specified or both lon and lat')
        if (lon is None) and (lat is not None): raise Exception('Either shape must be specified or both lon and lat')
        if (lon is not None) and (lat is None): raise Exception('Either shape must be specified or both lon and lat')
        # Determine shape
        if shape is not None:
            (nj,ni) = shape
        else: # Determine shape from lon and lat
            if (lon is None) or (lat is None): raise Exception('Either shape must be specified or both lon and lat')
            if len(lon.shape)==1: ni = lon.shape[0]-1
            elif len(lon.shape)==2: ni = lon.shape[1]-1
            else: raise Exception('lon must be 1D or 2D.')
            if len(lat.shape)==1 or len(lat.shape)==2: nj = lat.shape[0]-1
            else: raise Exception('lat must be 1D or 2D.')
        if from_cell_center: # Replace cell center coordinates with node coordinates
            ni,nj = ni+1, nj+1
            tmp = np.zeros(ni+1)
            tmp[1:-1] = 0.5 * ( lon[:-1] + lon[1:] )
            tmp[0] = 1.5 * lon[0] - 0.5 * lon[1]
            tmp[-1] = 1.5 * lon[-1] - 0.5 * lon[-2]
            lon = tmp
            tmp = np.zeros(nj+1)
            tmp[1:-1] = 0.5 * ( lat[:-1] + lat[1:] )
            tmp[0] = 1.5 * lat[0] - 0.5 * lat[1]
            tmp[-1] = 1.5 * lat[-1] - 0.5 * lat[-2]
            lat = tmp
        self.ni = ni
        self.nj = nj
        self.shape = (nj,ni)
        # Check shape of arrays and construct 2d coordinates
        if lon is not None and lat is not None:
            if len(lon.shape)==1:
                if len(lat.shape)>1: raise Exception('lon and lat must either be both 1d or both 2d')
                if lon.shape[0] != ni+1: raise Exception('lon has the wrong length')
            if len(lat.shape)==1:
                if len(lon.shape)>1: raise Exception('lon and lat must either be both 1d or both 2d')
                if lat.shape[0] != nj+1: raise Exception('lat has the wrong length')
            if len(lon.shape)==2 and len(lat.shape)==2:
                if lon.shape != lat.shape: raise Exception('lon and lat are 2d and must be the same size')
                if lon.shape != (nj+1,ni+1): raise Exception('lon has the wrong size')
                self.lon = lon
                self.lat = lat
            else:
                self.lon, self.lat = np.meshgrid(lon,lat)
        else: # Construct coordinates
            lon1d = np.linspace(-90.,90.,nj+1)
            lat1d = np.linspace(lon0,lon0+360.,ni+1)
            self.lon, self.lat = np.meshgrid(lon1d,lat1d)
        if area is not None:
            if area.shape != (nj,ni): raise Exception('area has the wrong shape or size')
            self.area = area
        else:
            self.area = None

        self.rfl = rfl #refining level

    # def read_center_coords(self, lont, latt):
    #     if lont.shape != (self.nj,self.ni): raise Exception('Cell center lon has the wrong size')
    #     if latt.shape != (self.nj,self.ni): raise Exception('Cell center lat has the wrong size')
    #     self.lont, self.latt = lont, latt

    def __copy__(self):
        return GMesh(shape = self.shape, lon=self.lon, lat=self.lat, area=self.area)
    def copy(self):
        """Returns new instance with copied values"""
        return self.__copy__()
    def __repr__(self):
        return '<%s nj:%i ni:%i shape:(%i,%i)>'%(self.__class__.__name__,self.nj,self.ni,self.shape[0],self.shape[1])
    def __getitem__(self, key):
        return getattr(self, key)
    def transpose(self):
        """Transpose data swapping i-j indexes"""
        self.ni, self.nj = self.nj, self.ni
        self.shape = (self.nj, self.ni)
        self.lat, self.lon = self.lon.T, self.lat.T
        if self.area is not None: self.area = self.area.T
    def dump(self):
        """Dump Mesh to tty."""
        print(self)
        print('lon = ',self.lon)
        print('lat = ',self.lat)

    def plot(self, axis, subsample=1, linecolor='k', **kwargs):
        for i in range(0,self.ni+1,subsample):
            axis.plot(self.lon[:,i], self.lat[:,i], linecolor, **kwargs)
        for j in range(0,self.nj+1,subsample):
            axis.plot(self.lon[j,:], self.lat[j,:], linecolor, **kwargs)

    def pcolormesh(self, axis, data, **kwargs):
        return axis.pcolormesh( self.lon, self.lat, data, **kwargs)

    def __lonlat_to_XYZ(lon, lat):
        """Private method. Returns 3d coordinates (X,Y,Z) of spherical coordiantes (lon,lat)."""
        deg2rad = np.pi/180.
        lonr,latr = deg2rad*lon, deg2rad*lat
        return np.cos( latr )*np.cos( lonr ), np.cos( latr )*np.sin( lonr ), np.sin( latr )

    def __XYZ_to_lonlat(X, Y, Z):
        """Private method. Returns spherical coordinates (lon,lat) of 3d coordinates (X,Y,Z)."""
        rad2deg = 180./np.pi
        lat = np.arcsin( Z ) * rad2deg # -90 .. 90
        # Normalize X,Y to unit circle
        sub_roundoff = 2./np.finfo(X[0,0]).max
        R = 1. / ( np.sqrt(X*X + Y*Y) + sub_roundoff )
        lon = np.arccos( R*X ) * rad2deg # 0 .. 180
        lon = np.where( Y>=0, lon, -lon ) # Handle -180 .. 0
        return lon,lat

    def interp_center_coords(self, work_in_3d=True):
        """Returns interpolated center coordinates from nodes"""

        def mean4(A):
            """Retruns a refined variable a with shape (2*nj+1,2*ni+1) by linearly interpolation A with shape (nj+1,ni+1)."""
            return 0.25 * ( ( A[:-1,:-1] + A[1:,1:] ) + ( A[1:,:-1] + A[:-1,1:] ) ) # Mid-point of cell on original mesh

        if work_in_3d:
            # Calculate 3d coordinates of nodes (X,Y,Z), Z points along pole, Y=0 at lon=0,180, X=0 at lon=+-90
            X,Y,Z = GMesh.__lonlat_to_XYZ(self.lon, self.lat)

            # Refine mesh in 3d and project onto sphere
            X,Y,Z = mean4(X), mean4(Y), mean4(Z)
            R = 1. / np.sqrt((X*X + Y*Y) + Z*Z)
            X,Y,Z = R*X, R*Y, R*Z

            # Normalize X,Y to unit circle
            #sub_roundoff = 2./np.finfo(X[0,0]).max
            #R = 1. / ( np.sqrt(X*X + Y*Y) + sub_roundoff )
            #X = R * X
            #Y = R * Y

            # Convert from 3d to spherical coordinates
            lon,lat = GMesh.__XYZ_to_lonlat(X, Y, Z)
        else:
            lon,lat = mean4(self.lon), mean4(self.lat)
        return lon, lat

    def coarsest_resolution(self, mask_idx=[]):
        """Returns the coarsest resolution at each grid"""
        def mdist(x1, x2):
            """Returns positive distance modulo 360."""
            return np.minimum(np.mod(x1 - x2, 360.0), np.mod(x2 - x1, 360.0))
        l, p = self.lon, self.lat
        del_lam = np.maximum(np.maximum(np.maximum(mdist(l[:-1,:-1], l[:-1,1:]), mdist(l[1:,:-1], l[1:,1:])),
                                        np.maximum(mdist(l[:-1,:-1], l[1:,:-1]), mdist(l[1:,1:], l[:-1,1:]))),
                             np.maximum(mdist(l[:-1,:-1], l[1:,1:]), mdist(l[1:,:-1], l[:-1,1:])))
        del_phi = np.maximum(np.maximum(np.maximum(np.abs(np.diff(p, axis=0))[:,1:], np.abs((np.diff(p, axis=0))[:,:-1])),
                                        np.maximum(np.abs(np.diff(p, axis=1))[1:,:], np.abs((np.diff(p, axis=1))[:-1,:]))),
                             np.maximum(np.abs(p[:-1,:-1]-p[1:,1:]), np.abs(p[1:,:-1]-p[:-1,1:])))
        if len(mask_idx)>0:
            for Js, Je, Is, Ie in mask_idx:
                jst, jed, ist, ied = Js*(2**self.rfl), Je*(2**self.rfl), Is*(2**self.rfl), Ie*(2**self.rfl)
                del_lam[jst:jed, ist:ied], del_phi[jst:jed, ist:ied] = 0.0, 0.0
        return del_lam, del_phi

    def refineby2(self, work_in_3d=True):
        """Returns new Mesh instance with twice the resolution"""

        def local_refine(A):
            """Retruns a refined variable a with shape (2*nj+1,2*ni+1) by linearly interpolation A with shape (nj+1,ni+1)."""
            nj,ni = A.shape
            a = np.zeros( (2*nj-1,2*ni-1) )
            a[::2,::2] = A[:,:] # Shared nodes
            a[::2,1::2] = 0.5 * ( A[:,:-1] + A[:,1:] ) # Mid-point along i-direction on original mesh
            a[1::2,::2] = 0.5 * ( A[:-1,:] + A[1:,:] ) # Mid-point along j-direction on original mesh
            a[1::2,1::2] = 0.25 * ( ( A[:-1,:-1] + A[1:,1:] ) + ( A[1:,:-1] + A[:-1,1:] ) ) # Mid-point of cell on original mesh
            return a

        if work_in_3d:
            # Calculate 3d coordinates of nodes (X,Y,Z), Z points along pole, Y=0 at lon=0,180, X=0 at lon=+-90
            X,Y,Z = GMesh.__lonlat_to_XYZ(self.lon, self.lat)

            # Refine mesh in 3d and project onto sphere
            X,Y,Z = local_refine(X), local_refine(Y), local_refine(Z)
            R = 1. / np.sqrt((X*X + Y*Y) + Z*Z)
            X,Y,Z = R*X, R*Y, R*Z

            # Normalize X,Y to unit circle
            #sub_roundoff = 2./np.finfo(X[0,0]).max
            #R = 1. / ( np.sqrt(X*X + Y*Y) + sub_roundoff )
            #X = R * X
            #Y = R * Y

            # Convert from 3d to spherical coordinates
            lon,lat = GMesh.__XYZ_to_lonlat(X, Y, Z)

        else:
            lon,lat = local_refine(self.lon), local_refine(self.lat)

        return GMesh(lon=lon, lat=lat, rfl=self.rfl+1)

    def rotate(self, y_rot=0, z_rot=0):
        """Sequentially apply a rotation about the Y-axis and then the Z-axis."""
        deg2rad = np.pi/180.
        # Calculate 3d coordinates of nodes (X,Y,Z), Z points along pole, Y=0 at lon=0,180, X=0 at lon=+-90
        X,Y,Z = GMesh.__lonlat_to_XYZ(self.lon, self.lat)
        # Rotate anti-clockwise about Y-axis
        C,S = np.cos( deg2rad*y_rot ), np.sin( deg2rad*y_rot )
        X,Z = C*X + S*Z, -S*X + C*Z
        # Rotate anti-clockwise about Y-axis
        C,S = np.cos( deg2rad*z_rot ), np.sin( deg2rad*z_rot )
        X,Y = C*X - S*Y, S*X + C*Y

        # Convert from 3d to spherical coordinates
        self.lon,self.lat = GMesh.__XYZ_to_lonlat(X, Y, Z)

        return self

    def coarsenby2(self, coarser_mesh, debug=False, timers=False):
        """Set the height for lower level Mesh by coarsening"""
        if(self.rfl == 0):
            raise Exception('Coarsest grid, no more coarsening possible!')

        if timers: gtic = GMesh._toc(None, "")
        coarser_mesh.height = 0.25 * ( ( self.height[:-1:2,:-1:2] + self.height[1::2,1::2] )
                                     + ( self.height[1::2,:-1:2] + self.height[:-1:2,1::2] ) )
        if timers: gtic = GMesh._toc(gtic, "Whole process")

    def find_nn_uniform_source(self, lon, lat, use_center=True, debug=False):
        """Returns the i,j arrays for the indexes of the nearest neighbor centers at (lon,lat) to the self nodes
        The option use_center=True is default so that lon,lat are cell-center coordinates."""

        assert is_mesh_uniform(lon,lat), 'Grid (lon,lat) is not uniform, this method will not work properly'
        assert len(lon.shape) == len(lat.shape), "lon and lat must both be either 1D or 2D"
        if len(lon.shape)==2:
            # Convert to 1D arrays
            lon,lat = lon[0,:],lat[:,0]
        sni,snj =lon.shape[0],lat.shape[0] # Shape of source cell centered data
        if use_center:
            # Searching for source cells that the self centers fall into
            lon_tgt, lat_tgt = self.interp_center_coords(work_in_3d=True)
        else:
            # Searching for source cells that the self nodes fall into
            lon_tgt, lat_tgt = self.lon, self.lat
        nn_i = reg_cell_index_of_x( lon, lon_tgt, debug=debug )
        nn_j = reg_cell_index_of_x( lat, lat_tgt, debug=debug )
        if debug:
            print('Self lon =',lon[0],'...',lon[-1])
            print('Self lat =',lat[0],'...',lat[-1])
            print('Target lon =',lon_tgt)
            print('Target lat =',lat_tgt)
            print('Source lon =',lon[nn_i])
            print('Source lat =',lat[nn_j])
            print('NN i =',nn_i)
            print('NN j =',nn_j)
        assert nn_j.min()>=0, 'Negative j index calculated! j='+str(nn_j.min())
        assert nn_j.max()<snj, 'Out of bounds j index calculated! j='+str(nn_j.max())
        assert nn_i.min()>=0, 'Negative i index calculated! i='+str(nn_i.min())
        assert nn_i.max()<sni, 'Out of bounds i index calculated! i='+str(nn_i.max())
        return nn_i,nn_j

    def source_hits(self, xs, ys, use_center=True, singularity_radius=0.25):
        """Returns an mask array of 1's if a cell with center (xs,ys) is intercepted by a node
           on the mesh, 0 if no node falls in a cell"""
        # Indexes of nearest xs,ys to each node on the mesh
        i,j = self.find_nn_uniform_source(xs,ys,use_center=use_center)
        sni,snj = xs.shape[0],ys.shape[0] # Shape of source
        hits = np.zeros((snj,sni))
        if singularity_radius>0: hits[np.abs(ys)>90-singularity_radius] = 1
        hits[j,i] = 1
        return hits

    def _toc(tic, label):
        if tic is not None:
            dt = time.time_ns() - tic
            print( '{:>10}ms : {}'.format( dt // 1000000, label) )
        return time.time_ns()

    def refine_loop(self, src_lon, src_lat, max_stages=32, max_mb=32000, fixed_refine_level=0, work_in_3d=True,
                    use_center=True, resolution_limit=True, mask_res=[], singularity_radius=0.25, verbose=True, timers=False):
        """Repeatedly refines the mesh until all cells in the source grid are intercepted by mesh nodes.
           Returns a list of the refined meshes starting with parent mesh."""
        if timers: gtic = GMesh._toc(None, "")
        GMesh_list, this = [self], self
        converged = False
        if fixed_refine_level<1:
            hits = this.source_hits(src_lon, src_lat, use_center=use_center, singularity_radius=singularity_radius)
            nhits, prev_hits = hits.sum().astype(int), 0
            converged = converged or np.all(hits) or (nhits==prev_hits)
        mb = 2*8*this.shape[0]*this.shape[1]/1024/1024
        if resolution_limit:
            sni,snj = src_lon.shape[0],src_lat.shape[0]
            dellon_s, dellat_s = (src_lon[-1]-src_lon[0])/(sni-1), (src_lat[-1]-src_lat[0])/(snj-1)
            del_lam, del_phi = this.coarsest_resolution(mask_idx=mask_res)
            dellon_t, dellat_t = del_lam.max(), del_phi.max()
            converged = converged or ( (dellon_t<=dellon_s) and (dellat_t<=dellat_s) )
        if timers: tic = GMesh._toc(gtic, "Set up")
        if verbose:
            print('Refine level', this.rfl, this, end=" ")
            if fixed_refine_level<1:
                print('Hit', nhits, 'out of', hits.size, 'cells', end=" ")
            if resolution_limit:
                print('dx~1/{} dy~1/{}'.format(int(1/dellon_t), int(1/dellat_t)), end=" ")
            print('(%.4f'%mb,'Mb)')
        # Conditions to refine
        # 1) Not all cells are intercepted
        # 2) A refinement intercepted more cells
        # 3) [if resolution_limit] Coarsest resolution in each direction is coarser than source.
        #    This avoids the excessive refinement which is essentially extrapolation.
        while ( (not converged) \
               and (len(GMesh_list)<max_stages) \
               and (4*mb<max_mb) \
               and (fixed_refine_level<1) \
              ) or (this.rfl < fixed_refine_level):
            if timers: tic = GMesh._toc(None, "")
            this = this.refineby2(work_in_3d=work_in_3d)
            if timers: stic = GMesh._toc(tic, "refine by 2")
            # Find nearest neighbor indices into source
            if fixed_refine_level<1:
                hits = this.source_hits(src_lon, src_lat, singularity_radius=singularity_radius)
                if timers: stic = GMesh._toc(stic, "calculate hits on topo grid")
                nhits, prev_hits = hits.sum().astype(int), nhits
                converged = converged or np.all(hits) or (nhits==prev_hits)
            mb = 2*8*this.shape[0]*this.shape[1]/1024/1024
            if resolution_limit:
                del_lam, del_phi = this.coarsest_resolution(mask_idx=mask_res)
                dellon_t, dellat_t = del_lam.max(), del_phi.max()
                converged = converged or ( (dellon_t<=dellon_s) and (dellat_t<=dellat_s) )
                if timers: stic = GMesh._toc(stic, "calculate resolution stopping criteria")
            GMesh_list.append( this )
            if timers: stic = GMesh._toc(stic, "extending list")
            if timers: tic = GMesh._toc(tic, "Total for loop")
            if verbose:
                print('Refine level', this.rfl, this, end=" ")
                if fixed_refine_level<1:
                    print('Hit', nhits, 'out of', hits.size, 'cells', end=" ")
                if resolution_limit:
                    print('dx~1/{} dy~1/{}'.format(int(1/dellon_t), int(1/dellat_t)), end=" ")
                print('(%.4f'%mb,'Mb)')

        if not converged and fixed_refine_level<1:
            print("Warning: Maximum number of allowed refinements reached without all source cells hit.")
        if timers: tic = GMesh._toc(gtic, "Total for whole process")

        return GMesh_list

    def project_source_data_onto_target_mesh(self, xs, ys, zs, use_center=True, timers=False):
        """Returns the array on target mesh with values equal to the nearest-neighbor source point data"""
        assert len(xs.shape) == len(ys.shape), "xs,ys must both be either 1D or 2D"
        # if xs.shape != ys.shape: raise Exception('xs and ys must be the same shape')
        if timers: gtic = GMesh._toc(None, "")
        nns_i,nns_j = self.find_nn_uniform_source(xs,ys,use_center=use_center)
        if timers: tic = GMesh._toc(gtic, "Calculate interpolation indexes")
        if use_center:
            self.height = np.zeros((self.nj,self.ni))
        else:
            self.height = np.zeros((self.nj+1,self.ni+1))
        if timers: tic = GMesh._toc(tic, "Allocate memory")
        self.height[:,:] = zs[nns_j[:,:],nns_i[:,:]]
        if timers: tic = GMesh._toc(tic, "indirect indexing")
        if timers: tic = GMesh._toc(gtic, "Whole process")
        return

    def find_source_spanning_slices(self, lon, lat, debug=False):
        """Finds the slice ranges of the cell data, described by lon and lat cell-center coordinates, that span the mesh."""

        assert is_mesh_uniform(lon,lat), 'Grid (lon,lat) is not uniform, this method will not work properly'
        assert len(lon.shape) == len(lat.shape), "lon and lat must both be either 1D or 2D"
        if len(lon.shape)==2:
            # Convert to 1D arrays
            lon,lat = lon[0,:],lat[:,0]

        nn_i, nn_j = self.find_nn_uniform_source(lon, lat, use_center=False, debug=debug)
        return slice( nn_i.min(), nn_i.max()+1 ), slice( nn_j.min(), nn_j.max() )

class RegularCoord:
    """Container for uniformly spaced global cell center coordinate parameters

    For use with uniformly gridded data that has cell center global coordinates"""
    def __init__( self, n, origin, periodic, degppi=180 ):
        """Create a RegularCoord
        n         is number of cells;
        origin    is the coordinate on the left edge (not first);
        periodic  distinguishes between longitude and latitude
        """
        self.n = n # Global parameter
        self.periodic = periodic # Global parameter
        if periodic: self.delta, self.rdelta = ( 2 * degppi ) / n, n / ( 2 * degppi )  # Global parameter
        else: self.delta, self.rdelta = degppi / n, n / degppi # Global parameter
        self.origin = origin # Global parameter
        self.offset = np.floor( self.rdelta * self.origin ).astype(int) # Global parameter
        self.rem = np.mod( self.rdelta * self.origin, 1 ) # Global parameter ( needed for odd n)
        self.start = 0 # Special for each subset
        self.stop = self.n # Special for each subset
    def __repr__( self ):
        return '<RegularCoord n={}, dx={}, rdx={}, x0={}, io={}, rem={}, is-ie={}-{}, periodic={}>'.format( \
            self.n, self.delta, self.rdelta, self.origin, self.offset, self.rem, self.start, self.stop, self.periodic)
    def subset( self, slc ):
        """Subset a RegularCoord with slice "slc" """
        Is, Ie = 0, self.n
        if slc.start is not None: Is = slc.start
        if slc.stop is not None: Ie = slc.stop
        S = RegularCoord( self.n, self.origin, self.periodic ) # This creates a copy of "self"
        S.start, S.stop = Is, Ie
        return S
    def indices( self, x ):
        """Return indices of cells that contain x

        If RegularCoord is non-periodic (i.e. latitude), out of range values of "x" will be clipped to -90..90 .
        If regularCoord is periodic, any value of x will be globally wrapped.
        If RegularCoord is a subset, then "x" will be clipped to the bounds of the subset (after periodic wrapping).
        """
        ind = np.floor( self.rdelta * np.array(x) - self.rem ).astype(int) - self.offset
        if self.periodic:
            ind = np.mod( ind, self.n )
        else:
            ind = np.maximum( 0, np.minimum( self.n - 1, ind ) )
        ind = np.maximum( self.start, np.minimum( self.stop - 1, ind ) ) - self.start
        return ind

class UniformEDS:
    """Container for a uniform elevation dataset"""
    def __init__( self, lon, lat, elevation=None ):
        """(lon,lat) are cell centers and 1D with combined shape equalt that of elevation."""
        assert len(lon.shape) == 1, "Longitude must be 1D"
        assert len(lat.shape) == 1, "Latitude must be 1D"
        self.ni = len(lon)
        self.nj = len(lat)
        # Store coordinates for posterity
        self.lonh, self.lath = lon, lat

        if elevation is None: # We allow the creation of a "blank" UniformEDS when subsetting
            self.lon_coord, self.lat_coord = None, None
            self.lonq, self.latq = None, None
            self.data = np.zeros((0))
        else: # This is the real constructor
            assert len(lon) == elevation.shape[1], "Inconsistent longitude shape"
            assert len(lat) == elevation.shape[0], "Inconsistent latitude shape"
            dlon, dlat = 360. / self.ni, 180 / self.nj
            assert np.abs( lon[-1] - lon[0] - 360 + dlon ) < 1.e-5 * dlon, "longitude does not appear to be global"
            assert np.abs( lat[-1] - lat[0] - 180 + dlat ) < 1.e-5 * dlat, "latitude does not appear to be global"
            lon0 = np.floor( lon[0] - 0.5 * dlon + 0.5 ) # Calculating the phase this way restricts ourselves to data starting on integer values
            print('LON0 =',lon0)
            assert np.abs( lon[0] - 0.5 * dlon - lon0 ) < 1.e-9 * dlon, "edge of longitude is not a round number"
            assert np.abs( lat[0] - 0.5 * dlat + 90 ) < 1.e-9 * dlat, "edge of latitude is not 90"
            self.lon_coord = RegularCoord( self.ni, lon0, True)
            self.lat_coord = RegularCoord( self.nj, -90, False)
            # Calculate node coordinates for convenient plotting
            self.lonq = lon0 + dlon * ( np.arange( self.ni + 1 ) )
            self.latq = dlat * ( np.arange( self.nj + 1 ) - 0.5 * self.nj )
            self.data = elevation
    def __repr__( self ):
        return '<UniformEDS {} x {}\nlon = {}\n{}\n{}\nlat = {}\n{}\n{}\ndata = {}>'.format( \
            self.ni, self.nj, self.lon_coord, self.lonh, self.lonq, self.lat_coord, self.lath, self.latq, self.data.shape )
    def subset( self, islice, jslice ):
        """Subset a UniformEDS as [jslice,islice]"""
        S = UniformEDS( self.lonh[islice], self.lath[jslice] )
        S.lon_coord = self.lon_coord.subset(islice)
        S.lat_coord = self.lat_coord.subset(jslice)
        S.lonq = self.lonq[ slice( islice.start, islice.stop + 1 ) ]
        S.latq = self.latq[ slice( jslice.start, jslice.stop + 1 ) ]
        S.data = self.data[jslice, islice]
        return S
    def indices( self, lon, lat ):
        """Return the i,j indices of cells in which (lon,lat) fall"""
        return self.lon_coord.indices( lon ), self.lat_coord.indices( lat )
    def bb_slices( self, lon, lat ):
        """Returns the slices defining the bounding box of data hit by (lon,lat)"""
        si, sj = self.indices( lon, lat )
        return slice( si.min(), si.max() +1 ), slice( sj.min(), sj.max() + 1 )