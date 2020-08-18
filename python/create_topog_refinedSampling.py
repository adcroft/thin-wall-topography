#!/usr/bin/env python

import sys, getopt
import datetime, os, subprocess
import GMesh
import imp
import netCDF4
import numpy as np

def break_array_to_blocks(a,xb=4,yb=1):
    a_win = []
    if(xb == 4 and yb ==1):
        i1 = a.shape[1]//xb
        i2 = 2*i1
        i3 = 3*i1
        i4 = a.shape[1]
        
        j1 = a.shape[0]//yb
        a_win.append(a[0:j1,0:i1])
        a_win.append(a[0:j1,i1:i2])
        a_win.append(a[0:j1,i2:i3])
        a_win.append(a[0:j1,i3:i4])
        return a_win
    else:    
        raise Exception('This rotuine can only make 2x2 blocks!')
        ##Niki: Implement a better algo and lift this restriction

def undo_break_array_to_blocks(a,xb=4,yb=1):    
    if(xb == 4 and yb ==1):        
        ao = np.append(a[0],a[1],axis=1)
        ao = np.append(ao,a[2],axis=1)
        ao = np.append(ao,a[3],axis=1)
        return ao
    else:    
        raise Exception('This rotuine can only make 2x2 blocks!')
        ##Niki: Implement a better algo and lift this restriction

def write_topog(h,hstd,hmin,hmax,xx,yy,fnam=None,format='NETCDF3_CLASSIC',description=None,history=None,source=None,no_changing_meta=None):
    import netCDF4 as nc

    if fnam is None:
      fnam='topog.nc'
    fout=nc.Dataset(fnam,'w',format=format)

    ny=h.shape[0]; nx=h.shape[1]
    print ('Writing netcdf file ',fnam,' with ny,nx= ',ny,nx)

    ny=fout.createDimension('ny',ny)
    nx=fout.createDimension('nx',nx)
    string=fout.createDimension('string',255)    
    tile=fout.createVariable('tile','S1',('string'))
    height=fout.createVariable('height','f8',('ny','nx'))
    height.units='meters'
    height[:]=h
    wet=fout.createVariable('wet','f8',('ny','nx'))
    wet.units='none'
    wet[:]=np.where(h<0.,1.0,0.0)
    
    h_std=fout.createVariable('h_std','f8',('ny','nx'))
    h_std.units='meters'
    h_std[:]=hstd
    h_min=fout.createVariable('h_min','f8',('ny','nx'))
    h_min.units='meters'
    h_min[:]=hmin
    h_max=fout.createVariable('h_max','f8',('ny','nx'))
    h_max.units='meters'
    h_max[:]=hmax
    x=fout.createVariable('x','f8',('ny','nx'))
    x.units='meters'
    x[:]=xx
    y=fout.createVariable('y','f8',('ny','nx'))
    y.units='meters'
    y[:]=yy
    #global attributes
    if(not no_changing_meta):
    	fout.history = history
    	fout.description = description
    	fout.source =  source

    fout.sync()
    fout.close()

def get_indices1D_old(lon_grid,lat_grid,x,y):
    """This function returns the j,i indices for the grid point closest to the input lon,lat coordinates."""
    """It returns the j,i indices."""
    lons=np.fabs(lon_grid-x)
    lonm=np.where(lons==lons.min())
    lats=np.fabs(lat_grid-y)
    latm=np.where(lats==lats.min())
    j0=latm[0][0]
    i0=lonm[0][0]
#    print("wanted: ",x,y)
#    print("got:    ",lon_grid[i0] , lat_grid[j0])
#    print(j0,i0)
    return j0,i0
def mdist(x1,x2):
    """Returns positive distance modulo 360."""
    return np.minimum( np.mod(x1-x2,360.), np.mod(x2-x1,360.) )
def get_indices1D(lon_grid,lat_grid,x,y):
    """This function returns the j,i indices for the grid point closest to the input lon,lat coordinates."""
    """It returns the j,i indices."""
#    lons=np.fabs(lon_grid-x)
    lons=np.fabs(mdist(lon_grid,x))
    lonm=np.where(lons==lons.min())
    lats=np.fabs(lat_grid-y)
    latm=np.where(lats==lats.min())
    j0=latm[0][0]
    i0=lonm[0][0]
    print(" wanted: ",x,y)
    print(" got:    ",lon_grid[i0] , lat_grid[j0])
    good=False
    if(abs(x-lon_grid[i0]) < abs(lon_grid[1]-lon_grid[0])):
        good=True
        print("  good")
    else:
        print("  bad")
    print(" j,i=",j0,i0)
    return j0,i0,good

def get_indices2D(lon_grid,lat_grid,x,y):
    """This function returns the j,i indices for the grid point closest to the input lon,lat coordinates."""
    """It returns the j,i indices."""
    lons=np.fabs(lon_grid-x)
    lonm=np.where(lons==lons.min())
    lats=np.fabs(lat_grid-y)
    latm=np.where(lats==lats.min())
    j0=latm[0][0]
    i0=lonm[1][0]
#    print("wanted: ",x,y)
#    print("got:    ",lon_grid[j0,i0] , lat_grid[j0,i0])
#    print(j0,i0)
    return j0,i0
#Gibraltar
#wanted:  32.0 -12.5
#got:     31.9958333333 -12.5041666667
#9299 25439
#Gibraltar
#wanted:  40.7 4.7
#got:     40.6958333333 4.69583333333
#11363 26483
#Black sea
#wanted:  44.0 36
#got:     43.9958333333 36.0041666667
#15120 26879

def plot():
    pl.clf()
    display.clear_output(wait=True)
    plt.figure(figsize=(10,10))
    lonc = (lon[0,0]+lon[0,-1])/2
    latc = 90
    ax = plt.subplot(111, projection=cartopy.crs.NearsidePerspective(central_longitude=lonc, central_latitude=latc))
    ax.set_global()
    ax.stock_img()
    ax.coastlines()
    ax.gridlines()
    target_mesh.plot(ax,subsample=100, transform=cartopy.crs.Geodetic())
    plt.show(block=False)
    plt.pause(1)
    display.display(pl.gcf())

def refine_by_repeat(x,rf):
    xrf=np.repeat(np.repeat(x[:,:],rf,axis=0),rf,axis=1) #refine by repeating values
    return xrf

def extend_by_zeros(x,shape):
    ext=np.zeros(shape)
    ext[:x.shape[0],:x.shape[1]] = x
    return ext

def do_block(part,lon,lat,topo_lons,topo_lats,topo_elvs, max_mb=8000):
    print("  Doing block number ",part)
    print("  Target sub mesh shape: ",lon.shape)

    target_mesh = GMesh.GMesh( lon=lon, lat=lat )

    #plot()

    # Indices in topographic data
    ti,tj = target_mesh.find_nn_uniform_source( topo_lons, topo_lats )

    #Sample every other source points 
    ##Niki: This is only for efficeincy and we want to remove the constraint for the final product.
    ##Niki: But in some cases it may not work!
    #tis,tjs = slice(ti.min(), ti.max()+1,2), slice(tj.min(), tj.max()+1,2)
    tis,tjs = slice(ti.min(), ti.max()+1,1), slice(tj.min(), tj.max()+1,1)
    print('  Slices j,i:', tjs, tis )

    # Read elevation data
    topo_elv = topo_elvs[tjs,tis]
    # Extract appropriate coordinates
    topo_lon = topo_lons[tis]
    topo_lat = topo_lats[tjs]

    print('  Topo shape:', topo_elv.shape)
    print('  topography longitude range:',topo_lon.min(),topo_lon.max())
    print('  topography latitude  range:',topo_lat.min(),topo_lat.max())

    print("  Target     longitude range:", lon.min(),lon.max())
    print("  Target     latitude  range:", lat.min(),lat.max())

    # Refine grid by 2 till all source points are hit 
    print("  Refining the target to hit all source points ...")
    Glist = target_mesh.refine_loop( topo_lon, topo_lat , max_mb=max_mb);
    hits = Glist[-1].source_hits( topo_lon, topo_lat )
    print("  non-hit ratio: ",hits.size-hits.sum().astype(int)," / ",hits.size)

    # Sample the topography on the refined grid
    print("  Sampling the source points on target mesh ...")
    Glist[-1].sample_source_data_on_target_mesh(topo_lon,topo_lat,topo_elv)
    print("  Sampling finished...")

    # Coarsen back to the original taget grid
    print("  Coarsening back to the original taget grid ...")
    for i in reversed(range(1,len(Glist))):   # 1, makes it stop at element 1 rather than 0
        Glist[i].coarsenby2(Glist[i-1])
    
    print("Roughness calculation via plane fit")
    #Roughness calculation by plane fitting
    #Calculate the slopes of the planes on the coarsest (model) grid cells
    G=Glist[0]
    denom=(G.xxm-G.xm*G.xm)*(G.yym-G.ym*G.ym)-(G.xym-G.xm*G.ym)*(G.xym-G.xm*G.ym)
    alphd=(G.xzm-G.xm*G.zm)*(G.yym-G.ym*G.ym)-(G.yzm-G.ym*G.zm)*(G.xym-G.xm*G.ym)
    betad=(G.yzm-G.ym*G.zm)*(G.xxm-G.xm*G.xm)-(G.xzm-G.xm*G.zm)*(G.xym-G.xm*G.ym)
    #alph = alphd/denom
    #beta = betad/denom

    rf=2**(len(Glist)-1) #refinement factor
    #Generate the refined arrays from coarse arrays by repeating the coarse elements rf times
    #These arrays have the same values on finest mesh points inside each coarse cell by construction.
    #They are being used to calculate the (least-square) distance of data points
    #inside that cell from the fitted plane in that cell.
    xmrf=refine_by_repeat(G.xm,rf)
    ymrf=refine_by_repeat(G.ym,rf)
    zmrf=refine_by_repeat(G.zm,rf)
    alphdrf=refine_by_repeat(alphd,rf)
    betadrf=refine_by_repeat(betad,rf)
    denomrf=refine_by_repeat(denom,rf)
    #The refined mesh has a shape of (2*nj-1,2*ni-1) rather than (2*nj,2*ni) and hence
    #is missing the last row/column by construction! 
    #So, the finest mesh does not have (rf*nj,rf*ni) points but is smaller by ...
    #Bring it to the same shape as (rf*nj,rf*ni) by padding with zeros.
    #This is for algorithm convenience and we remove the contribution of them later.
    xs=extend_by_zeros(Glist[-1].xm,zmrf.shape)
    ys=extend_by_zeros(Glist[-1].ym,zmrf.shape)
    zs=extend_by_zeros(Glist[-1].zm,zmrf.shape)
    #Calculate the vertical distance D of each source point from the least-square plane
    #Note that the least-square plane passes through the mean data point.
    #The last rf rows and columns are for padding and denom is not zero on them.
    #To avoid division by zero calculate denom*D instead
    D_times_denom=denomrf*(zs-zmrf) - alphdrf*(xs-xmrf) - betadrf*(ys-ymrf)
    #Calculate topography roughness as the standard deviation of D on each coarse (model) grid cell
    #This is why we wanted to have a (nj*rf,ni*rf) shape arrays and padded with zeros above.
    D_times_denom_coarse=np.reshape(D_times_denom,(G.xm.shape[0],rf,G.xm.shape[1],rf))
    D_times_denom_coarse_std = D_times_denom_coarse.std(axis=(1,3))
    D_std=np.zeros(G.zm.shape)
    epsilon=1.0e-20 #To avoid negative underflow
    D_std[:,:] = D_times_denom_coarse_std[:,:]/(denom[:,:]+epsilon)

    print("")
    #print("Writing ...")
    #filename = 'topog_refsamp_BP.nc'+str(b) 
    #write_topog(Glist[0].height,fnam=filename,no_changing_meta=True)
    #print("haigts shape:", lons[b].shape,Hlist[b].shape)
    return Glist[0].height,D_std,Glist[0].h_min,Glist[0].h_max, hits


def usage(scriptbasename):
    print(scriptbasename + ' --hgridfilename <input_hgrid_filepath> --outputfilename <output_topog_filepath>  [--plot --no_changing_meta --open_channels]')


def main(argv):
    import socket
    import time
    host = str(socket.gethostname())
    scriptpath = sys.argv[0]
    scriptbasename = subprocess.check_output("basename "+ scriptpath,shell=True).decode('ascii').rstrip("\n")
    scriptdirname = subprocess.check_output("dirname "+ scriptpath,shell=True).decode('ascii').rstrip("\n")

    plotem = False
    open_channels = False
    no_changing_meta = False
    # URL of topographic data, names of longitude, latitude and elevation variables
    url,vx,vy,ve = '/work/Niki.Zadeh/datasets/topography/GEBCO_2014_2D.nc','lon','lat','elevation'
    # url,vx,vy,ve = 'http://thredds.socib.es/thredds/dodsC/ancillary_data/bathymetry/MED_GEBCO_30sec.nc','lon','lat','elevation'
    # url,vx,vy,ve = 'http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NGDC/.ETOPO1/.z_bedrock/dods','lon','lat','z_bedrock'
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["hgridfilename=","outputfilename=","no_changing_meta","open_channels","source_file=","source_lon=","source_lat=","source_elv="])
    except getopt.GetoptError as err:
        print(err)
        usage(scriptbasename)
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-i", "--hgridfilename"):
            gridfilename = arg
        elif opt in ("-o", "--outputfilename"):
            outputfilename = arg
        elif opt in ("--source_file"):
            url = arg
        elif opt in ("--source_lon"):
            vx = arg
        elif opt in ("--source_lat"):
            vy = arg
        elif opt in ("--source_elv"):
            ve = arg
        elif opt in ("--plot"):
            plotem = True
        elif opt in ("--no_changing_meta"):
            no_changing_meta = True
        elif opt in ("--open_channels"):
            open_channels = True
        else:
            assert False, "unhandled option"


    print("")
    print("Generatin model topography for target grid ", gridfilename)
    #Information to write in file as metadata
    scriptgithash = subprocess.check_output("cd "+scriptdirname +";git rev-parse HEAD; exit 0",stderr=subprocess.STDOUT,shell=True).decode('ascii').rstrip("\n")
    scriptgitMod  = subprocess.check_output("cd "+scriptdirname +";git status --porcelain "+scriptbasename+" | awk '{print $1}' ; exit 0",stderr=subprocess.STDOUT,shell=True).decode('ascii').rstrip("\n")
    if("M" in str(scriptgitMod)):
        scriptgitMod = " , But was localy Modified!"

    hist = "This file was generated via command " + ' '.join(sys.argv)
    if(not no_changing_meta):
        hist = hist + " on "+ str(datetime.date.today()) + " on platform "+ host

    desc = "This is a model topography file generated by the refine-sampling method from source topography. "
    
    source =""
    if(not no_changing_meta):
        source =  source + scriptpath + " had git hash " + scriptgithash + scriptgitMod 
        source =  source + ". To obtain the grid generating code do: git clone  https://github.com/nikizadehgfdl/thin-wall-topography.git ; cd thin-wall-topography;  git checkout "+scriptgithash

    #Time it
    tic = time.perf_counter()
    # # Open and read the topographic dataset
    # Open a topography dataset, check that the topography is on a uniform grid.
    topo_data = netCDF4.Dataset(url)

    # Read coordinates of topography
    topo_lons = np.array( topo_data.variables[vx][:] )
    topo_lats = np.array( topo_data.variables[vy][:] )
    topo_elvs = np.array( topo_data.variables[ve][:,:] )

    #Fix the topography to open some channels
    if(open_channels):
        #Bosporus mouth at Marmara Sea (29.03,41.04)
        j0,i0=15724,39483 #get_indices1D(topo_lons, topo_lats ,29.03, 41.04)
        #One grid cell thick (not survived ice9)
        #topo_elvs[j0,i0]=topo_elvs[j0,i0-1]
        #topo_elvs[j0+1,i0+2]=topo_elvs[j0+1,i0+1]
        #topo_elvs[j0+3,i0+3]=topo_elvs[j0+3,i0+2]
        #wide channel
        j2,i2=15756, 39492 #get_indices1D(topo_lons, topo_lats ,29.1, 41.3)
        topo_elvs[j0-10:j2,i0-10:i2+10]=topo_elvs[j0,i0-1]

        #Dardanells' constrict
        j1,i1=15616, 39166 #get_indices1D(topo_lons, topo_lats ,26.39, 40.14)
        topo_elvs[j1+1,i1]=topo_elvs[j1,i1]
    #Read a target grid
    targ_grid =  netCDF4.Dataset(gridfilename)
    targ_lon = np.array(targ_grid.variables['x'])
    targ_lat = np.array(targ_grid.variables['y'])
    #x and y have shape (nyp,nxp). Topog does not need the last col for global grids (period in x). 
    targ_lon = targ_lon[:,:-1]
    targ_lat = targ_lat[:,:-1]
    print(" Target mesh shape: ",targ_lon.shape)
    #Translate topo data to start at target_mesh.lon_m[0]
    #Why/When?
    jllc,illc,status1=get_indices1D(topo_lons, topo_lats ,targ_lon[0,0] ,targ_lat[0,0])
    jurc,iurc,status2=get_indices1D(topo_lons, topo_lats ,targ_lon[0,-1],targ_lat[-1,0])
    if(not status1 or not status2):
        print(' shifting topo data to start at target lon')
        topo_lons = np.roll(topo_lons,-illc,axis=0) #Roll data longitude to right
        topo_lons = np.where(topo_lons>=topo_lons[0] , topo_lons-360, topo_lons) #Rename (0,60) as (-300,-180) 
        topo_elvs = np.roll(topo_elvs,-illc,axis=1) #Roll data depth to the right by the same amount.

    print(' topography grid array shapes: ' , topo_lons.shape,topo_lats.shape)
    print(' topography longitude range:',topo_lons.min(),topo_lons.max())
    print(' topography longitude range:',topo_lons[0],topo_lons[-1000])
    print(' topography latitude range:',topo_lats.min(),topo_lats.max())
    print(' Is mesh uniform?', GMesh.is_mesh_uniform( topo_lons, topo_lats ) )
    ### Partition the Target grid into non-intersecting blocks
    #This works only if the target mesh is "regular"! Niki: Find the mathematical buzzword for "regular"!!
    #Is this a regular mesh?
    # if( .NOT. is_mesh_regular() ) throw

    #Niki: Why 4,1 partition?
    xb=4
    yb=1
    lons=break_array_to_blocks(targ_lon,xb,yb)
    lats=break_array_to_blocks(targ_lat,xb,yb)

    #We must loop over the 4 partitions
    Hlist=[]
    Hstdlist=[]
    Hminlist=[]
    Hmaxlist=[]
    for part in range(0,xb):
        lon = lons[part]
        lat = lats[part]
        h,hstd,hmin,hmax,hits = do_block(part,lon,lat,topo_lons,topo_lats,topo_elvs)
        Hlist.append(h)
        Hstdlist.append(hstd)
        Hminlist.append(hmin)
        Hmaxlist.append(hmax)

    print(" Merging the blocks ...")
    height_refsamp = undo_break_array_to_blocks(Hlist,xb,yb)
    hstd_refsamp = undo_break_array_to_blocks(Hstdlist,xb,yb)
    hmin_refsamp = undo_break_array_to_blocks(Hminlist,xb,yb)
    hmax_refsamp = undo_break_array_to_blocks(Hmaxlist,xb,yb)
    write_topog(height_refsamp,hstd_refsamp,hmin_refsamp,hmax_refsamp,targ_lon,targ_lat,fnam=outputfilename,description=desc,history=hist,source=source,no_changing_meta=no_changing_meta)

    #Niki: Why isn't h periodic in x?  I.e., height_refsamp[:,0] != height_refsamp[:,-1]
    print(" Periodicity test  : ", height_refsamp[0,0] , height_refsamp[0,-1])
    print(" Periodicity break : ", (np.abs(height_refsamp[:,0]- height_refsamp[:,-1])).max() )
    toc = time.perf_counter()
    print("It took {toc - tic:0.4f} seconds on platform ",host)

    if(plotem):
        import matplotlib.pyplot as plt
        import pylab as pl
        from IPython import display
        import cartopy

        plt.figure(figsize=(10,10))
        ax = plt.subplot(111, projection=cartopy.crs.NearsidePerspective(central_latitude=90))
        ax.set_global()
        ax.stock_img()
        ax.coastlines()
        ax.gridlines()
        im = ax.pcolormesh(targ_lon,targ_lat,height_refsamp, transform=cartopy.crs.PlateCarree())
        plt.colorbar(im,ax=ax);


if __name__ == "__main__":
    main(sys.argv[1:])



#B ny,nx=  118 720 , st ny,nx=  119 721
#M ny,nx=  350 720 , st ny,nx=  351 721
#S ny,nx=   56 720 , st ny,nx=   57 721

#T ny,nx=  524 720
