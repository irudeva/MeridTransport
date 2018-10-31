from netCDF4 import Dataset
from collections import namedtuple
import numpy as np
import datetime as datetime


def leapyear(year):
    if year/400 == 0 :
        return False
    if year%100 == 0 :
        return False
    if year%4 == 0 :
        return True

Point = namedtuple('Point', 'x y')

def points_adjoin(p1, p2):
    # to accept diagonal adjacency, use this form
    #return -1 <= p1.x-p2.x <= 1 and -1 <= p1.y-p2.y <= 1
    return (-1 <= p1.x-p2.x <= 1 and p1.y == p2.y or
             p1.x == p2.x and -1 <= p1.y-p2.y <= 1)

def adjoins(pts, pt):
    return any(points_adjoin(p,pt) for p in pts)

def locate_regions(datastring):
    # data = map(list, datastring.splitlines())
    data = datastring
    regions = []
    datapts = [Point(x,y)
                for y,row in enumerate(data)
                    for x,value in enumerate(row) if value<=-2]
    for dp in datapts:
        # find all adjoining regions
        adjregs = [r for r in regions if adjoins(r,dp)]
        if adjregs:
            adjregs[0].add(dp)
            if len(adjregs) > 1:
                # joining more than one reg, merge
                regions[:] = [r for r in regions if r not in adjregs]
                regions.append(reduce(set.union, adjregs))
        else:
            # not adjoining any, start a new region
            regions.append(set([dp]))
    return regions

def region_index(regs, p):
    return next((i for i,reg in enumerate(regs) if p in reg), -1)

def print_regions(regs):
    maxx = max(p.x for r in regs for p in r)
    maxy = max(p.y for r in regs for p in r)
    allregionpts = reduce(set.union, regs)
    for y in range(-1,maxy+2):
        line = []
        for x in range(-1,maxx+2):
            p = Point(x, y)
            if p in allregionpts:
                line.append(str(region_index(regs, p)))
            else:
                line.append('.')
        print ''.join(line)
    print


# yr1 = 1983
# yr2 = 2017
#
# # clmon = np.array([3,9])
# # chmon = np.array(['Mar','Sep'], dtype=str)
# # clmon = np.array([3])
# # chmon = np.array(['Mar'], dtype=str)
# clmon = np.array([1])
# chmon = np.array(['Jan'], dtype=str)


yr = 1997
mm = 12
dd = 15
hr = 00
# t = 0
hs ="SH"


#----reading flux-------------------------------------------------------------
print("reading flux...")

fin = "/Users/irudeva/work/Projects/MeridTransport/data/flux.LH_Teddies.%d.nc"%(yr)
fout = "/Users/irudeva/work/Projects/MeridTransport/output/objflux.LH_Teddies.%d.nc"%(yr)


# for iyr in range(yr,yr+1):
    # print iyr
    # fout = "../output/wind/frwnd6h_gt%d.%s%d.erain.nc"%(XXwnd,chmon[0],iyr)


    # fin = "/Users/irudeva/work/DATA/ERAint/Plev/erain.hgt_air_wind.6h.%d.nc"%(iyr)
    # print "read file %s"%(fin)


dimnam=('time','lon','lat')
varnam=['time','lon','lat','Tr_eddies']

# Read zonal and meridional wind components from file using the netCDF4
# module. The components are defined on pressure levels and are in separate
# files.
nc = Dataset(fin, 'r')
v=0
for var in varnam:
    if nc.variables[varnam[v]].name != var:
        print "Variables don't agree", var, nc.variables[varnam[v]].name, v
        exit()
    v += 1

lons = nc.variables[varnam[1]][:]
lats = nc.variables[varnam[2]][:]
# if iyr == yr1:
#     # tmp = nc.variables[varnam[3]][:]
#     time = np.zeros((nyr,1464))
#     uwnd = np.zeros((nyr,1464,lev.size,lats.size,lons.size))
#     vwnd = np.zeros((nyr,1464,lev.size,lats.size,lons.size))

time = nc.variables[varnam[0]][:]
var = nc.variables[varnam[3]][:]
# del tmp

# nc.close()

print("flux data read")

# quit()
#----selecting timesteps-------------------------------------------------------------
print "creating time mask..."
btime = np.zeros_like(time,dtype = bool)
bvar = np.zeros_like(var,dtype = bool)

# for iyr in range(yr1,yr2+1):
#     print iyr

# tsize = 1460
# if leapyear(iyr):
#     tsize =1464
# print tsize
# dt_time = [datetime.datetime(1900, 1, 1) + datetime.timedelta(hours=int(t))\
#        for t in time[iyr-yr1,:tsize]]
dt_time = [datetime.datetime(1900, 1, 1) + datetime.timedelta(hours=int(t))\
       for t in time]


# print dt_time
# quit()
# datetime.date(iyr,m,1)

# for m in mm :
    # yr1 = iyr
for t in dt_time :
    if t.month == mm and t.day == dd and t.hour == hr:
        print 'selected time: ', t
        ind = dt_time.index(t)
        # print ind
        # btime[ind] = 1
        # bvar[ind,:,:,:] = 1

del dt_time


var1 = var[ind,:,:]*10**-8
# print(var1*10**-8)
# print(np.amax(var1*10**-8, axis=0))
# print(np.amin(var1*10**-8, axis=0))

regs = locate_regions(var1)
print len(regs)
print_regions(regs)

print(regs)


#---NetCDF write---------------------------------------------------------------
print("Start NetCDF writing")

ncout = Dataset(fout, 'w', format='NETCDF4')
ncout.description = "flux objects %s" % (fin)

nc = Dataset(fin, 'r')  # open last fin to copy attributes

# Using our previous dimension info, we can create the new time dimension
# Even though we know the size, we are going to set the size to unknown

ncout.createDimension(dimnam[1], lons.size)
ncout.createDimension(dimnam[2], lats.size)
ncout.createDimension(dimnam[0], None)

for nv in range(1, 3) :
    print(varnam[nv])
    ncout_var = ncout.createVariable(varnam[nv], nc.variables[varnam[nv]].dtype,dimnam[nv])
    for ncattr in nc.variables[varnam[nv]].ncattrs():
        ncout_var.setncattr(ncattr, nc.variables[varnam[nv]].getncattr(ncattr))
#print(nc.variables['latitude'].ncattrs())

ncout.variables[varnam[1]][:] = lons
ncout.variables[varnam[2]][:] = lats
#  for time:
# nv = 2  #
# ncout_var = ncout.createVariable(varnam[nv], 'i2',dimnam[nv])
# ncout_var.long_name = 'time'
# ncout_var.months = chmon[0]
# ncout_var[:] = np.arange(iyr,iyr+1)

ncout_var1 = ncout.createVariable('fluxobj', 'i2',dimnam[2:0:-1])
ncout_var1.long_name = 'objects_from_%sflux'%("LH")
# sf_scale = 1.e+7
# sf_add   = 0.
# ncout_sf.scale_factor = sf_scale
# ncout_sf.add_offset   = sf_add
# ncout_var1.units        = 'm s**-1'

print dimnam[2:0:-1]
print var1.shape
# print regs.shape

#!!!automatically takes scale and offset into account
#!!! no need for: ncout_sf[:] = (sf-sf_add)/sf_scale
ncout_var1[:] = var1

# ncout_var1 = ncout.createVariable('frwnd_vpos', 'i2',dimnam[2::-1])
# ncout_var1.long_name = 'frequency_wind_greater_%d_v_positive'%(XXwnd)
# ncout_var1.units        = 'm s**-1'
# ncout_var1[:] = frwndXX_vpos
#
# ncout_var1 = ncout.createVariable('frwnd_vneg', 'i2',dimnam[2::-1])
# ncout_var1.long_name = 'frequency_wind_greater_%d_v_negative'%(XXwnd)
# ncout_var1.units        = 'm s**-1'
# ncout_var1[:] = frwndXX_vneg


nc.close()
ncout.close()
