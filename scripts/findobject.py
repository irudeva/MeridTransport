data = """\
000000011111111000000000000000000000000
000000011111111000000000000000000000000
000000011111111000000000000000000000000
000000001111110000000001000000000110000
000000000111110000000011000000001111100
000000000011100000000000100000011111100
000000000000000000000000000000011111100
000000000000000000000000000000011111110
000000000000000000000000000000111111110
000000000000000000000000000000111111110
000000000000000000000000000000011111100
000000000000000000000000000000001000000
000000000000000000000000000000000000000"""

from collections import namedtuple
Point = namedtuple('Point', 'x y')

def points_adjoin(p1, p2):
    # to accept diagonal adjacency, use this form
    #return -1 <= p1.x-p2.x <= 1 and -1 <= p1.y-p2.y <= 1
    return (-1 <= p1.x-p2.x <= 1 and p1.y == p2.y or
             p1.x == p2.x and -1 <= p1.y-p2.y <= 1)

def adjoins(pts, pt):
    return any(points_adjoin(p,pt) for p in pts)

def locate_regions(datastring):
    data = map(list, datastring.splitlines())
    regions = []
    datapts = [Point(x,y)
                for y,row in enumerate(data)
                    for x,value in enumerate(row) if value=='1']
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


# test against data set
regs = locate_regions(data)
print len(regs)
print_regions(regs)
