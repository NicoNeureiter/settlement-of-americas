import numpy
import collections
import cartopy.geodesic as geodesic
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

boundingbox = collections.namedtuple("boundingbox",
    ["w", "e", "s", "n"])
point = collections.namedtuple("point",
    ["longitude", "latitude"])

australia = boundingbox(
    112.8708,
    153.6392,
    -43.8615,
    -9.6712)

continent = australia

bbox_centre = point(
    (continent.e + continent.w)/2,
    (continent.n + continent.s)/2)

area = 450000000 #m²
# Assuming ideal hexagons,
SQRT3 = 3 ** 0.5
hexagon_side = (2 * SQRT3 * area) ** 0.5 / 3
grid_point_distance = SQRT3 * hexagon_side
assert abs((grid_point_distance * hexagon_side * 1.5) - area) < 1e-2

geodesic = geodesic.Geodesic()

points = [bbox_centre]

# FIXME: This breaks when the bounding box crosses the date line.
# TODO: Write a regression test for that case.

# Neighbors east and west: Direct tiling
while points[-1].longitude < continent.e:
    next = geodesic.direct(points[-1], 90, grid_point_distance)
    points.append(point(*numpy.array(next)[0, :2]))
while points[0].longitude > continent.w:
    next = geodesic.direct(points[0], 270, grid_point_distance)
    points.insert(0, point(*numpy.array(next)[0, :2]))

grid = numpy.array([points])
while (grid[0, :, 1] < continent.n).any():
    next = geodesic.direct(grid[0], 0, 3 * hexagon_side)
    grid = numpy.vstack((
        [numpy.array(next)[:, :2]],
        grid))
while (grid[-1, :, 1] > continent.s).any():
    next = geodesic.direct(grid[-1], 180, 3 * hexagon_side)
    grid = numpy.vstack((
        grid,
        [numpy.array(next)[:, :2]]))

mid_grid = (grid[:-1, :-1, :] + grid[1:, 1:, :]) / 2

land_shp_fname = shpreader.natural_earth(
    resolution='10m', category='physical', name='land')

land_geom = unary_union(
    [record.geometry
     for record in shpreader.Reader(land_shp_fname).records()
     if record.attributes.get('featurecla') != "Null island"])

land = prep(land_geom)

def is_land(x, y):
   return land.contains(sgeom.Point(x, y))

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines("50m")
# ax.set_extent(continent)
plt.scatter(grid[..., 0], grid[..., 1],
            c="r", s=1)
plt.scatter(mid_grid[..., 0], mid_grid[..., 1],
            c="r", s=1)

plt.show()
