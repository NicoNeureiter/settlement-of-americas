import numpy
import collections
import cartopy.geodesic as geodesic
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.collections
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
import tifffile
import matplotlib.transforms as mtransforms

GEODESIC = geodesic.Geodesic()
SQRT3 = 3 ** 0.5

land_shp_fname = shpreader.natural_earth(
    resolution='50m', category='physical', name='land')

land_geom = unary_union(
    [record.geometry
     for record in shpreader.Reader(land_shp_fname).records()
     if record.attributes.get('featurecla') != "Null island"])

LAND = prep(land_geom)


# Define classes
BoundingBox = collections.namedtuple("boundingbox",
    ["w", "e", "s", "n"])
Point = collections.namedtuple("point",
    ["longitude", "latitude"])


def hexagon(area):
    """Calculate the edge length and midpoint distance of a plane hexagon.

    Given a hexagon of area `area`, return its edge length and the distance
    between closest neighbour mid points in a hexagonal grid.

    Parameters
    ==========
    area: number
        The area of a plane hexagon

    Returns
    =======
    hexagon_side: float
        The edge length of the hexagon
    grid_point_distance: float
        The distance beween the mid points of two adjacent hexagons in a hexagonal tiling
    """
    hexagon_side = (2 * SQRT3 * area) ** 0.5 / 3
    grid_point_distance = SQRT3 * hexagon_side
    assert abs((grid_point_distance * hexagon_side * 1.5) - area) < 1e-2
    return hexagon_side, grid_point_distance


def hexagonal_earth_grid(bbox, area):
    """Generate a hexagonal grid.

    Generate a spherical hexagonal grid inside the bounding box in which each
    hexagon has roughly the given average area.

    In the return value, mid_grid[k, l] has neighbours mid_grid[k-1, l],
    mid_grid[k+1, l], grid[k, l], grid[k+1, l], grid[k, l+1], grid[k+1, l+1]

    Parameters
    ==========
    area: float
        The desired hexagon surface area in m²
    bbox: BoundingBox
        The w/e/s/n bounding box to fill with the grid

    Returns
    =======
    grid: numpy.array, shape=(M, N, 2)
        one rectangular sub-grid of the hexagonal grid
    mid_grid:
        the complementary rectangular sub-grid of the hexagonal grid

    """
    bbox_centre = Point(
        (continent.e + continent.w)/2,
        (continent.n + continent.s)/2)
    hexagon_side, grid_point_distance = hexagon(area)

    points = [bbox_centre]

    # FIXME: This breaks when the bounding box crosses the date line.
    # TODO: Write a regression test for that case.

    # Neighbors east and west: Direct tiling
    while points[-1].longitude < continent.e:
        next = GEODESIC.direct(points[-1], 90, grid_point_distance)
        points.append(Point(*numpy.array(next)[0, :2]))
    while points[0].longitude > continent.w:
        next = GEODESIC.direct(points[0], 270, grid_point_distance)
        points.insert(0, Point(*numpy.array(next)[0, :2]))

    grid = numpy.array([points])
    while (grid[0, :, 1] < continent.n).any():
        next = GEODESIC.direct(grid[0], 0, 3 * hexagon_side)
        grid = numpy.vstack((
            [numpy.array(next)[:, :2]],
            grid))
    while (grid[-1, :, 1] > continent.s).any():
        next = GEODESIC.direct(grid[-1], 180, 3 * hexagon_side)
        grid = numpy.vstack((
            grid,
            [numpy.array(next)[:, :2]]))

    mid_grid = (grid[:-1, :-1, :] + grid[1:, 1:, :]) / 2

    return grid, mid_grid


# Generate a hexagonal grid over the continent
# Define continents
australia = BoundingBox(
    112.8708,
    153.7392,
    -43.8615,
    -9.6712)

continent = australia

area = 450000000 #m²

grid = hexagonal_earth_grid(continent, area)

def is_land(xy):
   return LAND.contains(sgeom.Point(*xy))

land = (
    numpy.apply_along_axis(is_land, axis=2, arr=grid[0]),
    numpy.apply_along_axis(is_land, axis=2, arr=grid[1]))

try:
    precipitation
except NameError:
    precipitation = tifffile.imread("../worldclim/wc2.0_bio_30s_12.tif").clip(0)

def coordinates_to_index(points, resolution=2 * 60):
    """Convert long,lat coordinate pairs into indices in a TIF

    Convert a [..., 2] ndarray, or a pair of coordinates, into the matching
    grid indices of a Mercator projection pixel map with a given resolution in
    pixels per degree.

    Paramaters
    ==========
    points: ndarray-like, shape=(..., 2)
        An array of longitude. latitude pairs to be converted to grid indices

    resolution:
        The resolution of the grid in indices per degree

    Returns
    =======
    ndarray(int), shape=(..., 2)
        An integer array of grid indices

    """
    points = numpy.asarray(points)
    return numpy.stack(
        (numpy.round((-points[..., 1] + 90) * resolution).astype(int),
         numpy.round((points[..., 0] + 180) * resolution).astype(int)),
        -1)

def index_to_coordinates(indices, resolution=2 * 60):
    """Convert grid indices into long,lat coordinate pairs

    Convert a (ndarray of) grid indices of a Mercator projection pixel map with
    a given resolution in pixels per degree into geocoordinate pairs (longitude,
    latitude).

    Paramaters
    ==========
    indices: ndarray(int), shape=(..., 2)
        An integer array of grid indices

    resolution:
        The resolution of the grid in indices per degree

    Returns
    =======
    ndarray, shape=(..., 2)
        An array of longitude. latitude pairs

    """
    indices = numpy.asarray(indices)
    return numpy.stack(
        (indices[..., 1] / resolution - 180,
         90 - indices[..., 0] / resolution),
        -1)

# TODO: Write tests for middle, random, each corner, forwards and backwards.

all_gridcells = {}
def gridcell(m, i, j):
    if i < 0 or i >= grid[m].shape[0]:
        return None
    if j < 0 or j >= grid[m].shape[1]:
        return None
    try:
        return all_gridcells[m, i, j]
    except KeyError:
        all_gridcells[m, i, j] = GridCell(m, i, j)
        return all_gridcells[m, i, j]

class GridCell():
    alpha = 10 ** -8.07
    beta = 2.64
    grid = hexagonal_earth_grid(continent, area)

    def __init__(self, m, i, j, grid=grid):
        self.m = m
        self.ij = i, j
        self.population = 0
        self.popcap = self.population_capacity() * area / 1000000 # area is in m², popcap is per km²
        self.language = None

    def polygon(self):
        try:
            return self._polygon
        except AttributeError:
            neighbors = numpy.array([n.point for n in self.neighbors(True, True)])
            self._polygon = (neighbors + numpy.roll(neighbors, 1, 0) + self.point) / 3
            return self._polygon

    @property
    def point(self):
        return Point(*grid[self.m][self.ij])

    def population_capacity(self):
        """Calculate the pop cap of a cell given its precipitation

        Return the carrying capacity K of a hexagonal cell of area AREA with
        the given mean yearly precipitation measured in mm

        In Gavin et al. 2017, the function used is

        K = α * P ** β

        with eg. α = 10 ** -8.07, β = 2.64 [Note: Their supplementary material
        lists α=-8.07, but that is a nonsensical number. The actual plot line
        does not exactly correspond to α=10^-8.07, but more closely to
        α=10^-7.96, but that suggests that this is at least close to the
        described behaviour.]

        Parameters
        ==========
        precipitation: float
            The cell's mean yearly precipitation P, in mm

        Returns
        =======
        float
            The cell's carrying capacity, in individuals/km^2
        """
        return self.alpha * self.precipitation() ** self.beta

    def __hash__(self):
        return hash((self.m, self.ij))

    def precipitation(self):
        index = tuple(coordinates_to_index(self.point))
        return precipitation[index]

    def neighbors(self, include_unlivable=False, include_foreign=False):
        i, j = self.ij
        if self.m==0:
            neighbors = [
                gridcell(0, i, j+1),
                gridcell(1, i, j),
                gridcell(1, i, j-1),
                gridcell(0, i, j-1),
                gridcell(1, i-1, j-1),
                gridcell(1, i-1, j),
            ]
        else:
            neighbors = [
                gridcell(1, i, j+1),
                gridcell(0, i, j+1),
                gridcell(0, i, j),
                gridcell(1, i, j-1),
                gridcell(0, i+1, j),
                gridcell(0, i+1, j+1),
            ]
        neighbors = [g or self for g in neighbors]
        return [g for g in neighbors
                if include_foreign or g.language == self.language or g.language is None
                if include_unlivable or g.popcap >= 1]

growth_factor = 1.1

import csv

popcaps = []
with open("../binford/data-raw/HG_Output_final.csv") as binford_data:
    for people in csv.DictReader(binford_data, skipinitialspace=True):
        if people["STATE"].startswith("AUSTRALIA"):
            # Exclude AUS populations
            continue
        if int(people["CLIM"]) < 3:
            # Exclude arctic and sub-arctic populations
            continue
        popcaps.append(float(people["TLPOP"]))

class Language:
    def __init__(self, color, cell):
        g = gridcell(*cell)
        g.population = 10
        g.language = self
        self.cells = {g}
        self.id = color
        self.popcap = popcaps[numpy.random.randint(len(popcaps))]

    def grow(self):
        grow_into = {}
        growth = 0
        for cell in self.cells:
            if cell not in grow_into:
                grow_into[cell] = cell.popcap - cell.population
            region_popcap = cell.popcap
            region_population = cell.population
            for n in cell.neighbors():
                if n not in grow_into:
                    grow_into[n] = n.popcap - n.population
                region_popcap += n.popcap
                region_population += n.population
            growth += growth_factor * cell.population * (1 - region_population / region_popcap)
        if not set(grow_into) - self.cells:
            raise StopIteration
        self.popcap -= growth
        if self.popcap < 0:
            raise StopIteration
        distribution = growth / sum(grow_into.values())
        if distribution > 1:
            print("Population will grow by more than its population capacity allows: {:}".format(distribution))
        for cell, proportion in grow_into.items():
            if cell.popcap < 1:
                continue
            cell.language = self
            cell.population += proportion * distribution
        self.cells = self.cells.union(grow_into)
        print(growth)

def plot_grid(function):
    polygons = []
    values = []
    for cell in list(all_gridcells.values()):
        polygons.append(cell.polygon())
        values.append(function(cell))
    collection = matplotlib.collections.PolyCollection(
        polygons, facecolors=values)

    return collection

def random_cell():
    m = numpy.random.randint(2)
    return (
        m,
        numpy.random.randint(grid[m].shape[0]),
        numpy.random.randint(grid[m].shape[1]))

g = random_cell()
while gridcell(*g).popcap < 1:
    g = random_cell()

l = Language((numpy.random.random(), numpy.random.random(), numpy.random.random(), 0.8),
             g)
while True:
    while True:
        try:
            l.grow()
        except StopIteration:
            break
    free = [i for i, g in all_gridcells.items()
            if g.language is None
            if g.popcap >= 1]
    if not free:
        break
    l = Language((numpy.random.random(), numpy.random.random(), numpy.random.random(), 0.8),
                 free[numpy.random.randint(len(free))])



ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines("50m")
ax.set_extent(continent)

ax.add_collection(plot_grid(lambda cell: (cell.language.id if cell.language else (0,0,0,1))))
plt.show()

