import numpy as np
from scipy.special import comb
from enum import Enum, auto
import matplotlib.pyplot as plt


class Point:
	def __init__(self, *args):
		msg = "Invalid arguments"
		self._point = np.array(args).flatten()

		if not len(self._point) == 2:
			raise ValueError(msg)
		elif not all([self._point.dtype.kind in ['i', 'f'] for e in self._point]):
			raise TypeError(msg)

	@property
	def point(self):
		return self._point

	@property
	def x(self):
		return self._point[0]

	@x.setter
	def x(self, _x):
		self._point[0] = _x

	@property
	def y(self):
		return self._point[1]

	@y.setter
	def y(self, _y):
		self._point[1] = _y
	

	def plot(self, ax = None, fmt = 'o'):
		if ax is None:
			_, ax = plt.subplots()
		ax.plot(self.x, self.y, fmt)


	def distance(self, point):
		"""
		Calcurate distance between two points in 2D.

		Parameters
		----------
		point : Point or number tuple or number list
			Target point coodinates.

		Returns
		-------
		two points distance : float
			Distance from this point to target point.
		"""

		if not isinstance(point, Point):
			try:
				point = Point(point[0], point[1])
			except:
				msg = "Convert failed"
				raise TypeError(msg)
		return np.linalg.norm(self.point - point.point)


class Bezier:
	def __init__(self, points):
		if not isinstance(points, list):
			msg = "We assume 'list' object but get {0}".format(type(points))
			raise TypeError(msg)

		if all(isinstance(e, Point) for e in points):
			pass
		else:
			points = [Point(e) for e in points]

		self._n = len(points)
		self._points = points

	@property
	def dims(self):
		return self._n - 1

	@property
	def plist(self):
		return [e.point for e in self._points]

	@property
	def points(self):
		return self._points

	@property
	def points4matplot(self):
		xs, ys = zip(*self.plist)
		return xs, ys


	def getConvexhullLines(self):
		"""
		Fetch Bezier curve's convex-hull lines list

		Parameters
		----------
		None.

		Returns
		-------
		convexhull_lines : list of PlaneLine object
			List of Convex-hull lines.
		
		Note
		-------
		Using gift wrapping algorithm.
		"""
		plst = self.plist
		## search most left and bottom point
		plst = sorted(plst, key=lambda x: x[0], reverse=True)
		plst = sorted(plst, key=lambda x: x[1])
		init_p, plst = plst[0], plst[1:]

		convexhull_lines = []
		current_p = init_p
		prev_p = None
		end = False
		while end == False:
			if convexhull_lines == []:
				# If first, make holizontal line
				base_line = PlaneLine([current_p, (current_p[0]+1, current_p[1])])
			else:
				base_line = PlaneLine([prev_p, current_p])

			rads = [base_line.radian(PlaneLine([current_p, p]),
			                         PlaneLine.DegRange.ZERO__2PI) for p in plst]
			
			last_index = None
			if not convexhull_lines == []:
				last_index = len(rads)
				rads.append(base_line.radian(PlaneLine([current_p, init_p]), PlaneLine.DegRange.ZERO__2PI))

			pop_index = np.array(rads).argmin()
			if pop_index == last_index:
				end = True
				nxt = init_p
			else:
				nxt = plst.pop(pop_index)
			convexhull_lines.append(PlaneLine([current_p, nxt]))
			prev_p = current_p
			current_p = nxt
		return convexhull_lines


	def _bernstein(self, n, i, t):
		return comb(n, i) * t**i * (1 - t)**(n-i)


	def bezier_point(self, t):
		"""
		Fetch the value of a point on a Bezier curve.

		Parameters
		----------
		t : float
			The value of the parameter t.

		Returns
		-------
		p : Point object
			(x, y) coodinates at parameter t.
		
		Note 
		-------
		t must be in the range 0 to 1.
		"""
		if t < 0 or 1 < t:
			msg = "t must be in the range 0 to 1"
			raise ValueError(msg)

		p = np.zeros(2)
		for i, f in enumerate(self._points):
			p += self._bernstein(self.dims, i, t) * f.point
		return Point(p)
	

	def create_Ux(self, dims=None):
		"""
		Make lower trianglar matrix for formula transformation.

		Parameters
		----------
		dims : int
			The value of the Bezier cirve dimensional.

		Returns
		-------
		Ux : ndarray
			Lower trianglar matrix.
		"""
		if dims is None:
			dims = self.dims

		U = np.zeros([dims + 1, dims + 1])
		lmbd = lambda n, i, j: comb(n, j) * comb(n-j, i-j) * (-1)**(i - j)

		for i in range(dims + 1):
			lst = list(range(i+1)) + [-1]*(dims-i)
			elems = [lmbd(dims, i, j) if j >= 0 else 0.0 for j in lst]
			mtx = np.array(elems)
			U[i] = mtx
		return U
	

	def split(self, t):
		"""
		Split Bezier curve at point 't'.

		Parameters
		----------
		t : float
			The value of the parameter t at split point.

		Returns
		-------
		(b1, b2) : Tuple of Bezier object
			Splited left Bezier curve and right Bezier curve.
		
		Note 
		-------
		t must be in the range 0 to 1.
		"""
		if t < 0 or 1 < t:
			msg = "t must be in the range 0 to 1"
			raise ValueError(msg)

		M = self.create_Ux()
		iM = np.linalg.inv(M)

		Z = np.eye(self.dims + 1)
		for i in range(self.dims + 1):
			Z[i] = Z[i] * t ** i 		
		Q = iM @ Z @ M
		
		xs, ys = self.points4matplot
		X = np.array(xs)
		Y = np.array(ys)

		left_bezier = Bezier(list(zip(Q @ X, Q @ Y)))

		_Q = np.zeros((self.dims + 1, self.dims + 1))
		lst = []
		for i in reversed(range(self.dims + 1)):
			l = [-1] * i + list(range(self.dims + 1 - i))
			lst.append(l)
		for i, l in enumerate(lst):
			mtx = [Q[i][e] if not e == -1 else 0 for e in l]
			_Q[i] = np.array(mtx)
		_Q = np.flipud(_Q)
		
		right_bezier = Bezier(list(zip(_Q @ X, _Q @ Y)))

		return left_bezier, right_bezier
			
	def plot(self, ax=None, with_ctrl_pt=True, bcolor="black", ccolor="gray", resolution=100):
		if ax is None:
			_, ax = plt.subplots()
		prev_point = None
		for t in np.linspace(0, 1, resolution):
			bp = self.bezier_point(t)
			if prev_point is None:
				prev_point = (bp.x, bp.y)
			xs, ys = zip(*(prev_point, (bp.x, bp.y)))
			ax.plot(xs, ys, '-', color=bcolor)
			prev_point = (bp.x, bp.y)

		if with_ctrl_pt:
			xs, ys = self.points4matplot
			ax.plot(xs, ys, 'x--', lw=2, color=ccolor, ms=10)
		return ax


class PlaneLine:
	class DegRange(Enum):
		ZERO__PI = auto()
		NEG_PI__POS_PI = auto()
		ZERO__2PI = auto()

	def __init__(self, points):
		if not isinstance(points, list):
			msg = "We assume 'list' object but get {0}".format(type(points))
			raise TypeError(msg)

		if all(isinstance(e, Point) for e in points):
			pass
		else:
			points = [Point(e) for e in points]

		if not len(points) == 2:
			msg = "You must be input 2 points"
			raise ValueError(msg)

		self._points = points

	@property
	def plist(self):
		return [e.point for e in self._points]

	@property
	def points(self):
		return self._points

	@property
	def points4matplot(self):
		xs, ys = zip(*self.plist)
		return xs, ys

	@property
	def equation_coefficient(self):
		a = self._points[1].y - self._points[0].y
		b = self._points[0].x - self._points[1].x
		c = self._points[1].x * self._points[0].y - \
			self._points[0].x * self._points[1].y
		return a, b, c

	@property
	def midpoint(self):
		return Point((self.plist[0] + self.plist[1]) / 2.0)

	@property
	def length(self):
		return self._points[0].distance(self._points[1])


	def intersection(self, pl, error=1e-12):
		"""
		Fetch the intersection point of two lines.

		Parameters
		----------
		pl : PlaneLine object
			Target line.
		error : float
			Small values for avoiding calculation errors.

		Returns
		-------
		p : Point object
			Intersection coodinates.
		
		Note 
		-------
		If no intersection, return None.
		"""
		if not isinstance(pl, PlaneLine):
			msg = "Invalid arguments"
			raise TypeError(msg)
		p1, p2 = self._points
		p3, p4 = pl.points

		Xp = ((p3.y * p4.x - p3.x * p4.y) * (p2.x - p1.x) - (p1.y * p2.x - p1.x * p2.y) *
		      (p4.x - p3.x)) / ((p2.y - p1.y) * (p4.x - p3.x) - (p2.x - p1.x) * (p4.y - p3.y) + error)
		Yp = ((p3.y * p4.x - p3.x * p4.y) * (p2.y - p1.y) - (p1.y * p2.x - p1.x * p2.y) *
		      (p4.y - p3.y)) / ((p2.y - p1.y) * (p4.x - p3.x) - (p2.x - p1.x) * (p4.y - p3.y) + error)

		res = (self.is_point_on_line((Xp, Yp)) and pl.is_point_on_line((Xp, Yp)))
		if res:
			return Point(Xp, Yp)
		return None


	def is_point_on_line(self, point, error=1e-3):
		"""
		Check the point is on the line.

		Parameters
		----------
		point : Point object or number tuple or number list
			The value of point to check.
		error : float
			Calculation precision.

		Returns
		-------
		result : bool
			Returns True if the point is on the line, False otherwise.
		
		Note 
		-------
		Using triangle inequality algorithm.
		"""
		p1 = self._points[0]
		p2 = self._points[1]
		if isinstance(point, Point):
			p3 = point
		else:
			p3 = Point(point)

		## Triangle inequality
		ac = p1.distance(p3)
		bc = p3.distance(p2)
		ab = p1.distance(p2)

		return ac + bc < ab + error


	def translation(self, move_dist=None):
		"""
		Create a new PlaneLine object by translating this PlaneLine.

		Parameters
		----------
		move_dist : Point object or ndarray
			Amount to translate.
		
		Returns
		-------
		pl : PlaneLine object
			PlaneLine after translation.
		
		Note 
		-------
		if 'move_dist' is none, move to based on (0,0).
		"""
		if move_dist is None:
			## if none, set first point to move (0,0)
			move_dist = -self.plist[0]
		if isinstance(move_dist, Point):
			move_dist = move_dist.point
		plst = self.plist
		return PlaneLine([p + move_dist for p in plst])


	def radian(self, line, rad_range=None):
		"""
		Calculate the angle between two straight lines.

		Parameters
		----------
		line : PlaneLine object
			Target line.
		rad_range : DegRange
			Range of radian.

		Returns
		-------
		radian : float
			Value of the radian.
		"""
		if rad_range is None:
			rad_range = self.DegRange.ZERO__PI
		base_line = self.translation()
		line = line.translation()

		a = np.dot(base_line.plist[1], line.plist[1])
		b = np.cross(base_line.plist[1], line.plist[1])
		cos = a / (base_line.length * line.length + 1e-12)
		theta = np.arccos(cos)

		if rad_range == self.DegRange.ZERO__PI:
			pass
		elif rad_range == self.DegRange.NEG_PI__POS_PI:
			if b < 0:
				theta = -theta
		elif rad_range == self.DegRange.ZERO__2PI:
			if b < 0:
				theta = 2 * np.pi - theta
		return theta


	def plot(self, ax=None, linestyle="o-", color=None):
		if ax is None:
			_, ax = plt.subplots()

		xs, ys = self.points4matplot
		ax.plot(xs, ys, linestyle, color=color)
		return ax


def main():
	b1 = Bezier([(0, -1), (1.0/3.0, 3), (2.0/3.0, -4), (1.0, 3)])
	base = PlaneLine([(0, 0), (1, 0)])
	ax = b1.plot()
	#base.plot(ax)

	nb1, nb2 = b1.split(0.5)
	nb1.plot(ax, bcolor = "red")
	nb2.plot(ax, bcolor = "blue")

	for p in nb2.points:
		p.plot(ax)

	plt.grid()
	#ax.set_aspect('equal')
	plt.show()


if __name__ == "__main__":
	main()
