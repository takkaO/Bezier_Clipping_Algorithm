import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
if __name__ == "__main__":
	## For run as main process
	from line_module import Point, Bezier, PlaneLine
else:
	## For call from submodules
	from .line_module import Point, Bezier, PlaneLine


class BezierClipping:
	def _check_arguments(self, target_bezier, target_line):
		if not isinstance(target_bezier, Bezier):
			return False
		if not isinstance(target_line, PlaneLine):
			return False
		return True
	

	def _flatten(self, data):
		if not isinstance(data, list) and not isinstance(data, tuple):
			return [data]
		return [element for item in data for element in (self._flatten(item) if hasattr(item, '__iter__') else [item])]


	def convert_to_distance_based_bezier(self, target_bezier, target_line):
		"""
		Create a non-parametric Bezier curve with a distance from a target Bezier curve to a target line.

		Parameters
		----------
		target_bezier : Bezier object
			Bezier curve object.
		target_line : PlaneLine object
			Straight line object.

		Returns
		-------
		bezier_curve : Bezier object
			New non-parametric Bezier curve.
		"""
		if not self._check_arguments(target_bezier, target_line):
			msg = "Invalid arguments."
			raise TypeError(msg.format('target_line', type(target_line)))

		error = 1e-12
		a, b, c = target_line.equation_coefficient
		control_point = []
		n = target_bezier.dims
		for i, point in enumerate(target_bezier.points):
			d = -(a * point.x + b * point.y + c) / (np.sqrt(a**2 + b**2) + error)
			control_point.append((i/n, d))
		return Bezier(control_point)
	

	def detect_intersection(self, bezier, line, precision=1e-3):
		"""
		Fetches the intersection of a straight line and a Bezier curve.

		Parameters
		----------
		bezier : Bezier object
			Target Bezier curve.
		line : PlaneLine object
			Target line.
		precision : float
			Precision of intersection value.

		Returns
		-------
		points : List of tuple
			[(t, Point), ...]. t is parameter value, Point is (x, y) coodinates.
		"""
		t_values = self.clipping(bezier, line, precision)

		on_line = lambda t: line.is_point_on_line(bezier.bezier_point(t))
		p_lst = [(t, bezier.bezier_point(t)) for t in t_values if on_line(t)]
		
		return p_lst

	def clipping(self, bezier, line, precision = 1e-3):
		"""
		Bezier clipping.

		Parameters
		----------
		bezier : Bezier object
			Target Bezier curve.
		line : PlaneLine object
			Target line.
		precision : float
			Precision of intersection value.

		Returns
		-------
		t_values : list
			Intersection 't' value. Not flatten list.
		
		Note 
		-------
		It is not checked whether the intersection is on the line segment.
		"""
		_bezier = self.convert_to_distance_based_bezier(bezier, line)
		_line = PlaneLine([(0,0), (1,0)])
		dic = {}
		dic["base_bezier"] = _bezier
		dic["current_bezier"] = _bezier
		dic["base_line"] = _line
		dic["current_line"] = _line
		dic["precision"] = precision

		ret = self._clipping(dic)
		return self._flatten(ret)
	

	def _clipping(self, dic):
		convexhull = dic["current_bezier"].getConvexhullLines()
		t = []
		for line in convexhull:
			p = line.intersection(dic["current_line"])
			if p is None:
				continue
			t.append(p)
		
		if len(t) == 0:
			## No intersection
			return []
		elif len(t) == 1:
			return dic["current_line"].midpoint.x
		else:
			t = sorted(t, key=lambda p: p.x)
		
		## Fetch t_min and t_max
		t_min = t[0].x
		if t_min < 0:
			t_min = 0.0
		t_max = t[-1].x
		if t_max < 0.0:
			t_max = 0.0

		next_line = PlaneLine([(t_min, 0), (t_max, 0)])
		if next_line.length <= dic["precision"]:
			return next_line.midpoint.x

		bez, _ = dic["base_bezier"].split(t_max)
		bez_t = t_min / (t_max + 1e-12)
		_, next_bezier = bez.split(bez_t)

		if np.abs(dic["current_line"].length - next_line.length) < 1e-6:
			## Target has more than 2 intersection
			b1, b2 = next_bezier.split(0.5)
			dic1 = deepcopy(dic)
			dic1["current_bezier"] = b1
			dic1["current_line"] = next_line
			dic2 = deepcopy(dic)
			dic2["current_bezier"] = b2
			dic2["current_line"] = next_line
			return self._clipping(dic1), self._clipping(dic2)

		dic["current_bezier"] = next_bezier
		dic["current_line"] = next_line
		return self._clipping(dic)


def main():
	bc = BezierClipping()
	
	## Example1
	b1 = Bezier([(0, -5), (1.0/3.0, 8), (2.0/3.0, 1), (1.2, -6), (1.5, 5)])
	l1 = PlaneLine([(0,-2), (1.5,1.5)])
	ax = b1.plot()
	l1.plot(ax)
	
	res = bc.detect_intersection(b1, l1)	
	print("Result1:")
	if not res == []:
		ts, ps = zip(*res)
		for t, p in zip(ts, ps):
			print("t:{:.3f} -> (x, y) = {}".format(t, p.point))
			p.plot(ax, fmt="or")
	plt.grid()
	
	## Example2
	b2 = Bezier([(3, -5), (5, 8), (8, -1)])
	l2 = PlaneLine([(2, -5), (9, 4)])
	ax = b2.plot()
	l2.plot(ax)

	res = bc.detect_intersection(b2, l2)
	print("Result2:")
	if not res == []:
		ts, ps = zip(*res)
		for t, p in zip(ts, ps):
			print("t:{:.3f} -> (x, y) = {}".format(t, p.point))
			p.plot(ax, fmt="or")

	plt.grid()
	plt.show()

if __name__ == "__main__":
	main()
