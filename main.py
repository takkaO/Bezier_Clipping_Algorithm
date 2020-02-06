import numpy as np
import matplotlib.pyplot as plt
from line_module import Point, Bezier, PlaneLine

class BezierClipping:
	def _check_arguments(self, target_bezier, target_line):
		if not isinstance(target_bezier, Bezier):
			return False
		if not isinstance(target_line, PlaneLine):
			return False
		return True

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
		
		Note 
		-------
		The standard for the distance between a Bezier curve and a straight line is a point on the Bezier curve.
		"""
		if not self._check_arguments(target_bezier, target_line):
			msg = "Invalid arguments."
			raise TypeError(msg.format('target_line', type(target_line)))

		error = 1e-12
		a, b, c = target_line.equation_coefficient
		control_point = []
		n = target_bezier.dims
		for i, point in enumerate(target_bezier.points):
			d = (a * point.x + b * point.y + c) / (np.sqrt(a**2 + b**2) + error)
			control_point.append((i/n, d))
		return Bezier(control_point)
	
	def clipping(self, bezier, line, precision = 1e-3):
		_bezier = self.convert_to_distance_based_bezier(bezier, line)
		_line = PlaneLine([(0,0), (1,0)])
		dic = {}
		dic["base_bezier"] = _bezier
		dic["current_bezier"] = _bezier
		dic["base_line"] = _line
		dic["current_line"] = _line
		dic["precision"] = precision

		self._clipping(dic)
	
	def _clipping(self, dic):
		convexhull = dic["current_bezier"].getConvexhullLines()
		t = []
		for line in convexhull:
			p = line.cross_point(dic["current_line"])
			if p is None:
				continue
			t.append(p)
		
		print(t)

def main():
	b1 = Bezier([(0, -1), (1.0/3.0, 3), (2.0/3.0, -4), (1.0, 3)])
	l1 = PlaneLine([(0,0), (1,0)])
	ax = b1.plot()
	l1.plot(ax)

	bc = BezierClipping()
	bc.clipping(b1, l1)

	plt.grid()
	plt.show()

if __name__ == "__main__":
	main()
