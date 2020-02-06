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

	def convert_to_nonparametric_bezier(self, target_bezier, target_line):	
		if not self._check_arguments(target_bezier, target_line):
			msg = "Invalid arguments."
			raise TypeError(msg.format('target_line', type(target_line)))

		error = 1e-12
		a, b, c = target_line.equation_coefficient
		print(a,b,c)
		control_point = []
		n = target_bezier.dims
		for i, point in enumerate(target_bezier.points):
			d = (a * point.x + b * point.y + c) / (np.sqrt(a**2 + b**2) + error)
			print(d)
			control_point.append((i/n, d))
		return Bezier(control_point)

def main():
	b1 = Bezier([(0, -1), (1.0/3.0, 3), (2.0/3.0, -4), (1.0, 3)])
	l1 = PlaneLine([(0,0), (1,0)])
	ax = b1.plot()
	l1.plot(ax)

	bc = BezierClipping()
	npb = bc.convert_to_nonparametric_bezier(b1, l1)

	npb.plot(ax, bcolor="blue")	
	plt.grid()
	plt.show()

if __name__ == "__main__":
	main()
