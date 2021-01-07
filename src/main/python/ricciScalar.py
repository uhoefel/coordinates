import sympy as sym

from metric import Metric
from ricciTensor import RicciTensor

class RicciScalar(object):
	"""
	Represents the Ricci scalar.
    See https://www.visus.uni-stuttgart.de/publikationen/catalogue-of-spacetimes

	Parameters
	----------
	R : RicciTensor
		The Ricci tensor corresponding from which the Ricci scalar should be calculated.
	g : Metric
		The metric corresponding to the Ricci tensor

	Returns
	-------
	RicciScalar
		the Ricci scalar
	"""
	def __init__(self, Ric: RicciTensor, g: Metric):
		self.Ric = Ric
		self.g = g

	def value(self):
		"""
		Calculates the Ricci scalar.

		Returns
		-------
		sympy expression
			the expression representing the Ricci scalar
		"""
		RS = 0
		for mu in range(self.g.dim()):
			for nu in range(self.g.dim()):
				RS += self.g.uu(mu,nu)*self.Ric.dd(mu,nu)
		return RS.simplify()