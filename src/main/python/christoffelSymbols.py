import sympy as sym
import numpy as np

from metric import Metric

class ChristoffelSymbol2ndKind(object):
	"""
	Represents the Christoffel symbols of the second kind, Γ^m_ij, of a metric.
    See https://www.visus.uni-stuttgart.de/publikationen/catalogue-of-spacetimes

	Parameters
	----------
	g : Metric
		The metric of which the Christoffel symbols should be calculated
	symbols : numpy array or array_like of sympy.core.symbol.Symbol
		The symbols representing the different dimensions, e.g. `sympy.symbols('r theta z', real=True)`.

	Returns
	-------
	ChristoffelSymbol2ndKind
		the class representing the Christoffel symbols
	"""
	def __init__(self,g: Metric,x):
		self.g = g
		self.x = x
		self.dim = g.dim()
		self.christoffels = np.zeros((self.dim, self.dim, self.dim), dtype='O')
		self.evaluated = np.full((self.dim, self.dim, self.dim), False)

	def udd(self, m: int, i: int, j: int):
		"""
		Calculates the Christoffel symbols of the second kind, Γ^m_ij.

		Parameters
		----------
		m : int
			the upper index
		i : int
			the first lower index
		j : int
			the second lower index

		Returns
		-------
		sympy expression
			the expression representing the specified Christoffel symbol
		"""
		if self.evaluated[m,i,j]: return self.christoffels[m,i,j]
		christoffel = 0
		for k in range(self.dim):
			christoffel += self.g.uu(m,k)/2 * (self.g.dd(k,i).diff(self.x[j]) + self.g.dd(k,j).diff(self.x[i]) - self.g.dd(i,j).diff(self.x[k]))
		self.christoffels[m,i,j] = christoffel.simplify()
		self.evaluated[m,i,j] = True
		return self.christoffels[m,i,j]

