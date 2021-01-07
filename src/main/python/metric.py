import numpy as np
import sympy as sym

from sympy.matrices import dotprodsimp 

class Metric(object):
	"""
	Represents a metric.
    See https://www.visus.uni-stuttgart.de/publikationen/catalogue-of-spacetimes

	Parameters
	----------
	symbols : numpy array or array_like of sympy.core.symbol.Symbol
		The symbols representing the different dimensions, e.g. `sympy.symbols('r theta z', real=True)`. 
		Note that it is good to provide as much information on the symbols as possible, e.g. in the 
		aforementioned example that 'r' is positive
	m : Matrix
		The (covariant) metric represented by a Matrix
	to_base_point : numpy array or array_like, optional
		The transformation equation to get from the current coordinate system into the base coordinate system. 
		Should be expressed using the `symbols` and sympy functions for, e.g., the cosine.
		Note that not providing this may limit the scope on subsequent operations with this metric.

	Returns
	-------
	Metric
		the (analytical) metric as derived from the transformation for a point into the base coordinate system
	"""
	def __init__(self,symbols,m: sym.Matrix,to_base_point=None):
		self.gdd = m  # covariant
		print("Inverting metric...")

		# necessary due to a regression in sympy
		with dotprodsimp(False):
			self.guu = m.inv()  # contravariant

		print("Metric successfully inverted")
		self.symbols = symbols
		self.to_base_point = to_base_point

	@staticmethod
	def fromTransformation(symbols,to_base_point):
		"""
		Creates a new metric from the given to_base_point transformation.

		Parameters
		----------
		symbols : numpy array or array_like of sympy.core.symbol.Symbol
			The symbols representing the different dimensions, e.g. `sympy.symbols('r theta z', real=True)`
		to_base_point : numpy array or array_like
			The transformation equation to get from the current coordinate system into the base coordinate system. 
			Should be expressed using the `symbols` and sympy functions for, e.g., the cosine.

		Returns
		-------
		Metric
			the (analytical) metric as derived from the transformation for a point into the base coordinate system
		"""
		if to_base_point is None:
			raise Exception("You need to specify the point-to-base-transformation!")

		print("Reconstructing metric from point-to-base-transformation...")
		dim = len(to_base_point)
		partial_derivatives = [[0 for i in range(dim)] for i in range(dim)]
		
		for i in range(dim):
			for j in range(dim):
				partial_derivatives[i][j] = sym.diff(to_base_point[j], symbols[i])
		
		g = [[0 for i in range(dim)] for i in range(dim)]
		for i in range(dim):
			for j in range(dim):
				for k in range(dim):
					g[i][j] += partial_derivatives[i][k] * partial_derivatives[j][k]
				g[j][i] = g[i][j] = sym.simplify(sym.expand_trig(g[i][j].simplify()))  # catch some zeros

		print("Metric successfully reconstructed!")

		return Metric(symbols,sym.Matrix(g),to_base_point)			
	
	def __str__(self):
		"""
		Creates a representative string of the (covariant) metric.

		Returns
		-------
		str
			the analytic metric wth covariant components represented by a string
		"""
		return "g_dd = " + str(self.gdd)
	
	def dd(self,i,j):
		"""
		Gets the covariant metric coefficient at the specified position.

		Parameters
		----------
		i : int
			The row index in the metric
		j : int
			The column index in the metric

		Returns
		-------
		sympy expression
			the (analytical) covariant metric coefficient g_ij
		"""
		return self.gdd[i,j]
	
	def uu(self,i,j):
		"""
		Gets the contravariant metric coefficient at the specified position.

		Parameters
		----------
		i : int
			The row index in the metric
		j : int
			The column index in the metric

		Returns
		-------
		sympy expression
			the (analytical) contravariant metric coefficient g^ij
		"""
		return self.guu[i,j]

	def dim(self) -> int:
		"""
		Gets the dimensionality of the metric.
		For example, a metric with size 3x3 will be of dimensionality 3.

		Returns
		-------
		int
			the dimension of the metric
		"""
		return int(np.sqrt(len(self.gdd)))