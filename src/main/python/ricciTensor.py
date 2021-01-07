import sympy as sym
import numpy as np

from riemannTensor import RiemannTensor

class RicciTensor(object):
	"""
	Represents the Ricci tensor R_mu,nu.
    See https://www.visus.uni-stuttgart.de/publikationen/catalogue-of-spacetimes

	Parameters
	----------
	R : RiemannTensor
		The Riemann tensor corresponding to which the Ricci tensor should be calculated.

	Returns
	-------
	RicciTensor
		the Ricci tensor
	"""
	def __init__(self, R: RiemannTensor):
		self.R = R
		self.dim = R.g.dim()
		self.Ric = np.zeros((self.dim, self.dim), dtype='O')
		self.evaluated = np.full((self.dim, self.dim), False)

	def dd(self, mu: int, nu: int):
		"""
		Calculates the specified component of the Ricci tensor R_mu,nu.

		Parameters
		----------
		mu : int
			The first lower index.
		nu : int
			The second lower index.

		Returns
		-------
		sympy expression
			the expression representing the specified component of the Ricci tensor
		"""
		if self.evaluated[mu,nu]: return self.Ric[mu,nu]
		Ric = 0
		for rho in range(self.dim):
			Ric += self.R.uddd(rho,mu,rho,nu)
		self.Ric[mu,nu] = Ric.simplify()
		self.evaluated[mu,nu] = True
		return self.Ric[mu,nu]