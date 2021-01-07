import sympy as sym
import numpy as np

from metric import Metric
from christoffelSymbols import ChristoffelSymbol2ndKind

class RiemannTensor(object):
	"""
	Represents the Riemann tensor R^mu_nu,rho,sigma.
    See https://www.visus.uni-stuttgart.de/publikationen/catalogue-of-spacetimes

	Parameters
	----------
	g : Metric
		The metric of which the Christoffel symbols should be calculated.
	G : ChristoffelSymbol2ndKind
		The class representing the Christoffel symbols of the second kind for the metric specified in g.
	symbols : numpy array or array_like of sympy.core.symbol.Symbol
		The symbols representing the different dimensions, e.g. `sympy.symbols('r theta z', real=True)`.

	Returns
	-------
	RiemannTensor
		the Riemann tensor
	"""
	def __init__(self, g: Metric, G: ChristoffelSymbol2ndKind, x):
		self.g = g
		self.G = G
		self.x = x
		self.dim = g.dim()
		self.uddd_values = np.zeros((self.dim, self.dim, self.dim, self.dim), dtype='O')
		self.dddd_values = np.zeros((self.dim, self.dim, self.dim, self.dim), dtype='O')
		self.uddd_evaluated = np.full((self.dim, self.dim, self.dim, self.dim), False)
		self.dddd_evaluated = np.full((self.dim, self.dim, self.dim, self.dim), False)

	def uddd(self, mu: int, nu: int, rho: int, sigma: int):
		"""
		Calculates the component of the Riemann tensor R^mu_nu,rho,sigma.

		Parameters
		----------
		mu : int
			The upper index
		nu : int
			The first lower index
		rho : int
			The second lower index
		sigma : int
			The third lower index

		Returns
		-------
		sympy expression
			the expression representing the specified component of the Riemann tensor
		"""
		if self.uddd_evaluated[mu,nu,rho,sigma]: return self.uddd_values[mu,nu,rho,sigma]
		R = self.G.udd(mu,nu,sigma).diff(self.x[rho]) - self.G.udd(mu,nu,rho).diff(self.x[sigma])
		for lam in range(self.g.dim()):
			R += self.G.udd(mu,rho,lam)*self.G.udd(lam,nu,sigma) - self.G.udd(mu,sigma,lam)*self.G.udd(lam,nu,rho)
		self.uddd_values[mu,nu,rho,sigma] = R.simplify()
		self.uddd_values[mu,nu,rho,sigma] = sym.simplify(sym.expand_trig(self.uddd_values[mu,nu,rho,sigma]))
		self.uddd_evaluated[mu,nu,rho,sigma] = True
		return self.uddd_values[mu,nu,rho,sigma]

	def dddd(self, mu: int, nu: int, rho: int, sigma: int):
		"""
		Calculates the component of the Riemann tensor R_mu,nu,rho,sigma.

		Parameters
		----------
		mu : int
			The first lower index
		nu : int
			The second lower index
		rho : int
			The third lower index
		sigma : int
			The fourth lower index

		Returns
		-------
		sympy expression
			the expression representing the specified component of the Riemann tensor
		"""
		if self.dddd_evaluated[mu,nu,rho,sigma]: return self.dddd_values[mu,nu,rho,sigma]
		R = 0
		for lam in range(self.g.dim()):
			R += self.g.dd(mu,lam)*self.uddd(lam,nu,rho,sigma)
		self.dddd_values[mu,nu,rho,sigma] = sym.simplify(sym.expand_trig(R.simplify()))
		self.dddd_evaluated[mu,nu,rho,sigma] = True
		return self.dddd_values[mu,nu,rho,sigma]