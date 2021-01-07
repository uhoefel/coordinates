import sympy as sym

from metric import Metric
from coordinate_system_implementation_generator import JavaCoordinateSystemCreator

sigma = sym.symbols('sigma', real=True)
tau   = sym.symbols('tau',   real=True)

a = sym.symbols('a', real=True, constant=True)
extra_parms = {"a": "Constant"}

factor = a/(sym.cosh(tau)-sym.cos(sigma))
m = Metric.fromTransformation([sigma, tau], to_base_point=[factor*sym.sinh(tau),factor*sym.sin(sigma)])

axes = {0: "SiDerivedUnit.RADIAN", 1: "SiDerivedUnit.RADIAN"}

code = JavaCoordinateSystemCreator("BipolarCoordinates", m, axes, extra_parameters=extra_parms)
code.coordinate_system_symbols = ["bipolar", "bipol"]

x = sym.symbols('x', real=True)
y = sym.symbols('y', real=True)

code.base_symbols = [x, y]

code.from_base_point = [sym.atan2(2*a*y, x**2+y**2-a**2),0.5*sym.log(((x+a)**2+y**2)/((x-a)**2+y**2))]

code.variable_lower_limit['phi'] = '-Math.PI'
code.variable_upper_limit['phi'] =  'Math.PI'

code.save_to("../java")