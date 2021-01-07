import sympy as sym

from metric import Metric
from coordinate_system_implementation_generator import JavaCoordinateSystemCreator

sigma = sym.symbols('sigma', real=True)
tau   = sym.symbols('tau',   real=True)
phi   = sym.symbols('phi',   real=True)

a = sym.symbols('a', real=True, constant=True)
extra_parms = {"a": "Constant"}

factor = a/(sym.cosh(tau)-sym.cos(sigma))
m = Metric.fromTransformation([sigma, tau, phi], to_base_point=[factor*sym.sin(sigma)*sym.cos(phi),factor*sym.sin(sigma)*sym.sin(phi),factor*sym.sinh(tau)])

axes = {0: "SiDerivedUnit.RADIAN", 1: "SiDerivedUnit.RADIAN", 2: "SiDerivedUnit.RADIAN"}

code = JavaCoordinateSystemCreator("BisphericalCoordinates", m, axes, extra_parameters=extra_parms)
code.coordinate_system_symbols = ["bispherical", "bisph"]

x = sym.symbols('x', real=True)
y = sym.symbols('y', real=True)
z = sym.symbols('z', real=True)

code.base_symbols = [x, y, z]

R = sym.sqrt(x**2+y**2+z**2)
Q = sym.sqrt((R**2+a**2)**2 - (2*a*z)**2)

code.from_base_point = [sym.acos((R**2-a**2)/Q),sym.asinh(2*a*z/Q),sym.atan2(y,x)]

code.save_to("../java")