import sympy as sym

from metric import Metric
from coordinate_system_implementation_generator import JavaCoordinateSystemCreator

xi  = sym.symbols('xi',  real=True, positive=True)
eta = sym.symbols('eta', real=True)
phi = sym.symbols('phi', real=True, positive=True)

R = sym.symbols('R', real=True, positive=True, constant=True)
extra_parms = {"R": "Constant"}

factor = R / (sym.cosh(xi) - sym.cos(eta))

trafo = []
trafo.append(factor*sym.sinh(xi)*sym.cos(phi))
trafo.append(factor*sym.sinh(xi)*sym.sin(phi))
trafo.append(factor*sym.sin(eta))

m = Metric.fromTransformation([xi, eta, phi], to_base_point=trafo)

axes = {0: "SiDerivedUnit.RADIAN", 1: "SiDerivedUnit.RADIAN", 2: "SiDerivedUnit.RADIAN"}

code = JavaCoordinateSystemCreator("ToroidalCoordinates", m, axes, extra_parameters=extra_parms)
code.coordinate_system_symbols = ["toroidal", "tor"]
code.variable_lower_limit['R'] = 0

x = sym.symbols('x', real=True)
y = sym.symbols('y', real=True)
z = sym.symbols('z', real=True)

code.base_symbols = [x, y, z]

r = sym.sqrt(x**2+y**2+z**2)
rho = sym.sqrt(x**2+y**2)
first_foci_distance = sym.sqrt((rho+R)**2 + z**2)
second_foci_distance = sym.sqrt((rho-R)**2 + z**2)

inv_trafo = []
inv_trafo.append(sym.log(first_foci_distance/second_foci_distance))
inv_trafo.append(sym.sign(x)*sym.acos((-4*R**2+first_foci_distance**2+second_foci_distance**2)/(2*first_foci_distance*second_foci_distance)))
inv_trafo.append(sym.atan2(y,x))

code.from_base_point = inv_trafo

code.variable_lower_limit['eta'] = '-Math.PI'
code.variable_upper_limit['eta'] = 'Math.PI'

code.variable_lower_limit['xi'] = 0

code.variable_lower_limit['phi'] = 0
code.variable_upper_limit['phi'] = '2*Math.PI'

code.save_to("../java")