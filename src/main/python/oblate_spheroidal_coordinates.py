import sympy as sym

from metric import Metric
from coordinate_system_implementation_generator import JavaCoordinateSystemCreator

mu  = sym.symbols('mu',  real=True, positive=True)
nu  = sym.symbols('nu',  real=True)
phi = sym.symbols('phi', real=True)

a = sym.symbols('a', real=True, positive=True, constant=True)
extra_parms = {"a": "Constant"}

trafo = []
trafo.append(a*sym.cosh(mu)*sym.cos(nu)*sym.cos(phi))
trafo.append(a*sym.cosh(mu)*sym.cos(nu)*sym.sin(phi))
trafo.append(a*sym.sinh(mu)*sym.sin(nu))

m = Metric.fromTransformation([mu, nu, phi], to_base_point=trafo)

axes = {0: "SiDerivedUnit.RADIAN", 1: "SiDerivedUnit.RADIAN", 2: "SiDerivedUnit.RADIAN"}

code = JavaCoordinateSystemCreator("OblateSpheroidalCoordinates", m, axes, extra_parameters=extra_parms)
code.coordinate_system_symbols = ["oblate_spheroidal"]

x = sym.symbols('x', real=True)
y = sym.symbols('y', real=True)
z = sym.symbols('z', real=True)

code.base_symbols = [x, y, z]

r = sym.sqrt(x**2+y**2)
first_foci_distance  = sym.sqrt((r+a)**2+z**2)
second_foci_distance = sym.sqrt((r-a)**2+z**2)

inv_trafo = []
inv_trafo.append(sym.acosh((first_foci_distance+second_foci_distance)/(2*a)))
inv_trafo.append(sym.sign(z)*sym.acos((first_foci_distance-second_foci_distance)/(2*a)))
inv_trafo.append(sym.atan2(y,x))

code.from_base_point = inv_trafo

code.variable_lower_limit['mu'] = 0

code.variable_lower_limit['nu'] = '-Math.PI/2'
code.variable_upper_limit['nu'] =  'Math.PI/2'

code.variable_lower_limit['phi'] = '-Math.PI'
code.variable_upper_limit['phi'] =  'Math.PI'

code.save_to("../java")