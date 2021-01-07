import sympy as sym

from metric import Metric
from coordinate_system_implementation_generator import JavaCoordinateSystemCreator

sigma = sym.symbols('sigma',real=True, positive=True)
tau   = sym.symbols('tau',  real=True, positive=True)
phi   = sym.symbols('phi',  real=True, positive=True)

a = sym.symbols('a', real=True, positive=True, constant=True)
extra_parms = {"a": "Constant"}

trafo = []
trafo.append(a*sym.sinh(sigma)*sym.sin(tau)*sym.cos(phi))
trafo.append(a*sym.sinh(sigma)*sym.sin(tau)*sym.sin(phi))
trafo.append(a*sym.cosh(sigma)*sym.cos(tau))

m = Metric.fromTransformation([sigma, tau, phi], to_base_point=trafo)
#m = Metric.fromTransformation([sigma, tau, phi], to_base_point=[a*sym.sqrt((sigma**2-1)*(1-tau**2))*sym.cos(phi),a*sym.sqrt((sigma**2-1)*(1-tau**2))*sym.sin(phi),a*sigma*tau])

axes = {0: "SiDerivedUnit.RADIAN", 1: "SiDerivedUnit.RADIAN", 2: "SiDerivedUnit.RADIAN"}

code = JavaCoordinateSystemCreator("ProlateSpheroidalCoordinates", m, axes, extra_parameters=extra_parms)
code.coordinate_system_symbols = ["prolate_spheroidal"]

x = sym.symbols('x', real=True)
y = sym.symbols('y', real=True)
z = sym.symbols('z', real=True)

code.base_symbols = [x, y, z]

r = sym.sqrt(x**2+y**2)
first_foci_distance  = sym.sqrt(r**2+(z+a)**2)
second_foci_distance = sym.sqrt(r**2+(z-a)**2)

inv_trafo = []
inv_trafo.append(sym.acosh((first_foci_distance+second_foci_distance)/(2*a)))
inv_trafo.append(sym.sign(z)*sym.acos((first_foci_distance-second_foci_distance)/(2*a)))
inv_trafo.append(sym.atan2(y,x))

code.from_base_point = inv_trafo
#code.from_base_point = [(first_foci_distance+second_foci_distance)/(2*a), (first_foci_distance-second_foci_distance)/(2*a),sym.atan2(y,x)]

code.variable_lower_limit['mu'] = 0

code.variable_lower_limit['nu'] = 0
code.variable_upper_limit['nu'] = 'Math.PI'

code.variable_lower_limit['phi'] = 0
code.variable_upper_limit['phi'] = '2*Math.PI'

code.save_to("../java")