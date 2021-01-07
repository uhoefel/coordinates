import sympy as sym

from metric import Metric
from coordinate_system_implementation_generator import JavaCoordinateSystemCreator

r = sym.symbols('r', real=True, positive=True)
theta = sym.symbols('theta', real=True)
z = sym.symbols('z', real=True)

m = Metric.fromTransformation([r, theta, z], to_base_point=[r*sym.cos(theta),r*sym.sin(theta),z])

axes = {0: "SiBaseUnit.METER", 1: "SiDerivedUnit.RADIAN", 2: "SiBaseUnit.METER"}

code = JavaCoordinateSystemCreator("CylindricalCoordinates", m, axes)

x = sym.symbols('x', real=True)
y = sym.symbols('y', real=True)
z = sym.symbols('z', real=True)

code.base_symbols = [x, y, z]
code.from_base_point = [sym.sqrt(x**2 +y**2), sym.atan2(y, x), z]

code.variable_lower_limit['r'] = 0
code.variable_lower_limit['theta'] = 0
code.variable_upper_limit['theta'] = '2*Math.PI'

code.save_to("../java")
