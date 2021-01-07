import sympy as sym

from metric import Metric
from coordinate_system_implementation_generator import JavaCoordinateSystemCreator

r = sym.symbols('r', real=True, positive=True)
theta = sym.symbols('theta', real=True)

m = Metric.fromTransformation([r, theta], to_base_point=[r*sym.cos(theta),r*sym.sin(theta)])

axes = {0: "SiBaseUnit.METER", 1: "SiDerivedUnit.RADIAN"}

code = JavaCoordinateSystemCreator("PolarCoordinates", m, axes)
code.coordinate_system_symbols = ["polar", "pol"]

x = sym.symbols('x', real=True)
y = sym.symbols('y', real=True)

code.base_symbols = [x, y]
code.from_base_point = [sym.sqrt(x**2 +y**2), sym.atan2(y, x)]

code.variable_lower_limit['r'] = 0
code.variable_lower_limit['theta'] = 0
code.variable_upper_limit['theta'] = '2*Math.PI'

code.save_to("../java")