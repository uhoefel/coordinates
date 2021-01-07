import sympy as sym

from metric import Metric
from coordinate_system_implementation_generator import JavaCoordinateSystemCreator

u = sym.symbols('u', real=True)
v = sym.symbols('v', real=True)
w = sym.symbols('w', real=True)

denominator = u**2 + v**2 + w**2
m = Metric.fromTransformation([u, v, w], to_base_point=[u/denominator, v/denominator, w/denominator])

axes = {-1: "m^-1"}

code = JavaCoordinateSystemCreator("SixSphereCoordinates", m, axes)
code.coordinate_system_symbols = ["6sphere", "6sph"]

x = sym.symbols('x', real=True)
y = sym.symbols('y', real=True)
z = sym.symbols('z', real=True)

code.base_symbols = [x, y, z]

denominator = x**2 + y**2 + z**2
code.from_base_point = [x/denominator, y/denominator, z/denominator]

code.save_to("../java")