import sympy as sym

from metric import Metric
from coordinate_system_implementation_generator import JavaCoordinateSystemCreator

x = sym.symbols('x', real=True)
y = sym.symbols('y', real=True)
z = sym.symbols('z', real=True)

m = Metric.fromTransformation([x, y, z], to_base_point=[x, y, z])

axes = {-1: "SiBaseUnit.METER"}

code = JavaCoordinateSystemCreator("CartesianCoordinates", m, axes)
code.coordinate_system_symbols = ["cartesian", "cart"]

code.base_symbols = [x, y, z]
code.from_base_point = [x, y, z]

code.save_to("../java")
