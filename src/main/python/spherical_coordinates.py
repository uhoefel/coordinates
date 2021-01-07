import sympy as sym

from metric import Metric
from coordinate_system_implementation_generator import JavaCoordinateSystemCreator

use_4D = False

r      = sym.symbols('r',      real=True, positive=True)
theta1 = sym.symbols('theta1', real=True)
theta2 = sym.symbols('theta2', real=True)
theta3 = sym.symbols('theta3', real=True)

trafo = []

if use_4D:
    symbols = [r, theta1, theta2, theta3]
    trafo.append(r*sym.cos(theta1))
    trafo.append(r*sym.sin(theta1)*sym.cos(theta2))
    trafo.append(r*sym.sin(theta1)*sym.sin(theta2)*sym.cos(theta3))
    trafo.append(r*sym.sin(theta1)*sym.sin(theta2)*sym.sin(theta3))
else:
    # 3D
    symbols = [r, theta1, theta2]
    trafo.append(r*sym.cos(theta1))
    trafo.append(r*sym.sin(theta1)*sym.cos(theta2))
    trafo.append(r*sym.sin(theta1)*sym.sin(theta2))

m = Metric.fromTransformation(symbols, to_base_point=trafo)

axes = {-1: "SiDerivedUnit.RADIAN", 0: "SiBaseUnit.METER"}

code = JavaCoordinateSystemCreator("SphericalCoordinates", m, axes)
code.coordinate_system_symbols = ["spherical", "sph"]

x1 = sym.symbols('x1', real=True)
x2 = sym.symbols('x2', real=True)
x3 = sym.symbols('x3', real=True)
x4 = sym.symbols('x4', real=True)

inv_trafo = []

if use_4D:
    inv_symbols = [x1,x2,x3,x4]
    inv_trafo.append(sym.sqrt(x1**2+x2**2+x3**2+x4**2))
    inv_trafo.append(sym.atan2(sym.sqrt(x2**2+x3**2+x4**2),x1))
    inv_trafo.append(sym.atan2(sym.sqrt(x3**2+x4**2),x2))
    inv_trafo.append(sym.atan2(x4,x3))
else:
    # 3D
    inv_symbols = [x1,x2,x3]
    inv_trafo.append(sym.sqrt(x1**2+x2**2+x3**2))
    inv_trafo.append(sym.atan2(sym.sqrt(x2**2+x3**2),x1))
    inv_trafo.append(sym.atan2(x3,x2))

code.base_symbols = inv_symbols
code.from_base_point = inv_trafo

code.save_to("../java")