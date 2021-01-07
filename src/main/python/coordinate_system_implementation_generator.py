import re
import os
import sys
import random

import numpy as np
import sympy as sym

from collections import OrderedDict
from metric import Metric
from christoffelSymbols import ChristoffelSymbol2ndKind
from riemannTensor import RiemannTensor
from ricciTensor import RicciTensor
from ricciScalar import RicciScalar

class JavaCoordinateSystemCreator(object):
	"""
	Class for creating an implementation of the eu.hoefel.coordinates.CoordinateSystem Java interface via sympy.

	Parameters
	----------
	java_class_name : str
		The class name to be used for the Java class.
	metric : Metric
		The metric of the coordinate system.
	axes : dict
		The dictionary containing the dimension as key and the corresponding unit as a string. Note that 
		the key of the default dimension should be -1 and specific implementations like 'SiBaseUnit.METER' 
		are also allowed.
	java_name : str, optional
		The name of the coordinate system. By default, if the java_class_name ends in 'Coordinates', it will 
		use the preceding characters.
	extra_parameters : dict, optional
		Additional parameters required for the coordinate system. Their name should be supplied as key 
		and their type as value.

	Returns
	-------
	JavaCoordinateSystemCreator
		the class allowing to create a Java implementation of the specified coordinate system
	"""
	def __init__(self, java_class_name: str, metric: Metric, axes: dict, java_name=None, extra_parameters=None):
		self.java_class_name = str(java_class_name)
		self.g = metric
		self.dim = self.g.dim()
		
		self.symbols = self.g.symbols

		if self.g.to_base_point is not None:
			self.to_base_point = self.g.to_base_point
			self.from_base_point = None
			self.base_symbols = None
		else:
			self.to_base_point = None
			self.from_base_point = None
			self.base_symbols = None

		self.gamma = ChristoffelSymbol2ndKind(self.g, self.g.symbols)
		self.riem = RiemannTensor(self.g, self.gamma, self.g.symbols)
		self.ric = RicciTensor(self.riem)
		self.ric_scalar = RicciScalar(self.ric,self.g)

		self.code = ""
		
		if java_name == '' or java_name is None:
			index = java_class_name.index('Coordinates')
			if index == -1:
				raise Exception("Cannot guess java_name, please supply explicitly")
			self.java_name = java_class_name[:index]
		else:
			self.java_name = java_name
		
		if axes is None:
			raise Exception("Cannot create default axes!")
		else:
			self.axes = axes

		self.coordinate_system_symbols = []
		
		# needs to be a dictionary with type info
		self.extra_parameters = extra_parameters
		
		self.extra_imports = []
		self.__extra_symbol_check = {}
		
		self.unable_to_invert = False
		self.inv_transformation = None
		
		self.variable_lower_limit = {}
		self.variable_upper_limit = {}
		
		self.__use_validate_position = False
	
	def save_to(self, save_folder: str, package='eu.hoefel.coordinates'):
		"""
		Saves the generated Java implementation to the specified folder.
		
		Parameters
		----------
		save_folder : str
			The path to save to, except that the package gets appended.
		package : str, optional
			The package to use. Defaults to 'eu.hoefel.coordinates'.
		"""
		print("Trying to simplify input.")
		if self.to_base_point is not None:
			for i in range(len(self.to_base_point)):
				self.to_base_point[i] = self.to_base_point[i].simplify()
		if self.from_base_point is not None:
			for i in range(len(self.from_base_point)):
				self.from_base_point[i] = self.from_base_point[i].simplify()
		
		print("Creating implementation. This may take a while, please be patient.")
		print("Please note in the meantime that the resulting file may contain TODOs.")
		print("Most, but potentially not all of them will be related to documentation.")
		print("Please address the TODOs.")
		# run once to collect all things to be imported
		self.__build(True)
		
		# throw away the generated code, but not the collected extra imports
		self.code = ""
		
		# throw away duplicate imports
		self.extra_imports = list(dict.fromkeys(self.extra_imports))

		# Rebuild with the extra imports in place
		self.__build(False)

		if not save_folder.endswith("/"): save_folder += "/"
		save_folder += package.replace(".", "/")
		if not save_folder.endswith("/"): save_folder += "/"

		print("Saving to " + save_folder)
		if not os.path.exists(save_folder):
			os.makedirs(save_folder)
		filename = self.java_class_name + ".java"
		
		self.code = "package " + package + ";\n\n" + self.code
		
		with open(save_folder + filename, "w+") as file: 
			file.write(self.code)

	def __build(self, mute):
		"""
		Builds the code for a Java class implementing the eu.hoefel.coordinates.CoordinateSystem interface.

		Parameters
		----------
		mute : bool
			Print (most) stuff only if False
		"""
		self.code = ""
		
		if not mute: print("Build header...")
		self.__build_header()

		if not mute: print("Build Constructor...")
		self.__build_constructor()

		if not mute: print("Build update position...")
		self.__build_update_position()

		if not mute: print("Build dimension...")
		self.__build_dimension()

		if not mute: print("Build isOrthogonal...")
		self.__build_isOrthogonal()

		if not mute: print("Build toBaseUnits...")
		self.__build_to_base_transformation_units()

		if not mute: print("Build toBasePoint...")
		self.__build_to_base_transformation_point()

		if not mute: print("Build fromBasePoint...")
		self.__build_from_base_transformation_point()

		if not mute: print("Build toBaseVector...")
		self.__build_to_base_transformation_vector()

		if not mute: print("Build fromBaseVector...")
		self.__build_from_base_transformation_vector()

		if not mute: print("Build metric coefficient...")
		self.__build_metric_coefficient()

		if not mute: print("Build metric tensor...")
		self.__build_metric_tensor()

		if not mute: print("Build Jacobian...")
		self.__build_jacobian_determinant()

		if not mute: print("Build Christoffel Symbol 1st kind...")
		self.__build_christoffel_symbols_1st_kind()

		if not mute: print("Build Christoffel Symbol 2nd kind...")
		self.__build_christoffel_symbols_2nd_kind()

		if not mute: print("Build Riemann tensor...")
		self.__build_riemann_tensor(mute=mute)

		if not mute: print("Build isFlat...")
		self.__build_isFlat()

		if not mute: print("Build Ricci tensor...")
		self.__build_ricci_tensor()

		if not mute: print("Build Ricci scalar...")
		self.__build_ricci_scalar()
		
		self.code += "}\n"

	def __is_java_import(self, string: str) -> bool:
		"""
		Checks whether the specified string is a java import.
		
		Parameters
		----------
		string : str
			The string to check.

		Returns
		-------
		Boolean
			true if the string is for a Java import
		"""
		return string.startswith("import java.")

	def __is_eu_hoefel_import(self, string: str) -> bool:
		"""
		Checks whether the specified string is a eu.hoefel.* import.
		
		Parameters
		----------
		string : str
			The string to check.

		Returns
		-------
		Boolean
			true if the string is for a eu.hoefel.* import
		"""
		return string.startswith("import eu.hoefel.")

	def __build_header(self):
		"""
		Builds the header of the Java implementation of the eu.hoefel.coordinates.CoordinateSystem interface 
		and attaches it to the variable representing the Java code.
		"""
		found = False
		for extra_import in self.extra_imports:
			if self.__is_java_import(extra_import):
				self.code += extra_import + "\n"
				found = True
		if found: self.code += "\n"

		found = False
		for extra_import in self.extra_imports:
			if self.__is_eu_hoefel_import(extra_import):
				self.code += extra_import + "\n"
				found = True
		if found: self.code += "\n"

		found = False
		for extra_import in self.extra_imports:
			if not self.__is_java_import(extra_import) and not self.__is_eu_hoefel_import(extra_import):
				self.code += extra_import + "\n"
				found = True
		if found: self.code += "\n"


		self.extra_imports.append("import java.util.NavigableSet;")
		if self.extra_parameters is not None:
			self.extra_imports.append("import java.util.Objects;")
		else:
			self.extra_imports.append("import java.util.function.Consumer;")
		self.extra_imports.append("import eu.hoefel.coordinates.axes.Axes;")
		self.extra_imports.append("import eu.hoefel.coordinates.axes.Axis;")
		self.extra_imports.append("import eu.hoefel.coordinates.CoordinateSystemSymbols;")

		self.code += "/**\n"
		self.code += " * TODO\n"
		self.code += " * \n"
		self.code += " * @param axes the axes defining the coordinate system, not null\n"
		if self.extra_parameters is not None:
			for name, _ in self.extra_parameters.items():
				self.code += " * @param " + name + " TODO\n"
		self.code += " */\n"
		if self.coordinate_system_symbols != []:
			self.code += "@CoordinateSystemSymbols({\"" + "\", \"".join(self.coordinate_system_symbols) + "\"})\n"
		else:
			self.code += "@CoordinateSystemSymbols({\"\"}) // TODO: add symbols representing the coordinate system, e.g. 'cartesian' and 'cart' for the cartesian coordinate system\n"
		self.code += "public final record " + self.java_class_name + "(NavigableSet<Axis> axes"
		if self.extra_parameters is not None:
			for name, param_type in self.extra_parameters.items():
				self.code += ", " + param_type + " " + name
				
		self.code += ") implements CoordinateSystem {\n\n"

		self.code += "	/** The default axes. */\n"
		self.code += "	public static final NavigableSet<Axis> DEFAULT_AXES = Axes.of(\n"
		
		counter = 0
		for dimension, unit in sorted(self.axes.items()):
			axis = ""
			if dimension == -1:  # default axis
				axis += "new Axis("
			else:
				axis += "new Axis(" + str(dimension) + ", "
			if "." in unit:  # probably smth like SiBaseUnit.METER
				axis += unit
				# handle some standard units specially
				if unit.startswith("SiBaseUnit"): self.extra_imports.append("import eu.hoefel.unit.si.SiBaseUnit;")
				if unit.startswith("SiDerivedUnit"): self.extra_imports.append("import eu.hoefel.unit.si.SiDerivedUnit;")
				if unit.startswith("SiCommonUnit"): self.extra_imports.append("import eu.hoefel.unit.si.SiCommonUnit;")
			else:
				axis += "Unit.of(\"" + unit + "\")"
			axis += ")"
			self.code += "			" + axis
			counter += 1
			if counter < len(self.axes): self.code += ",\n"
		self.code += ");\n\n"

		if self.extra_parameters is None:
			self.code += "	/**\n"
			self.code += "	 * The consumer useful for checking the arguments in\n"
			self.code += " 	 * {@link #" + self.java_class_name + "(Object...)}.\n"
			self.code += "	 */\n"
			self.code += "	private static final Consumer<Object[]> ARG_CHECK = args -> Axes.DEFAULT_ARG_CHECK.accept(\"" + self.java_name + "\", args);\n\n"

	def __build_constructor(self):
		"""
		Attaches the Java code for constructors to the code.
		"""
		self.code += "	/** Constructs a new " + self.java_name.lower() + " coordinate system. */\n"
		self.code += "	public " + self.java_class_name + " {\n"
		
		self.extra_imports.append("import java.util.Objects;")

		self.code += "		Objects.requireNonNull(axes, \"Axes may not be null. \"\n"
		self.code += "				+ \"Use the DEFAULT_AXES"
		if self.extra_parameters is not None:
			self.code += " or the constructor that just requires " + str([*self.extra_parameters.keys()]) + " to get a reasonable default."
		else:
			self.code += " or the empty constructor to get a reasonable default."
		self.code += "\");\n"
		
		for dimension in sorted(self.__extra_symbol_check.keys()):
			self.code += "\n"
			self.code += self.__extra_symbol_check[dimension]

		if self.extra_parameters is not None:
			for extra_symbol in sorted(self.extra_parameters.keys()):
				expression_in_equation = extra_symbol
				is_constant = self.extra_parameters[extra_symbol] == 'Constant'
				if is_constant: expression_in_equation+= ".value()"
				
				if extra_symbol in self.variable_lower_limit and extra_symbol in self.variable_upper_limit:
					self.code += "\n"
					self.code += "		if (" + expression_in_equation + " < " + str(self.variable_lower_limit[extra_symbol]) + " || " + expression_in_equation + " > " + str(self.variable_upper_limit[extra_symbol]) + ") {\n"
					self.code += "			throw new IllegalArgumentException(\"'" + extra_symbol + "' needs to be between " + str(self.variable_lower_limit[extra_symbol]) + " and " + str(self.variable_upper_limit[extra_symbol]) + ", but is \" + " + extra_symbol + ");\n"
					self.code += "		}\n"
				elif extra_symbol in self.variable_lower_limit:
					self.code += "\n"
					self.code += "		if (" + expression_in_equation + " < " + str(self.variable_lower_limit[extra_symbol]) + ") {\n"
					self.code += "			throw new IllegalArgumentException(\"'" + extra_symbol + "' needs to be above " + str(self.variable_lower_limit[extra_symbol]) + ", but is \" + " + extra_symbol + ");\n"
					self.code += "		}\n"
				elif extra_symbol in self.variable_upper_limit:
					self.code += "\n"
					self.code += "		if (" + expression_in_equation + " > " + str(self.variable_upper_limit[extra_symbol]) + ") {\n"
					self.code += "			throw new IllegalArgumentException(\"'" + extra_symbol + "' needs to be below " + str(self.variable_upper_limit[extra_symbol]) + ", but is \" + " + extra_symbol + ");\n"
					self.code += "		}\n"
		self.code += "	}\n\n"

		self.code += "	/**\n"
		self.code += "	 * Constructs a new " + self.java_name.lower() + " coordinate system.\n"
		self.code += "	 * \n"
		self.code += "	 * @param args the arguments, must be either {@link Axes} or {@link Axis}, which\n"
		self.code += "	 *			   will take precedence over the {@link #DEFAULT_AXES} if given. If\n"
		self.code += "	 *			   no arguments are given, the default axes will be used.\n"
		if self.extra_parameters is not None:
			self.code += "	 *			   Further required arguments may also be passed in here, but (if \n"
			self.code += "	 *			   they are doubles or {@link Constant}) they have to be in the \n"
			self.code += "	 *			   same order in which they are specified in the record declaration.\n"
			
		self.code += "	 */\n"
		if self.extra_parameters is None:
			self.code += "	public " + self.java_class_name + "(Object... args) {\n"
			self.code += "		this(Axes.fromArgs(DEFAULT_AXES, ARG_CHECK, args));\n"
			self.code += "	}\n\n"
		else:	
			self.code += "	public " + self.java_class_name + "(Object... args) {\n"
			
			add_warning = 'double' in self.extra_parameters.values() and 'Constant' in self.extra_parameters.values()
			if add_warning:
				self.code += "		// Note that using both doubleFromArgs and constantFromArgs is not save if the\n"
				self.code += "		// values are input as strings that only contain a number and not an additional\n"
				self.code += "		// unit (for the constants)\n"

			self.code += "		this(Axes.fromArgs(DEFAULT_AXES, args)"

			doubleIndex = 0
			constantIndex = 0
			for name, param_type in self.extra_parameters.items():
				if param_type == "double":
					self.extra_imports.append("import eu.hoefel.coordinates.CoordinateSystems;")
					self.code += ", CoordinateSystems.doubleFromArgs(" + str(doubleIndex) + ", args)"
					doubleIndex += 1
				elif param_type == 'Constant':
					self.extra_imports.append("import eu.hoefel.unit.constant.Constant;")
					self.code += ", CoordinateSystems.constantFromArgs(" + str(constantIndex) + ", args)"
					constantIndex += 1
				else:
					raise Exception("Don't know how to handle " + param_type + " for variable '" + name + "'!")
			self.code += ");\n"
			self.code += "	}\n\n"
			
			self.code += "	/**\n"
			self.code += "	 * Constructs a new " + self.java_name.lower() + " coordinate system.\n"
			self.code += "	 * \n"
			if self.extra_parameters is not None:
				for name, _ in self.extra_parameters.items():
					self.code += "	 * @param " + name + " TODO\n"
			self.code += "	 */\n"
			self.code += "	public " + self.java_class_name + "("
			
			parms = ""
			for name, param_type in self.extra_parameters.items():
				parms += ", " + param_type + " " + name

			self.code += parms[2:]
			self.code += ") {\n"
			self.code += "		this(DEFAULT_AXES, "

			parms = ""
			for name, param_type in self.extra_parameters.items():
				parms += ", " + name

			self.code += parms[2:]
			self.code += ");\n"
			self.code += "	}\n\n"

	def __build_update_position(self):
		"""
		Attaches the code to validate the position (in the current coordinate system), if limits are given.
		"""
		symbols_str = [str(symbol) for symbol in self.symbols]
		self.__use_validate_position = True  # we probably can always use it

		if self.__use_validate_position:
			self.code += "	/**\n"
			self.code += "	 * Validates the position, i.e. it throws an exception if a dimension of the \n"
			self.code += "	 * position is out of range.\n"
			self.code += "	 * \n"
			self.code += "	 * @param position the position to validate\n"
			self.code += "	 * @throw IllegalArgumentException if the assumptions about the dimensionality \n"
			self.code += "	 *		  or the valid range of any dimension of the input are violated.\n"
			self.code += "	 */\n"
			self.code += "	private void validatePosition(double[] position) {\n"
            self.code += "		Objects.requireNonNull(position);\n"
            self.code += "\n"
			self.code += "		if (position.length > dimension()) {\n"
			self.code += "			throw new IllegalArgumentException(\n"
			self.code += "					\"The given dimensionality exceeds the maximum supported dimensionality (%d vs %d)\"\n"
			self.code += "							.formatted(position.length, dimension()));\n"
			self.code += "		}\n"

			for symbol in symbols_str:
				index = symbols_str.index(symbol)
				if symbol in self.variable_lower_limit and symbol in self.variable_upper_limit:
					self.code += "\n"
					self.code += "		if (position[" + str(index) + "] < " + str(self.variable_lower_limit[symbol]) + " || position[" + str(index) + "] > " + str(self.variable_upper_limit[symbol]) + ") {\n"
					self.code += "			throw new IllegalArgumentException(\"position[" + str(index) + "] needs to be between " + str(self.variable_lower_limit[symbol]) + " and " + str(self.variable_upper_limit[symbol]) + ", but is \" + position[" + str(index) + "]);\n"
					self.code += "		}\n"
				elif symbol in self.variable_lower_limit:
					self.code += "\n"
					self.code += "		if (position[" + str(index) + "] < " + str(self.variable_lower_limit[symbol]) + ") {\n"
					self.code += "			throw new IllegalArgumentException(\"position[" + str(index) + "] needs to be above " + str(self.variable_lower_limit[symbol]) + ", but is \" + position[" + str(index) + "]);\n"
					self.code += "		}\n"
				elif symbol in self.variable_upper_limit:
					self.code += "\n"
					self.code += "		if (position[" + str(index) + "] > " + str(self.variable_upper_limit[symbol]) + ") {\n"
					self.code += "			throw new IllegalArgumentException(\"position[" + str(index) + "] needs to be below " + str(self.variable_upper_limit[symbol]) + ", but is \" + position[" + str(index) + "]);\n"
					self.code += "		}\n"

			self.code += "	}\n\n"

	@staticmethod
	def __simplify(expr):
		# first simplify, if the expression is not too long
		if len(str(expr)) > 75: return expr
		
		# first simplify
		expr = expr.simplify()

		# then expand trigonometric funcs
		expr = sym.expand_trig(expr)

		# then simplify again
		expr = expr.simplify()

		# then try to rewrite it in terms of exp funcs
		expr_exp = expr.rewrite(sym.exp)

		# simplify again
		expr_exp = expr_exp.simplify()

		if len(str(expr_exp)) < len(str(expr)) and not 'I' in str(expr_exp):
			return expr_exp
		return expr
		
		
	def __build_to_base_transformation_units(self):
		"""
		Attaches the Java code for a method that provides the corresponding units in the corresponding base coordinate system.
		"""
		self.extra_imports.append("import java.util.NavigableMap;")
		self.extra_imports.append("import eu.hoefel.unit.Unit;")
		self.extra_imports.append("import eu.hoefel.unit.Units;")
		
		self.code += "	@Override\n"
		self.code += "	public NavigableMap<Integer, Unit> toBaseUnits() {\n"
		self.code += "		NavigableMap<Integer, Unit> map = new TreeMap<>();\n"
		
		self.extra_imports.append("import java.util.TreeMap;")

		#try:
		for dim, transformation in enumerate(self.to_base_point):
			# potentially add check for radian equality of some units
			self.__check_for_func_args(transformation)
			
			# remove specific funcs that do not change the unit
			expr = self.__remove_specific_funcs(transformation)

			# remove funcs (but leave their arguments in place) that do not change the unit
			expr = self.__remove_irrelevant_funcs(expr)

			# now we should have a 'clean' equation, i.e. only +* and ** should be left
			# get rid of + and keep only one term as summands don't change the unit
			expr = self.__keep_one_summand(expr).simplify()

			self.dimension_exponents = {}

			# so we have only * and ** left
			self.__collect_exponents(expr, self.symbols)

			base_unit = ""
			for counter, parm in enumerate(self.dimension_exponents.keys()):
				if not self.dimension_exponents[parm].is_integer:
					raise Exception("Can only handle integer exponents at the moment.")
				if counter > 0:
					base_unit += " + "

				if type(parm) is str:
					# this an extra param of type Constant
					base_unit += str(parm) + ".unit().symbols().get(0)"
				else:
					base_unit += "axis(" + str(parm) + ").unit().symbols().get(0)"

				separator = " " if counter < len(self.dimension_exponents) - 1 else ""
				if self.dimension_exponents[parm] != 1:
					base_unit += " + \"^" + str(int(self.dimension_exponents[parm])) + separator + "\""
				elif separator != "":
					base_unit += " + \" \""
		
			if len(self.dimension_exponents) == 1 and not "^" in base_unit:
				if type(parm) is str:
					self.code += "		map.put(" + str(dim) + ", " + str(parm) + ".unit());\n"
				else:
					self.code += "		map.put(" + str(dim) + ", axis(" + str(parm) + ").unit());\n"
			else:
				if base_unit == "": base_unit = "\"\""
				self.code += "		map.put(" + str(dim) + ", Unit.of(Units.simplify(" + base_unit + ")));\n"
			
		self.code += "		return Collections.unmodifiableNavigableMap(map);\n"
		self.extra_imports.append("import java.util.Collections;")
		#except Exception as e:
		#	print("error: " + str(e))
		#	self.code += "		// TODO: Implement\n"
		#	self.code += "		throw new UnsupportedOperationException(\"Failed to automatically derive the corresponding base units!\")"
		
		self.code += "	}\n"
		self.code += "\n"

	def __check_for_func_args(self, expr):
		"""
		Checks for specific functions that do not change the units of the expression, 
		which potentially allows to add a check in the constructor.

		Parameters
		----------
		expr : sympy expression
			The sympy expression to remove specific functions in
		"""
		trafo = str(expr)

		# try to recognize trigonometric functions
		# we know that the arguments of them need to be dimensionless, i.e. convertible to radians
		funcs = []
		
		# add math symbols
		funcs.append('asin')
		funcs.append('acos')
		funcs.append('atan')
		funcs.append('atan2')
		funcs.append('acot')
		funcs.append('asinh')
		funcs.append('acosh')
		funcs.append('sinh')
		funcs.append('cosh')
		funcs.append('tanh')
		funcs.append('sin')
		funcs.append('cos')
		funcs.append('tan')
		funcs.append('exp')
		funcs.append('log')
		funcs.append('factorial')

		symbols_str = [str(symbol) for symbol in self.symbols]

		for func in funcs:
			end_indices = [m.end() for m in re.finditer(func, trafo)]
			while len(end_indices) != 0:
				num_brackets = 0
				for counter, c in enumerate(trafo[end_indices[0]:]):
					if c == '(':
						num_brackets += 1
					elif c == ')':
						num_brackets -= 1
					if num_brackets == 0:
						bracket_expression = trafo[end_indices[0]+1:end_indices[0]+counter]
						if bracket_expression in symbols_str:
							index = symbols_str.index(bracket_expression)
							check = "		if (!Units.convertible(Axis.fromSet(axes, " + str(index) + ").unit(), Units.EMPTY_UNIT)) {\n"
							check += "			throw new IllegalArgumentException(\"The unit of dimension " + str(index) + " (%s) needs to be effectively dimensionless.\"\n"
							check += "					.formatted(Axis.fromSet(axes, " + str(index) + ").unit()));\n"
							check += "		}\n"
							self.__extra_symbol_check[index] = check
							self.extra_imports.append("import eu.hoefel.unit.Units;")
						trafo = trafo[:max(0, end_indices[0]-len(func)-1)] + trafo[end_indices[0]+counter+1:]  # -1 due to the preceding operator sign (will fail if the preceding operator is a '**'
						break
				end_indices = [m.end() for m in re.finditer(func, trafo)]

	def __remove_specific_funcs(self, expr):
		"""
		Removes specific functions (including their arguments) that do not change the units of the expression.

		Parameters
		----------
		expr : sympy expression
			The sympy expression to remove specific functions in

		Returns
		-------
		sympy expression
			the sympy expression with specific functions removed
		"""
		irrelevant_funcs_for_units = []
		
		irrelevant_funcs_for_units.append(sym.functions.elementary.trigonometric.asin)
		irrelevant_funcs_for_units.append(sym.functions.elementary.trigonometric.acos)
		irrelevant_funcs_for_units.append(sym.functions.elementary.trigonometric.atan)
		irrelevant_funcs_for_units.append(sym.functions.elementary.trigonometric.atan2)
		irrelevant_funcs_for_units.append(sym.functions.elementary.trigonometric.acot)
		irrelevant_funcs_for_units.append(sym.functions.elementary.hyperbolic.asinh)
		irrelevant_funcs_for_units.append(sym.functions.elementary.hyperbolic.acosh)
		irrelevant_funcs_for_units.append(sym.functions.elementary.hyperbolic.sinh)
		irrelevant_funcs_for_units.append(sym.functions.elementary.hyperbolic.cosh)
		irrelevant_funcs_for_units.append(sym.functions.elementary.hyperbolic.tanh)
		irrelevant_funcs_for_units.append(sym.functions.elementary.trigonometric.sin)
		irrelevant_funcs_for_units.append(sym.functions.elementary.trigonometric.cos)
		irrelevant_funcs_for_units.append(sym.functions.elementary.trigonometric.tan)
		irrelevant_funcs_for_units.append(sym.functions.elementary.exponential.exp)
		irrelevant_funcs_for_units.append(sym.functions.elementary.exponential.log)
		irrelevant_funcs_for_units.append(sym.functions.combinatorial.factorials.factorial)
		
		for irrelevant_func_for_units in irrelevant_funcs_for_units:
			# no idea why instanceof fails (╯°□°)╯︵ ┻━┻
			if str(expr.func) == str(irrelevant_func_for_units):
				return None
		args = []
		for arg in expr.args:
			if self.__remove_specific_funcs(arg) is not None:
				args.append(self.__remove_specific_funcs(arg))
		if args == []: return expr
		return expr.func(*args)

	def __remove_irrelevant_funcs(self, expr):
		"""
		Removes functions (not their arguments, though!) that do not change the units of the variables within.

		Parameters
		----------
		expr : sympy expression
			The sympy expression to remove irrelevant functions in

		Returns
		-------
		sympy expression
			the sympy expression with all irrelevant functions (with respect to units) removed
		"""
		irrelevant_funcs_for_units = []
		irrelevant_funcs_for_units.append(sym.functions.elementary.complexes.Abs)
		irrelevant_funcs_for_units.append(sym.functions.elementary.integers.floor)
		irrelevant_funcs_for_units.append(sym.functions.elementary.integers.ceiling)
		irrelevant_funcs_for_units.append(sym.functions.elementary.miscellaneous.Min)
		irrelevant_funcs_for_units.append(sym.functions.elementary.miscellaneous.Max)
		
		for irrelevant_func_for_units in irrelevant_funcs_for_units:
			# no idea why instanceof fails (╯°□°)╯︵ ┻━┻
			if str(expr.func) == str(irrelevant_func_for_units):
				return self.__remove_irrelevant_funcs(expr.args[0])
		
		args = []
		for arg in expr.args:
			args.append(self.__remove_irrelevant_funcs(arg))
		if args == []: return expr
		return expr.func(*args)
		
	def __collect_exponents(self, expr, symbols):
		"""
		Collects the exponents for the symbols in the given expression and adds them to the class-wide dict.

		Parameters
		----------
		expr : sympy expression
			The sympy expression to analyze.
		symbols : 
			The list of sympy symbols
		"""
		# most likely we have a mul object on the outside now
		# no idea why isinstance fails (╯°□°)╯︵ ┻━┻
		if str(expr.func) == str(sym.core.mul.Mul):
			for arg in expr.args:
				if str(arg.func) == str(sym.core.power.Pow):
					string_rep = sym.srepr(arg)
					# we should have smth akin to Pow(Symbol('x'),Integer(2))
					symbol, exponent = arg.args
					index = symbols.index(symbol)
					if '.' in str(exponent):
						if index in self.dimension_exponents:
							self.dimension_exponents[index] += float(str(exponent))
						else: 
							self.dimension_exponents[index] = float(str(exponent))
					else:
						if index in self.dimension_exponents:
							self.dimension_exponents[index] += float(str(exponent) + '.')
						else: 
							self.dimension_exponents[index] = float(str(exponent) + '.')
				elif str(arg.func) == str(sym.core.symbol.Symbol):
					symbols_str = [str(symbol) for symbol in symbols]
					if not str(arg) in symbols_str and str(arg) in self.extra_parameters and self.extra_parameters[str(arg)] == 'Constant':
						key = str(arg)
					else:
						key = symbols_str.index(str(arg))

					if key in self.dimension_exponents.keys():
						self.dimension_exponents[key] += 1.0
					else:
						self.dimension_exponents[key] = 1.0
		elif str(expr.func) == str(sym.core.power.Pow):
			string_rep = sym.srepr(expr)
			# we should have smth akin to Pow(Symbol('x'),Integer(2))
			symbol, exponent = expr.args
			if '.' in str(exponent):
				self.dimension_exponents[symbols.index(symbol)] = float(str(exponent))
			else:
				self.dimension_exponents[symbols.index(symbol)] = float(str(exponent) + '.')
		elif str(expr.func) == str(sym.core.symbol.Symbol):
			symbols_str = [str(symbol) for symbol in symbols]
			if not str(expr) in symbols_str and str(expr) in self.extra_parameters and self.extra_parameters[str(expr)] == 'Constant':
				self.dimension_exponents[str(expr)] = 1.0
			else:
				self.dimension_exponents[symbols_str.index(str(expr))] = 1.0
		else:
			print("---------")
			print("expr is " + str(expr))
			print("expr func is " + str(expr.func))
			print("Found type: " + str(type(expr.func)) + " for " + str(expr))
			print("---------")
			raise Exception("Cannot handle unit transformation (maybe too complex equation?)")

	def __keep_one_summand(self, expr):
		"""
		Keeps one summand only.

		Parameters
		----------
		expr : sympy expression
			The sympy expression to keep one summand only.

		Returns
		-------
		sympy expression
			the sympy expression with all except one summand removed
		"""
		args = []
		if isinstance(expr, sym.core.add.Add):
			return expr.args[0]
		else:
			for arg in expr.args:
				if isinstance(arg, sym.core.add.Add):
					args.append(arg.args[0])
				else:
					args.append(self.__keep_one_summand(arg))
		if args == []: return expr
		return expr.func(*args)
	
	def __build_to_base_transformation_point(self):
		"""
		Attaches the Java code for a method allowing to transform a point in the current coordinates to the corresponding base coordinate system.
		"""
		self.code += "	@Override\n"
		self.code += "	public double[] toBasePoint(double[] position) {\n"
		
		if self.__use_validate_position:
			self.code += "		validatePosition(position);\n\n"

		if self.to_base_point is None:
			self.code += "		throw new UnsupportedOperationException(\"Transformation not implemented!\")"
		else:
			self.code += "		double[] pointInBase = new double[" + str(self.g.dim()) + "];\n"
			for i in range(self.g.dim()):
				code, extra_imports = self.symbolic_to_java(self.to_base_point[i],self.symbols, extra_symbols=self.extra_parameters)
				self.code += "		pointInBase[" + str(i) + "] = " + code + ";\n"
				self.extra_imports.extend(extra_imports)
			self.code += "		return pointInBase;\n"
		self.code += "	}\n"
		self.code += "\n"

	def __build_from_base_transformation_point(self):
		"""
		Attaches the Java code for a method allowing to transform a point in corresponding base coordinates to the current coordinate system.
		"""
		self.code += "	@Override\n"
		self.code += "	public double[] fromBasePoint(double[] position) {\n"

		if not self.unable_to_invert and self.from_base_point is None:
			self.unable_to_invert, self.from_base_point, self.base_symbols = self.__invert_forward_transformation()

		if self.to_base_point is None or self.unable_to_invert:
			self.code += "		throw new UnsupportedOperationException(\"Transformation not implemented!\")\n"
		else:
			self.code += "		double[] pointInCurrentSystem = new double[" + str(self.dim) + "];\n"
			for i in range(self.dim):
				code, extra_imports = self.symbolic_to_java(self.from_base_point[i],self.base_symbols, extra_symbols=self.extra_parameters)
				self.code += "		pointInCurrentSystem[" + str(i) + "] = " + code + ";\n"
				self.extra_imports.extend(extra_imports)
			self.code += "		return pointInCurrentSystem;\n"
		self.code += "	}\n"
		self.code += "\n"

	def __invert_forward_transformation(self):
		"""
		Inverts the given to_base_point transformation, i.e. it finds the transformation for a point given 
		in base coordinates to the current coordinate system.

		Returns
		-------
		unable_to_invert
			True if there was any failure in inverting
		self.from_base_point
			The resulting inverse transformation or None if not successful
		self.base_symbols
			The symbols corresponding to the inverse transformation
		"""
		if self.from_base_point is None:
			unable_to_invert = False
			try:
				if self.base_symbols is None:
					print("changing base symbols...")
					self.base_symbols = [i.as_dummy() for i in self.g.symbols]
				print("Inverting forward transformation...")
				inv_results = sym.solve([t[0] - t[1] for t in zip(self.base_symbols, self.to_base_point)], list(self.symbols), dict=True)[0]
				self.from_base_point = [sym.trigsimp(inv_results[s].simplify()) for s in self.symbols]
			except Exception as e:
				unable_to_invert = True
				self.from_base_point = None
				self.base_symbols = None
				print("inversion failed: " + str(e))
		else:
			unable_to_invert = False  # effectively false as given by user or previously calculated

		return unable_to_invert, self.from_base_point, self.base_symbols

	def __build_to_base_transformation_vector(self):
		"""
		Attaches the Java code to transform a vector from the current coordinate system to the corresponding base coordinate system.
		"""
		self.code += "	@Override\n"
		self.code += "	public double[] toBaseVector(double[] position, double[] vector) {\n"
		
		if self.__use_validate_position:
			self.code += "		validatePosition(position);\n\n"
		
		if self.to_base_point is None:
			self.code += "		throw new UnsupportedOperationException(\"Transformation not implemented!\")"
		else:
			self.code += "		double[] vectorInBaseSys = new double[vector.length];\n"
			
			vec = []
			for i in range(self.dim):			
				for j in range(self.dim):
					expr = self.to_base_point[j].diff(self.symbols[i]).simplify()
					if expr == 0: continue
					
					if expr == 1:
						expr_as_string = "vector[" + str(j) + "]"
					else:
						code, extra_imports = self.symbolic_to_java(expr, self.symbols, extra_symbols=self.extra_parameters)
						expr_as_string = code + "*vector[" + str(j) + "]"
						self.extra_imports.extend(extra_imports)

					if len(vec) <= i:
						vec.append(expr_as_string)
					else:
						vec[i] += "+" + expr_as_string
				if vec[i] != 0: 
					self.code += "		vectorInBaseSys[" + str(i) + "] = " + vec[i] + ";\n"
			self.code += "		return vectorInBaseSys;\n"
		self.code += "	}\n"
		self.code += "\n"

	def __build_from_base_transformation_vector(self):
		"""
		Attaches the Java code to transform a vector from the corresponding base coordinate system to the current coordinate system.
		"""
		self.code += "	@Override\n"
		self.code += "	public double[] fromBaseVector(double[] position, double[] vector) {\n"
		
		if not self.unable_to_invert and self.from_base_point is None:
			self.unable_to_invert, self.from_base_point, self.base_symbols = self.__invert_forward_transformation()
		
		if self.to_base_point is None or self.unable_to_invert:
			self.code += "		throw new UnsupportedOperationException(\"Transformation not implemented!\")"
		else:
			self.code += "		double[] vectorInCurrentSys = new double[vector.length];\n"
			
			vec = []
			for i in range(self.g.dim()):			
				for j in range(self.g.dim()):
					expr = self.from_base_point[j].diff(self.base_symbols[i])
					if len(str(expr)) < 75:
						expr = expr.simplify()  # sympy might get stuck (or just take forever) for too long expressions
					if expr == 0: continue
					
					if expr == 1:
						expr_as_string = "vector[" + str(j) + "]"
					else:
						code, extra_imports = self.symbolic_to_java(expr,self.base_symbols, extra_symbols=self.extra_parameters)
						expr_as_string = code + "*vector[" + str(j) + "]"
						self.extra_imports.extend(extra_imports)

					if len(vec) <= i:
						vec.append(expr_as_string)
					elif not expr_as_string.startswith('-'):
						vec[i] += "+" + expr_as_string
					else:
						vec[i] += expr_as_string
				if vec[i] != 0: 
					self.code += "		vectorInCurrentSys[" + str(i) + "] = " + vec[i] + ";\n"
			self.code += "		return vectorInCurrentSys;\n"
		self.code += "	}\n"
		self.code += "\n"

	def __build_metric_tensor(self):
		"""
		Attaches the Java code for the method to calculate the metric tensor to the code.
		"""
		self.code += "	@Override\n"
		self.code += "	public double[][] metricTensor(double[] position, TensorTransformation behavior) {\n"
		if self.__use_validate_position:
			self.code += "		validatePosition(position);\n\n"

		self.code += "		int dim = " + str(self.g.dim()) + ";\n"
		self.code += "		double[][] g = new double[dim][dim];\n"
		self.code += "\n"
		self.code += "		// Note that we skip all elements that are zero anyways\n"
		
		# check whether the same elements in covariant and covariant entries are nonnull
		sameElementsAreNonNull = True
		for i in range(self.g.dim()):
				for j in range(self.g.dim()):
					if (self.g.gdd[i,j] != 0 and self.g.guu[i,j] == 0) or (self.g.gdd[i,j] == 0 and self.g.guu[i,j] != 0):
						sameElementsAreNonNull = False
				else:
					# Continue if the inner loop wasn't broken.
					continue
				# Inner loop was broken, break the outer.
				break
		
		if self.g.gdd != self.g.guu and not sameElementsAreNonNull:
			self.code += "		return switch (behavior) {\n"
			self.code += "			case COVARIANT -> {\n"
			for i in range(self.g.dim()):
				for j in range(self.g.dim()):
					if self.g.gdd[i,j] == 0: continue
					self.code += "				g[" + str(i) + "][" + str(j) + "] = metricCoefficient(position, behavior, " + str(i) + ", " + str(j) +");\n"
			self.code += "				return g;\n"
			self.code += "			}\n"
			
			self.code += "			case CONTRAVARIANT -> {\n"
			for i in range(self.g.dim()):
				for j in range(self.g.dim()):
					if self.g.guu[i,j] == 0: continue
					self.code += "				g[" + str(i) + "][" + str(j) + "] = metricCoefficient(position, behavior, " + str(i) + ", " + str(j) +");\n"
			self.code += "				return g;\n"
			self.code += "			}\n"
			self.code += "		};\n"
			self.code += "	}\n"
		else:
			for i in range(self.g.dim()):
				for j in range(self.g.dim()):
					if self.g.gdd[i,j] == 0: continue
					self.code += "		g[" + str(i) + "][" + str(j) + "] = metricCoefficient(position, behavior, " + str(i) + ", " + str(j) +");\n"
			self.code += "		return g;\n"
			self.code += "	}\n"
		self.code += "\n"

	def __build_metric_coefficient(self):
		"""
		Attaches the Java code for a method to get the metric coefficients to the code.
		"""
		self.extra_imports.append("import eu.hoefel.tensors.TensorTransformation;")
		self.extra_imports.append("import eu.hoefel.tensors.TensorIndexType;")
		self.extra_imports.append("import java.util.function.Function;")
	
		self.code += "	@Override\n"
		self.code += "	public double metricCoefficient(double[] position, TensorTransformation behavior, int i, int j) {\n"
		
		if self.__use_validate_position:
			self.code += "		validatePosition(position);\n\n"
		
		self.code += "		if (i < 0 || j < 0) {\n"
		self.code += "			throw new IllegalArgumentException((\"Metric coefficient not available for i=%d, j=%d \"\n"
		self.code += "					+ \"(too low dimension, only " + str(self.g.dim()) + " dimensions " + str([*range(self.g.dim())]).replace(" ","") + " are supported for " + self.java_name.lower() + " coordinates)\")\n"
		self.code += "					.formatted(i, j));\n"
		self.code += "		} else if (i > " + str(self.g.dim()-1) + " || j > " + str(self.g.dim()-1) + ") {\n"
		self.code += "			throw new IllegalArgumentException((\"Metric coefficient not available for i=%d, j=%d \"\n"
		self.code += "					+ \"(too high dimension, only " + str(self.g.dim()) + " dimensions " + str([*range(self.g.dim())]).replace(" ","") + " are supported for " + self.java_name.lower() + " coordinates)\")\n"
		self.code += "					.formatted(i, j));\n"
		self.code += "		}\n\n"
		
		self.code += "		if (behavior instanceof TensorIndexType tit) {\n"
		if self.g.gdd != self.g.guu:
			self.code += "			return switch (tit) {\n"
			self.code += "				case COVARIANT -> {\n"
			for i in range(self.g.dim()):
				for j in range(self.g.dim()):
					if self.g.gdd[i,j] == 0: continue
					code, extra_imports = self.symbolic_to_java(self.g.gdd[i,j], self.symbols, extra_symbols=self.extra_parameters)
					self.code += "					if (i == " + str(i) + " && j == " + str(j) + ") yield " + code + ";\n"
					self.extra_imports.extend(extra_imports)
			self.code += "					yield 0;\n"
			self.code += "				}\n"
			
			self.code += "				case CONTRAVARIANT -> {\n"
			for i in range(self.g.dim()):
				for j in range(self.g.dim()):
					if self.g.guu[i,j] == 0: continue
					code, extra_imports = self.symbolic_to_java(self.g.guu[i,j], self.symbols, extra_symbols=self.extra_parameters)
					self.code += "					if (i == " + str(i) + " && j == " + str(j) + ") yield " + code + ";\n"
					self.extra_imports.extend(extra_imports)
			self.code += "					yield 0;\n"
			self.code += "				}\n"
			self.code += "			};\n"
		else:
			for i in range(self.g.dim()):
				for j in range(self.g.dim()):
					if self.g.gdd[i,j] == 0: continue
					code, extra_imports = self.symbolic_to_java(self.g.gdd[i,j], self.symbols, extra_symbols=self.extra_parameters)
					self.code += "			if (i == " + str(i) + " && j == " + str(j) + ") return " + code + ";\n"
					self.extra_imports.extend(extra_imports)
			self.code += "			return 0;\n"
		self.code += "		}\n\n"
		
		self.code += "		// Fall back to metric tensor (which might fall back to this method) for mixed behavior\n"
		self.code += "		Function<double[], double[][]> metricTensor = pos -> metricTensor(pos, TensorIndexType.COVARIANT);\n"
		self.code += "		return TensorIndexType.COVARIANT.transform(this, metricTensor, behavior).apply(position)[i][j];\n"
		self.code += "	}\n\n"

	def __build_jacobian_determinant(self):
		"""
		Attaches the Java code for a method determining the Jacobi determinant to the code.
		"""
		self.code += "	@Override\n"
		self.code += "	public double jacobianDeterminant(double[] position) {\n"
		
		if self.__use_validate_position:
			self.code += "		validatePosition(position);\n"
		
		det = sym.sqrt(self.g.gdd.det()).simplify()

		code, extra_imports = self.symbolic_to_java(det, self.symbols, extra_symbols=self.extra_parameters)
		self.extra_imports.extend(extra_imports)

		self.code += "		return " + code + ";\n"
		self.code += "	}\n"
		self.code += "\n"

	def __build_dimension(self):
		"""
		Attaches the Java code for the method returning the maximum allowed dimensionality to the code.
		"""
		self.code += "	@Override\n"
		self.code += "	public int dimension() {\n"
		self.code += "		return " + str(self.g.dim()) + ";\n"
		self.code += "	}\n"
		self.code += "\n"

	def __build_isOrthogonal(self):
		"""
		Attaches the Java code for the method to determine whether the CoordinateSystem is orthogonal to the code.
		"""
		self.code += "	@Override\n"
		self.code += "	public boolean isOrthogonal() {\n"
		
		isOrthogonal = True
		for i in range(self.g.dim()):
			for j in range(self.g.dim()):
				if i == j: continue
				if self.g.gdd[i,j] != 0:
					isOrthogonal = False
					break
			else:
				# Continue if the inner loop wasn't broken.
				continue
			# Inner loop was broken, break the outer.
			break
		
		self.code += "		return " + str(isOrthogonal).lower() + ";\n"
		self.code += "	}\n"
		self.code += "\n"

	def __build_christoffel_symbols_1st_kind(self):
		"""
		Attaches the Java code of the method calculating the Christoffel symbols of the first kind to the code.
		"""
		self.code += "	@Override\n"
		self.code += "	public double christoffelSymbol1stKind(double[] position, int i, int j, int k) {\n"
		
		if self.__use_validate_position:
			self.code += "		validatePosition(position);\n\n"
			
		self.code += "		int dim = position.length;\n"
		self.code += "		if (i < 0 || j < 0 || k < 0 || i >= dim || j >= dim || k >= dim) {\n"
		self.code += "			throw new IllegalArgumentException(\n"
		self.code += "					\"i, j and k may not be <0 or exceed %d, but they were i=%d, j=%d and k=%d\"\n"
		self.code += "							.formatted(dim, i, j, k));\n"
		self.code += "		}\n\n"


		for i in range(self.g.dim()):
			for j in range(self.g.dim()):
				for k in range(self.g.dim()):
					christoffelSymbol1stKind = 0
					for u in range(self.g.dim()):
						christoffelSymbol1stKind += self.g.gdd[u,k] * self.gamma.udd(u,i,j)
					christoffelSymbol1stKind = christoffelSymbol1stKind.simplify()
				
					if christoffelSymbol1stKind == 0: continue
					code, extra_imports = self.symbolic_to_java(christoffelSymbol1stKind, self.symbols, extra_symbols=self.extra_parameters)
					self.code += "		if (i == " + str(i) + " && j == " + str(j) + " && k == " + str(k) + ") return " + code + ";\n"
					self.extra_imports.extend(extra_imports)
		self.code += "		return 0;\n"
		self.code += "	}\n"
		self.code += "\n"

	def __build_christoffel_symbols_2nd_kind(self):
		"""
		Attaches the Java code of the method calculating the Christoffel symbols of the second kind to the code.
		"""
		self.code += "	@Override\n"
		self.code += "	public double christoffelSymbol2ndKind(double[] position, int m, int i, int j) {\n"
		
		if self.__use_validate_position:
			self.code += "		validatePosition(position);\n\n"
			
		self.code += "		int dim = position.length;\n"
		self.code += "		if (m < 0 || i < 0 || j < 0 || m >= dim || i >= dim || j >= dim) {\n"
		self.code += "			throw new IllegalArgumentException(\n"
		self.code += "					\"m, i and j may not be <0 or exceed %d, but they were m=%d, i=%d and j=%d\"\n"
		self.code += "							.formatted(dim, m, i, j));\n"
		self.code += "		}\n\n"

		for m in range(self.g.dim()):
			for i in range(self.g.dim()):
				for j in range(self.g.dim()):
					if self.gamma.udd(m,i,j) == 0: continue
					code, extra_imports = self.symbolic_to_java(self.gamma.udd(m,i,j), self.symbols, extra_symbols=self.extra_parameters)
					self.code += "		if (m == " + str(m) + " && i == " + str(i) + " && j == " + str(j) + ") return " + code + ";\n"
					self.extra_imports.extend(extra_imports)
		self.code += "		return 0;\n"
		self.code += "	}\n"
		self.code += "\n"

	def __build_riemann_tensor(self, mute=True):
		"""
		Attaches the Java code that calculates a component of the Riemann tensor R^mu_nu,rho,sigma to the code.
		
		Parameters
		----------
		mute : bool, optional
			Whether to print to sys out
		"""
		self.code += "	@Override\n"
		self.code += "	public double riemannTensor(double[] position, int mu, int nu, int rho, int sigma) {\n"
		
		if self.__use_validate_position:
			self.code += "		validatePosition(position);\n\n"

		self.code += "		int dim = position.length;\n"
		self.code += "		if (mu < 0 || nu < 0 || rho < 0 || sigma < 0 || mu >= dim || nu >= dim || rho >= dim || sigma >= dim) {\n"
		self.code += "			throw new IllegalArgumentException(\n"
		self.code += "					\"mu, nu, rho and sigma may not be <0 or exceed %d, but they were mu=%d, nu=%d, rho=%d and sigma=%d\"\n"
		self.code += "							.formatted(dim, mu, nu, rho, sigma));\n"
		self.code += "		}\n\n"

		self.__Riemann_is_zero = True
		for i in range(self.dim):
			for j in range(self.dim):
				for k in range(self.dim):
					for l in range(self.dim):
						expr = self.__simplify(self.riem.uddd(i,j,k,l))
						if expr == 0:
							continue
						else:
							# So it seems we have a nonzero expression
							# However, in most instances sympy probably just fails to realize that the expression is zero
							# To try to realize this, we will try a 1000 random points within the specified limits
							# If they are all zero, most likely the expression is just 0
							# In this case we keep the original expressions added as a comment and add a TODO to allow
							# manual checking by the user (e.g. via wolfram alpha or by hand), but by default return just 0.
							if not mute: 
								print("Detected a nonzero component of the Riemann tensor (for R^" + str(i) + "_" + str(j) + "," + str(k) + "," + str(l) + ")")
								print("Often, this is due to the automatic routines failing to recognize that the expression is in fact zero.")
								print("For this reason, a number of points are checked numerically, if all of them are zero we assume that it is indeed a simplification problem.")
								print("In this case, a TODO is added to the generated class, though.")
								print("Note that if a term of the Riemann tensor turns out to be nonzero the isFlat method should return false.")
							
							conc_symbols = self.__concatenate_symbols()
							lambda_expr = sym.lambdify(conc_symbols, expr, "numpy")

							points = self.__generate_points(conc_symbols)
							deviates_notable_from_zero = False
							for point in self.progressbar(points, prefix="Checking nullness of R^" + str(i) + "_" + str(j) + "," + str(k) + "," + str(l) + ": ", mute=mute):
								if lambda_expr(*point) > 1e-6:
									deviates_notable_from_zero = True
									break
							
							code, extra_imports = self.symbolic_to_java(expr, self.symbols, extra_symbols=self.extra_parameters)
							
							if deviates_notable_from_zero:
								if not mute: print("Found a clear deviation from zero, keeping expression.")
								self.extra_imports.extend(extra_imports)
								self.code += "		if (mu == " + str(i) + " && nu == " + str(j) + " && rho == " + str(k) + " && sigma == " + str(l) + ") return " + code + ";\n"
								self.__Riemann_is_zero = False
							else:
								if not mute: print("Expression appears to be zero. Adding TODO and continuing.\n")
								self.code += "		// TODO: The expression listed below appears to be zero.\n"
								self.code += "		// However, to be sure please check manually.\n"
								self.code += "		// Expression in sympy: " + str(expr) + "\n"
								self.code += "		// The Java code is listed below (incl. indices): \n"
								self.code += "		// if (mu == " + str(i) + " && nu == " + str(j) + " && rho == " + str(k) + " && sigma == " + str(l) + ") return " + code + ";\n\n"
						
		self.code += "		return 0;\n"
		self.code += "	}\n"
		self.code += "\n"

	def __concatenate_symbols(self):
		"""
		Concatenates the symbols with potentially available extra symbols.

		Returns
		-------
		list
			the concatenated sympy symbols
		"""
		conc_symbols = self.symbols.copy()
		for extra_param in self.extra_parameters:
			conc_symbols.append(sym.sympify(extra_param))
			
		return conc_symbols

	def __generate_points(self, symbols, n=1000):
		"""
		Generates points, respecting the specified limits.
		By default, if no limits a re given for a symbol, a value will 
		be taken from the range -100 to 100.

		Parameters
		----------
		n : int, optional
			The number of points to create

		Returns
		-------
		list
			the created points, i.e. a list of tuples
		"""
		points = []
		for i in range(n):
			point = []
			for symbol in symbols:
				symbol_str = str(symbol)
				low = 100
				high = 100
				if symbol_str in self.variable_lower_limit:
					low = self.__pythonize(self.variable_lower_limit[symbol_str])
				if symbol_str in self.variable_upper_limit:
					high = self.__pythonize(self.variable_upper_limit[symbol_str])
				point.append(float(random.uniform(low, high)))
			points.append(tuple(point))
		return points

	def __pythonize(self, value, default=100):
		"""
		Small method to get strings that potentially contain pi or e 
		expressed in java code to python understandable values.

		Parameters
		----------
		value : str or float
			The value to sanitize for python
		default : float
			The default to use for infinity (and -infinity with a 
			prepended minus sign)

		Returns
		-------
		float
			the value as a python float
		"""
		try:
			val = float(value)
			return val
		except Exception:
			# not a float, obviously
			value_str = str(value)
			value_str = value_str.replace("Math.PI", "pi")
			value_str = value_str.replace("Math.E", "e")
			value_str = value_str.replace("Double.POSITIVE_INFINITY", str(default))
			value_str = value_str.replace("Double.NEGATIVE_INFINITY", str(-default))
			return float(sym.sympify(value_str).evalf())

	@staticmethod
	def progressbar(it, prefix="", size=50, mute=False):
		"""
		Small method to provide a progressbar. 
		See https://stackoverflow.com/a/34482761

		Parameters
		----------
		it : iterable
			The iterable to iterate over
		prefix : str, optional
			The string to write before the progressbar.
		size : int, optional
			The width of the progressbar
		mute : bool, optional
			Whether to print anything to sys out

		Returns
		-------
		item
			the current item of the iterable
		"""
		count = len(it)
		def show(j):
			x = int(size*j/count)
			if not mute:
				sys.stdout.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
				sys.stdout.flush()		
		show(0)
		for i, item in enumerate(it):
			yield item
			show(i+1)
		if not mute:
			sys.stdout.write("\n")
			sys.stdout.flush()

	def __build_isFlat(self):
		"""
		Attaches the Java code that checks whether the current coordinate system is flat to the code.
		"""
		self.code += "	@Override\n"
		self.code += "	public boolean isFlat() {\n"
		self.code += "		return " + str(self.__Riemann_is_zero).lower() + ";\n"
		self.code += "	}\n"
		self.code += "\n"

	def __build_ricci_tensor(self, mute=True):
		"""
		Attaches the Java code that calculates a component of the Ricci tensor R_mu,nu to the code.

		Parameters
		----------
		mute : bool, optional
			Whether to print anything to sys out
		"""
		self.code += "	@Override\n"
		self.code += "	public double ricciTensor(double[] position, int mu, int nu) {\n"
		
		if self.__use_validate_position:
			self.code += "		validatePosition(position);\n\n"
		
		self.code += "		int dim = position.length;\n"
		self.code += "		if (mu < 0 || nu < 0 || mu >= dim || nu >= dim) {\n"
		self.code += "			throw new IllegalArgumentException(\n"
		self.code += "					\"mu and nu may not be <0 or exceed %d, but they were mu=%d and nu=%d\".formatted(dim, mu, nu));\n"
		self.code += "		}\n\n"
		  
		for i in range(self.dim):
			for j in range(self.dim):
				expr = self.ric.dd(i,j)
				
				if expr == 0:
					continue
				else:
					# So it seems we have a nonzero expression
					# However, in most instances sympy probably just fails to realize that the expression is zero
					# To try to realize this, we will try a 1000 random points within the specified limits
					# If they are all zero, most likely the expression is just 0
					# In this case we keep the original expressions added as a comment and add a TODO to allow
					# manual checking by the user (e.g. via wolfram alpha or by hand), but by default return just 0.
					if not mute: 
						print("Detected a nonzero component of the Ricci tensor (for R_" + str(i) + "," + str(j) + ")")
						print("Often, this is due to the automatic routines failing to recognize that the expression is in fact zero.")
						print("For this reason, a number of points are checked numerically, if all of them are zero we assume that it is indeed a simplification problem.")
						print("In this case, a TODO is added to the generated class, though.")
					
					conc_symbols = self.__concatenate_symbols()
					lambda_expr = sym.lambdify(conc_symbols, expr, "numpy")

					points = self.__generate_points(conc_symbols)
					deviates_notable_from_zero = False
					for point in self.progressbar(points, prefix="Checking nullness of R_" + str(i) + "," + str(j) + ": ", mute=mute):
						if lambda_expr(*point) > 1e-6:
							deviates_notable_from_zero = True
							break
					
					code, extra_imports = self.symbolic_to_java(self.ric.dd(i,j), self.symbols, extra_symbols=self.extra_parameters)
					
					if deviates_notable_from_zero:
						if not mute: print("Found a clear deviation from zero, keeping expression.")
						self.extra_imports.extend(extra_imports)
						self.code += "		if (mu == " + str(i) + " && nu == " + str(j) + ") return " + code + ";\n"
					else:
						if not mute: print("Expression appears to be zero. Adding TODO and continuing.\n")
						self.code += "		// TODO: The expression listed below appears to be zero.\n"
						self.code += "		// However, to be sure please check manually.\n"
						self.code += "		// Expression in sympy: " + str(expr) + "\n"
						self.code += "		// The Java code is listed below (incl. indices): \n"
						self.code += "		// if (mu == " + str(i) + " && nu == " + str(j) + ") return " + code + ";\n\n"

		self.code += "		return 0;\n"
		self.code += "	}\n"
		self.code += "\n"

	def __build_ricci_scalar(self, mute=True):
		"""
		Attachs the Java code to calculate the Ricci scalar.
		
		Parameters
		----------
		mute : bool, optional
			Whether to print anything to sys out
		"""
		self.code += "	@Override\n"
		self.code += "	public double ricciScalar(double[] position) {\n"
		
		if self.__use_validate_position:
			self.code += "		validatePosition(position);\n"

		expr = self.ric_scalar.value()
		
		if expr == 0:
			self.code += "		return 0;\n"
		else:
			# So it seems we have a nonzero expression
			# However, in most instances sympy probably just fails to realize that the expression is zero
			# To try to realize this, we will try a 1000 random points within the specified limits
			# If they are all zero, most likely the expression is just 0
			# In this case we keep the original expressions added as a comment and add a TODO to allow
			# manual checking by the user (e.g. via wolfram alpha or by hand), but by default return just 0.
			if not mute: 
				print("Detected a nonzero Ricci scalar.")
				print("Often, this is due to the automatic routines failing to recognize that the expression is in fact zero.")
				print("For this reason, a number of points are checked numerically, if all of them are zero we assume that it is indeed a simplification problem.")
				print("In this case, a TODO is added to the generated class, though.")
			
			conc_symbols = self.__concatenate_symbols()
			lambda_expr = sym.lambdify(conc_symbols, expr, "numpy")

			points = self.__generate_points(conc_symbols)
			deviates_notable_from_zero = False
			for point in self.progressbar(points, prefix="Checking nullness of R: ", mute=mute):
				if lambda_expr(*point) > 1e-6:
					deviates_notable_from_zero = True
					break
			
			code, extra_imports = self.symbolic_to_java(expr, self.symbols, extra_symbols=self.extra_parameters)
			
			if deviates_notable_from_zero:
				if not mute: print("Found a clear deviation from zero, keeping expression.")
				self.extra_imports.extend(extra_imports)
				self.code += "		return " + code + ";\n"
			else:
				if not mute: print("Expression appears to be zero. Adding TODO and continuing.\n")
				self.code += "		// TODO: The expression listed below appears to be zero.\n"
				self.code += "		// However, to be sure please check manually.\n"
				self.code += "		// Expression in sympy: " + str(expr) + "\n"
				self.code += "		// The Java code is listed below (incl. indices): \n"
				self.code += "		// return " + code + ";\n\n"
				self.code += "		return 0;\n"

		self.code += "	}\n"

	@staticmethod
	def symbolic_to_java(expr, symbols, java_name='position', outermost=True, extra_symbols=None):
		"""
		Converts the given mathematical expression to the corresponding Java expression.

		Parameters
		----------
		expr : symypy expression
			The expression to be converted, which may contain symbols 
			and, e.g., trigonometric functions.
		symbols : numpy array or array_like
			The symbols representing the different dimensions, e.g. `sympy.symbols('r theta z')`.
		java_name : str
			The name of the array that will be used to replace the known symbols, 
			with the index of the array corresponding to the index of the symbol in 
			the list of symbols.
		outermost : bool, optional
			The boolean deciding if brackets have to be added in some circumstances. If the expression 
			is the outermost one, brackets can be avoided.
		extra_symbols : numpy array or array_like
			Further symbols. Note that any occurences of these symbols will use the exact string 
			representation of the symbol.

		Returns
		-------
		str, list
			the expression using standard Java where possible and eu.hoefel.utils.Maths where necessary. 
			Furthermore, the list containing additional imports is provided. 
		"""
		stringified = ""
		all_imports = []

		# no idea why isinstance fails (╯°□°)╯︵ ┻━┻
		if str(expr.func) == str(sym.core.mul.Mul):
			factors = []
			for arg in expr.args:
				factor, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(arg, symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
				all_imports.extend(extra_imports)
				factors.append(factor)
			stringified = "*".join(factors)
		elif str(expr.func) == str(sym.core.add.Add):
			summands = []
			for counter, arg in enumerate(expr.args):
				summand, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(arg, symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
				all_imports.extend(extra_imports)
				if counter > 0 and not summand.startswith("-"): stringified += "+"
				stringified += summand
				summands.append(summand)
			if not outermost: stringified = "(" + stringified + ")"
		elif str(expr.func) == str(sym.core.power.Pow):
			base, exponent = expr.args
			base_symbol, base_imports = JavaCoordinateSystemCreator.symbolic_to_java(base, symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			exponent_symbol, exponent_imports = JavaCoordinateSystemCreator.symbolic_to_java(exponent, symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(base_imports)
			all_imports.extend(exponent_imports)
			stringified = "Math.pow(" +  base_symbol + "," + exponent_symbol + ")"
		elif str(expr.func) == str(sym.core.symbol.Symbol):
			symbols_str = [str(symbol) for symbol in symbols]
			if str(expr) in symbols_str:
				index = symbols_str.index(str(expr))
				stringified = java_name + "[" + str(index) + "]"
			elif extra_symbols is not None and str(expr) in extra_symbols:
				stringified = str(expr)
				if extra_symbols[str(expr)] == 'Constant': stringified += ".value()"
			else:
				raise Exception("Unknown symbol: " + str(expr))
		elif str(expr.func) == str(sym.core.numbers.Integer):
			stringified = str(expr)
		elif str(expr.func) == str(sym.functions.elementary.complexes.Abs):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.abs(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.complexes.sign):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.signum(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.sin):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.sin(" + x + ")"	
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.cos):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.cos(" + x + ")"		
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.tan):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.tan(" + x + ")"	  
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.cot):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "1/Math.tan(" + x + ")"			
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.sec):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "1/Math.cos(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.csc):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "1/Math.sin(" + x + ")"   
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.sinc):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.sinc(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.asin):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.asin(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.acos):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.acos(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.atan):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.atan(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.acot):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.acot(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.asec):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.acos(1.0/(" + x + "))"
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.acsc):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.asin(1.0/(" + x + "))"
		elif str(expr.func) == str(sym.functions.elementary.trigonometric.atan2):
			y, y_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			x, x_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(y_extra_imports)
			all_imports.extend(x_extra_imports)
			stringified = "Math.atan2(" + y + "," + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.sinh):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.sinh(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.cosh):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.cosh(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.tanh):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.tanh(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.coth):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.coth(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.sech):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "1/Math.cosh(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.csch):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "1/Math.sinh(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.asinh):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.arsinh(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.acosh):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.arcosh(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.atanh):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.artanh(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.acoth):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.arcoth(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.asech):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.arsech(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.hyperbolic.acsch):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Maths.arcsch(" + x + ")"
		elif str(expr.func) == str(sym.functions.elementary.integers.ceiling):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.ceil(" + x + ")"	
		elif str(expr.func) == str(sym.functions.elementary.integers.floor):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.floor(" + x + ")"  
		elif str(expr.func) == str(sym.functions.elementary.integers.RoundFunction):
			raise Exception("Roundfunction is currently not supported")
		elif str(expr.func) == str(sym.functions.elementary.integers.floor):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = x + "-Math.floor(" + x + ")" 
			if not outermost: stringified = "(" + stringified + ")" 
		elif str(expr.func) == str(sym.functions.elementary.exponential.exp):
			x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.exp(" + x + ")" 
		elif str(expr.func) == str(sym.functions.elementary.exponential.LambertW):
			raise Exception("LambertW is currently not supported")
		elif str(expr.func) == str(sym.functions.elementary.exponential.log):
			if len(expr.args) == 2:
				x, x_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
				b, b_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
				all_imports.extend(x_extra_imports)
				all_imports.extend(b_extra_imports)
				stringified = "Math.log(" + x + ")/Math.log(" + b + ")"
			else:
				x, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
				all_imports.extend(extra_imports)
				stringified = "Math.log(" + x + ")" 
		elif str(expr.func) == str(sym.functions.elementary.exponential.exp_polar):
			raise Exception("exp_polar is currently not supported")
		elif str(expr.func) == str(sym.functions.elementary.miscellaneous.IdentityFunction):
			stringified, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
		elif str(expr.func) == str(sym.functions.elementary.miscellaneous.Min):
			a, a_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			b, b_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(a_extra_imports)
			all_imports.extend(b_extra_imports)
			stringified = "Math.min(" + a + "," + b + ")"
		elif str(expr.func) == str(sym.functions.elementary.miscellaneous.Max):
			a, a_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			b, b_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(a_extra_imports)
			all_imports.extend(b_extra_imports)
			stringified = "Math.max(" + a + "," + b + ")"
		elif str(expr.func) == str(sym.functions.elementary.miscellaneous.root):
			if len(expr.args) != 2: raise Exception("root is currently only supported with 2 arguments")
			a, a_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			b, b_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(a_extra_imports)
			all_imports.extend(b_extra_imports)
			stringified = "Math.pow(" + a + ",1.0/" + b + ")"
		elif str(expr.func) == str(sym.functions.elementary.miscellaneous.sqrt):
			sqrt, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.sqrt(" + sqrt + ")"
		elif str(expr.func) == str(sym.functions.elementary.miscellaneous.cbrt):
			cbrt, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = "Math.pow(" + cbrt + ",1.0/3.0)"
		elif str(expr.func) == str(sym.functions.elementary.miscellaneous.real_root):
			raise Exception("real_root is currently not supported")
		elif str(expr.func) == str(sym.core.numbers.Float) or str(expr.func) == str(sym.core.numbers.RealNumber):
			stringified = str(expr)
		elif str(expr.func) == str(sym.core.numbers.Rational):
			stringified = str(expr.p) + ".0/" + str(expr.q)
		elif str(expr.func) == str(sym.core.numbers.Zero):
			stringified = "0"
		elif str(expr.func) == str(sym.core.numbers.One):
			stringified = "1"
		elif str(expr.func) == str(sym.core.numbers.NegativeOne):
			stringified = "-1"
		elif str(expr.func) == str(sym.core.numbers.Half):
			stringified = "0.5"
		elif str(expr.func) == str(sym.core.numbers.NaN):
			stringified = "Double.NaN"
		elif str(expr.func) == str(sym.core.numbers.Infinity):
			stringified = "Double.POSITIVE_INFINITY"
		elif str(expr.func) == str(sym.core.numbers.NegativeInfinity):
			stringified = "Double.NEGATIVE_INFINITY"
		elif str(expr.func) == str(sym.core.numbers.Exp1):
			stringified = "Math.E"
		elif str(expr.func) == str(sym.core.numbers.Pi):
			stringified = "Math.PI"
		elif str(expr.func) == str(sym.functions.special.delta_functions.DiracDelta):
			a, extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(extra_imports)
			stringified = a + " == 0 ? 1 : 0"
			if not outermost: stringified = "(" + stringified + ")"
		elif str(expr.func) == str(sym.functions.elementary.piecewise.Piecewise):
			pieces = expr.args
			for counter, piece in enumerate(pieces):
				expression, expression_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(piece.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
				condition, condition_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(piece.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
				all_imports.extend(expression_extra_imports)
				all_imports.extend(condition_extra_imports)
				
				if condition == 'true':
					# skip unnecessary stuff
					stringified += expression
					break
				elif condition != 'false':
					try:
						float(expression)
						stringified += condition + " ? " + expression + " : "
					except:
						stringified += condition + " ? (" + expression + ") : "
				if counter < len(pieces) - 1:
					stringified += "("
				else:
					stringified += "Double.NaN"
			stringified += (len(pieces)-1) * ")"
			stringified = "(" + stringified + ")"
		elif str(expr.func) == str(sym.core.relational.Eq) or str(expr.func) == str(sym.core.relational.Equality):
			lhs, lhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			rhs, rhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(lhs_extra_imports)
			all_imports.extend(rhs_extra_imports)
			stringified = lhs + "==" + rhs
			#if not outermost: stringified = "(" + stringified + ")"
		elif str(expr.func) == str(sym.core.relational.Ne) or str(expr.func) == str(sym.core.relational.Unequality):
			lhs, lhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			rhs, rhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(lhs_extra_imports)
			all_imports.extend(rhs_extra_imports)
			stringified = lhs + "!=" + rhs
			if not outermost: stringified = "(" + stringified + ")"
		elif str(expr.func) == str(sym.core.relational.Lt) or str(expr.func) == str(sym.core.relational.StrictLessThan):
			lhs, lhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			rhs, rhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(lhs_extra_imports)
			all_imports.extend(rhs_extra_imports)
			stringified = lhs + "<" + rhs
			if not outermost: stringified = "(" + stringified + ")"
		elif str(expr.func) == str(sym.core.relational.Le) or str(expr.func) == str(sym.core.relational.LessThan):
			lhs, lhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			rhs, rhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(lhs_extra_imports)
			all_imports.extend(rhs_extra_imports)
			stringified = lhs + "<=" + rhs
			if not outermost: stringified = "(" + stringified + ")"
		elif str(expr.func) == str(sym.core.relational.Gt) or str(expr.func) == str(sym.core.relational.StrictGreaterThan):
			lhs, lhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			rhs, rhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(lhs_extra_imports)
			all_imports.extend(rhs_extra_imports)
			stringified = lhs + ">" + rhs
			if not outermost: stringified = "(" + stringified + ")"
		elif str(expr.func) == str(sym.core.relational.Ge) or str(expr.func) == str(sym.core.relational.GreaterThan):
			lhs, lhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[0], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			rhs, rhs_extra_imports = JavaCoordinateSystemCreator.symbolic_to_java(expr.args[1], symbols, java_name=java_name, outermost=False, extra_symbols=extra_symbols)
			all_imports.extend(lhs_extra_imports)
			all_imports.extend(rhs_extra_imports)
			stringified = lhs + ">=" + rhs
			if not outermost: stringified = "(" + stringified + ")"
		elif str(expr.func) == str(sym.logic.boolalg.BooleanFalse):
			stringified = "false"
		elif str(expr.func) == str(sym.logic.boolalg.BooleanTrue):
			stringified = "true"
		else:
			print("Found no implementation for " + str(expr))

		if 'Maths' in stringified:
			all_imports.append("import eu.hoefel.utils.Maths;")

		return stringified, all_imports

if __name__ == "__main__":
	# TODO vector trafo

	r = sym.symbols('r', real=True, positive=True)
	theta = sym.symbols('theta', real=True)
	syms = sym.symbols('r theta z', real=True)
	
	axes = {0: "SiBaseUnit.METER", 1: "rad", 2: "rad"}
	extra_parms = {"a": "double"}

	m = Metric.fromTransformation(syms, to_base_point=[syms[0]*sym.cos(syms[1]),syms[0]*sym.sin(syms[1]),syms[2]])
	
	code = JavaCoordinateSystemCreator("CylindricalCoordinates", m, axes, extra_parameters=extra_parms)
	code.variable_lower_limit['a'] = 2.0
	code.variable_upper_limit['a'] = 3.0
	

	code.variable_lower_limit['r'] = 0
	code.variable_upper_limit['theta'] = 'Math.PI'
	
	print(code)
	print("----")
	
	#m = Metric.fromTransformation(syms, to_base_point=[syms[0],syms[1],syms[2]])
	#print(JavaCoordinateSystemCreator("CartesianCoordinates", m, axes))
	#print("----")
	#
	#denominator = (syms[0]**2 + syms[1]**2 + syms[2]**2)
	#m = Metric.fromTransformation(syms, to_base_point=[syms[0] / denominator, syms[1] / denominator, syms[2] / denominator])
	#print(JavaCoordinateSystemCreator("SixSphereCoordinates", m, axes))
	#print("----")
	