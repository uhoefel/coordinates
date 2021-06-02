package eu.hoefel.coordinates;

import java.lang.reflect.AnnotatedType;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Stream;

import eu.hoefel.coordinates.axes.Axes;
import eu.hoefel.coordinates.axes.Axis;
import eu.hoefel.coordinates.tensors.TensorIndexType;
import eu.hoefel.coordinates.tensors.TensorTransformation;
import eu.hoefel.unit.Unit;
import eu.hoefel.unit.Units;
import eu.hoefel.unit.constant.Constant;
import eu.hoefel.utils.Maths;
import eu.hoefel.utils.Types;

/**
 * Helper class for coordinate systems, i.e. it contains method to deal with
 * coordinate transformations, and also provides convenience methods to extract
 * e.g. doubles from objects given to a constructor.
 * 
 * @author Udo Hoefel
 */
public final class CoordinateSystems {

	/**
	 * The "identity" coordinate system, i.e. a coordinate system that keeps the
	 * coordinates as given.
	 */
	public static final CoordinateSystem IDENTITY_COORDINATE_SYSTEM;

	/**
	 * The default coordinates used within {@link CoordinateSystem}. These do not
	 * include the {@link #IDENTITY_COORDINATE_SYSTEM}. Might also be useful within
	 * (some) implementations of {@link CoordinateSystem}.
	 */
	public static final Set<Class<? extends CoordinateSystem>> DEFAULT_COORDINATE_SYSTEMS;

	/** A set of axes containing only one unitless default axis. */
	private static final NavigableSet<Axis> IDENTITY_AXES = Axes.of(new Axis(Axes.DEFAULT_DIMENSION, Units.EMPTY_UNIT, ""));
	
	static {
		IDENTITY_COORDINATE_SYSTEM = new @CoordinateSystemSymbols("") CoordinateSystem() {
			@Override public int dimension() { return Maths.MAX_ARRAY_SIZE; }
			@Override public NavigableSet<Axis> axes() { return IDENTITY_AXES; }
			@Override public boolean isBasic() { return true; }
			@Override public Class<? extends CoordinateSystem> baseCoordinates() { return this.getClass(); }
			@Override public double[] toBasePoint(double[] coords) { return coords.clone(); }
			@Override public double[] fromBasePoint(double[] coords) { return coords.clone(); }
			@Override public boolean isOrthogonal() { return true; }
			@Override public double[][] metricTensor(double[] coords, TensorTransformation behavior) { return Maths.eye(coords.length); }
			@Override public double metricCoefficient(double[] coords, TensorTransformation behavior, int i, int j) { return i == j ? 1 : 0; }
			@Override public boolean isFlat() { return true; }
			@Override
			public NavigableMap<Integer, Unit> toBaseUnits() {
				var map = Map.of(Axes.DEFAULT_DIMENSION, Units.EMPTY_UNIT);
				return Collections.unmodifiableNavigableMap(new TreeMap<>(map));
			}
			
		};

		Set<Class<? extends CoordinateSystem>> defCoordinates = new LinkedHashSet<>();
		defCoordinates.add(CartesianCoordinates.class);
		defCoordinates.add(PolarCoordinates.class);
		defCoordinates.add(CylindricalCoordinates.class);
		defCoordinates.add(SphericalCoordinates.class);
		defCoordinates.add(ToroidalCoordinates.class);
		defCoordinates.add(SixSphereCoordinates.class);
		defCoordinates.add(BipolarCoordinates.class);
		defCoordinates.add(BisphericalCoordinates.class);
		defCoordinates.add(OblateSpheroidalCoordinates.class);
		defCoordinates.add(ProlateSpheroidalCoordinates.class);
		DEFAULT_COORDINATE_SYSTEMS = Collections.unmodifiableSet(defCoordinates);
	}

	/** Class is only a utility class! */
	private CoordinateSystems() {
	    throw new IllegalStateException("This is a pure utility class!");
	}

	/**
	 * Transforms the given coordinates from the original coordinate system to the
	 * target coordinate system. Note that the base of the original coordinate
	 * system must be identical to the base of the target coordinate system. If one
	 * of the coordinate systems requires state and has no default constructor you
	 * have to wrap it in {@link CoordinateSystem#from(String, Object...)} or
	 * {@link CoordinateSystem#from(String, Set, Object...)} and supply it to
	 * {@link #transform(double[], CoordinateSystem, CoordinateSystem)}.
	 * 
	 * @param coords the coordinates in the original coordinate system
	 * @param origin the original coordinate system
	 * @param target the target coordinate system
	 * @return the coordinates in the target coordinate system
	 */
	public static final double[] transform(double[] coords, String origin, String target) {
		return transform(coords, origin, target, DEFAULT_COORDINATE_SYSTEMS);
	}

	/**
	 * Transforms the given coordinates from the original coordinate system to the
	 * target coordinate system. Note that the base of the original coordinate
	 * system must be identical to the base of the target coordinate system.
	 * 
	 * @param coords      the coordinates in the original coordinate system
	 * @param origin      the original coordinate system
	 * @param target      the target coordinate system
	 * @param extraCoords additional coordinates needed to parse the original and/or
	 *                    target coordinates
	 * @return the coordinates in the target coordinate system
	 */
	public static final double[] transform(double[] coords, String origin, String target,
			Set<Class<? extends CoordinateSystem>> extraCoords) {
		// I would prefer a vararg method here, but that would lead to potential heap
		// pollution
		CoordinateSystem originalCoords = CoordinateSystem.from(origin, extraCoords);
		CoordinateSystem targetCoords = CoordinateSystem.from(target, extraCoords);
		return transform(coords, originalCoords, targetCoords);
	}

	/**
	 * Transforms the given coordinates from the original coordinate system to the
	 * target coordinate system. Note that the base of the original coordinate
	 * system must be identical to the base of the target coordinate system.
	 * 
	 * @param position the coordinates in the original coordinate system, not null
	 * @param origin   the original coordinate system, not null
	 * @param target   the target coordinate system, not null
	 * @return the coordinates in the target coordinate system
	 */
	public static final double[] transform(double[] position, CoordinateSystem origin, CoordinateSystem target) {
		Objects.requireNonNull(position);
		Objects.requireNonNull(origin);
		Objects.requireNonNull(target);

		if (origin.equals(target)) {
			// shortcut for performance
			return position;
		}

		checkDimensionality(position.length, origin, target);
		checkBaseCoordinates(origin, target);
		checkCoordinateSystemUnits(origin, target);

		double[] positionInOriginBase = origin.toBasePoint(position);
		
		// convert to base units of target system
		double[] positionInTargetBase = new double[positionInOriginBase.length];
		for (int i = 0; i < positionInTargetBase.length; i++) {
			// due to the check made above we know that we can convert and either a specific
			// or the default unit is given
			var originBaseUnits = origin.toBaseUnits();
			var targetBaseUnits = target.toBaseUnits();
			Unit originUnit = originBaseUnits.getOrDefault(i, originBaseUnits.get(Axes.DEFAULT_DIMENSION));
			Unit targetUnit = targetBaseUnits.getOrDefault(i, targetBaseUnits.get(Axes.DEFAULT_DIMENSION));
			positionInTargetBase[i] = Units.convert(positionInOriginBase[i], originUnit, targetUnit);
		}
		
		return target.fromBasePoint(positionInTargetBase);
	}

	/**
	 * Checks whether the given dimensionality is compatible with the specified
	 * coordinate systems.
	 * 
	 * @param dimension the dimensionality of e.g. a point
	 * @param origin    the coordinate system of origin
	 * @param target    the target coordinate system
	 */
	private static void checkDimensionality(int dimension, CoordinateSystem origin, CoordinateSystem target) {
		if (origin.dimension() < dimension) {
			throw new IllegalArgumentException("Dimension of given position exceeds maximum dimension "
					+ "handleable by the coordinate system of origin (%d vs %d)!".formatted(dimension,
							origin.dimension()));
		} else if (target.dimension() < dimension) {
			throw new IllegalArgumentException("Dimension of given position exceeds maximum dimension "
					+ "handleable by the target coordinate system (%d vs %d)!".formatted(dimension,
							target.dimension()));
		}
	}

	/**
	 * Checks whether the base coordinate systems of the two given coordinate
	 * systems are the same.
	 * 
	 * @param origin the coordinate system of origin
	 * @param target the target coordinate system
	 * @throws IllegalArgumentException if the base coordinate are not the same
	 */
	private static final void checkBaseCoordinates(CoordinateSystem origin, CoordinateSystem target) {
		if (origin.baseCoordinates() != target.baseCoordinates()) {
			throw new IllegalArgumentException(
					("Incompatible base coordinate systems: coordinate system of origin is %s (base: %s) "
							+ "vs. target coordinate system %s (base: %s).")
					.formatted(origin.getClass().getSimpleName(), origin.baseCoordinates().getSimpleName(), 
							   target.getClass().getSimpleName(), target.baseCoordinates().getSimpleName()));
		}
	}

	/**
	 * Checks whether the coordinate system units are convertible at each step.
	 * 
	 * @param origin the coordinate system of origin
	 * @param target the target coordinate system
	 * @throws IllegalArgumentException if the dimensions of the coordinate systems
	 *                                  are not the same (and cannot be recovered by
	 *                                  the default axes), or the units of a
	 *                                  specific dimension are not convertible, or
	 *                                  the units of the default dimension are not
	 *                                  convertible
	 */
	private static final void checkCoordinateSystemUnits(CoordinateSystem origin, CoordinateSystem target) {
		var originBaseUnits = origin.toBaseUnits();
		var targetBaseUnits = target.toBaseUnits();

		// the default axes units must match if default axes are present in both systems
		boolean originHasDefault = originBaseUnits.containsKey(Axes.DEFAULT_DIMENSION);
		boolean targetHasDefault = targetBaseUnits.containsKey(Axes.DEFAULT_DIMENSION);

		if (originHasDefault && targetHasDefault && !Units.convertible(originBaseUnits.get(Axes.DEFAULT_DIMENSION),
						targetBaseUnits.get(Axes.DEFAULT_DIMENSION))) {
			throw new IllegalArgumentException(("Units of default axes do not match (%s vs %s). "
					+ "A match is required to ensure that arbitrarily-dimensional coordinate systems work.").formatted(
							originBaseUnits.get(Axes.DEFAULT_DIMENSION), targetBaseUnits.get(Axes.DEFAULT_DIMENSION)));
		}

		// all the non-default axes units must match
		int minNondefaultDimension = Axes.DEFAULT_DIMENSION + 1;
		var originNondefaultBaseUnits = originBaseUnits.tailMap(minNondefaultDimension);
		var targetNondefaultBaseUnits = targetBaseUnits.tailMap(minNondefaultDimension);

		int maxNondefaultDimension = Math.max(originNondefaultBaseUnits.size(), targetNondefaultBaseUnits.size());

		for (int i = 0; i < maxNondefaultDimension; i++) {
			Unit originUnit;
			if (i < originNondefaultBaseUnits.size() || originHasDefault) {
				originUnit = originNondefaultBaseUnits.getOrDefault(i, originBaseUnits.get(Axes.DEFAULT_DIMENSION));
			} else {
				throw new IllegalArgumentException(
						"%dth dimension is not defined in the original coordinate system, but is required to convert to the target coordinate system."
								.formatted(i));
			}

			Unit targetUnit;
			if (i < targetNondefaultBaseUnits.size() || targetHasDefault) {
				targetUnit = targetNondefaultBaseUnits.getOrDefault(i, targetBaseUnits.get(Axes.DEFAULT_DIMENSION));
			} else {
				throw new IllegalArgumentException(
						"%dth dimension is not defined in the target coordinate system, but is required to convert to from the original coordinate system."
								.formatted(i));
			}

			if (!Units.convertible(originUnit, targetUnit)) {
				throw new IllegalArgumentException("Units of axes at %dth dimension do not match (%s vs %s)."
						.formatted(i, originUnit.symbols().get(0), targetUnit.symbols().get(0)));
			}
		}
	}

	/**
	 * Gets the symbols representing the coordinate system from the class.
	 * 
	 * @param clazz the class of the coordinate system
	 * @return the symbols
	 */
	static final List<String> symbolsFromClass(Class<? extends CoordinateSystem> clazz) {
		CoordinateSystemSymbols cs = clazz.getAnnotation(CoordinateSystemSymbols.class);
		if (cs == null) {
			// maybe an anonymous inner class? we check parent superclass and interfaces
			AnnotatedType as = clazz.getAnnotatedSuperclass();
			cs = as.getAnnotation(CoordinateSystemSymbols.class);
			if (cs == null) {
				AnnotatedType[] ais = clazz.getAnnotatedInterfaces();
				for (AnnotatedType ai : ais) {
					cs = ai.getAnnotation(CoordinateSystemSymbols.class);
					if (cs != null) break; 
				}
			}
		}
		if (cs == null) return List.of();
		return List.of(cs.value());
	}
	
	/**
	 * Gets the <code>n</code>th double from the given <code>args</code>.
	 * <p>
	 * This is intended to be used to extract doubles from arguments passed on to an
	 * implementation of a {@link CoordinateSystem coordinate system}. If the
	 * implementation requires multiple doubles, the <code>n</code> parameter allows
	 * to extract the <code>n</code>th double from the given <code>args</code>, even
	 * if, e.g., an {@link Axis} object is given between the doubles.
	 * 
	 * @param n    the index of the double requested, starting with 0
	 * @param args the arguments to search for the nth double
	 * @return the double
	 */
	public static final Optional<Double> doubleFromArgs(int n, Object... args) {
	    return Stream.of(args)
	                 .filter(CoordinateSystems::isConvertibleToDouble)
	                 .skip(n)
	                 .map(CoordinateSystems::toDouble)
	                 .findFirst();
	}

	/**
	 * Checks whether the given object is a double, or can be widened to a double,
	 * or is a String that can be parsed to a double.
	 * 
	 * @param o the object
	 * @return true if the object can be sensibly converted into a double
	 */
	public static final boolean isConvertibleToDouble(Object o) {
		return Types.isCompatible(double.class, o) || (o instanceof String s && Maths.isDouble(s));
	}

	/**
	 * Gets the double form the given object if the given object is a double, or can
	 * be widened to a double, or is a String that can be parsed to a double.
	 * 
	 * @param o the object
	 * @return the double
	 */
	public static final double toDouble(Object o) {
	    Objects.requireNonNull(o);
        
        if (Types.isCompatible(double.class, o)) {
            if (Types.boxedClass(o.getClass()) == Character.class) {
                return (double) (char) o;
            }

            Object boxed = Types.box(o);
            if (boxed instanceof Number number) {
                return number.doubleValue();
            }
        } else if (o instanceof String s && Maths.isDouble(s)) {
            return Double.parseDouble(s);
        }

        throw new IllegalArgumentException(o + " is not an double and cannot be converted to one!");
	}

	/**
	 * Gets the <code>n</code>th int from the given <code>args</code>.
	 * <p>
	 * This is intended to be used to extract ints from arguments passed on to an
	 * implementation of a {@link CoordinateSystem coordinate system}. If the
	 * implementation requires multiple ints, the <code>n</code> parameter allows to
	 * extract the <code>n</code>th int from the given <code>args</code>, even if,
	 * e.g., an {@link Axis} object is given between the ints. Cannot be combined
	 * with {@link #doubleFromArgs(int, Object...)}, as it recognizes ints as
	 * doubles as well.
	 * 
	 * @param n    the index of the int requested, starting with 0
	 * @param args the arguments to search for the nth int
	 * @return the int
	 */
	public static final Optional<Integer> intFromArgs(int n, Object... args) {
		int counter = 0;
		for (Object o : args) {
			if (isConvertibleToInt(o) && counter++ == n) return Optional.of(toInt(o));
		}
		return Optional.empty();
	}

	/**
	 * Checks whether the given object is an int, or can be widened to an int,
	 * or is a String that can be parsed to an int.
	 * 
	 * @param o the object
	 * @return true if the object can be sensibly converted into an int
	 */
	public static final boolean isConvertibleToInt(Object o) {
		return Types.isCompatible(int.class, o) || (o instanceof String s && Maths.isInteger(s));
	}

	/**
	 * Gets the int form the given object if the given object is an int, or can
	 * be widened to an int, or is a String that can be parsed to an int.
	 * 
	 * @param o the object
	 * @return the int
	 */
	public static final int toInt(Object o) {
	    Objects.requireNonNull(o);
	    
	    if (Types.isCompatible(int.class, o)) {
	        if (Types.boxedClass(o.getClass()) == Character.class) {
                return (int) (char) o;
            }

	        Object boxed = Types.box(o);
	        if (boxed instanceof Number number) {
	            return number.intValue();
	        }
        } else if (o instanceof String s && Maths.isInteger(s)) {
            return Integer.parseInt(s);
        }

		throw new IllegalArgumentException(o + " is not an int and cannot be converted to one!");
	}

	/**
	 * Gets the <code>n</code>th Constant from the given <code>args</code>.
	 * <p>
	 * This is intended to be used to extract a {@link Constant} from arguments
	 * passed on to an implementation of a {@link CoordinateSystem coordinate
	 * system}. If the implementation requires multiple constants, the
	 * <code>n</code> parameter allows to extract the <code>n</code>th constant from
	 * the given <code>args</code>, even if, e.g., an {@link Axis} object is given
	 * between the constants.
	 * <p>
	 * <em>Important note</em>:<br>
	 * If a coordinate system requires both extra parameters that are of type double
	 * and Constant, using {@link #doubleFromArgs(int, Object...)} in conjunction
	 * with {@link #constantFromArgs(int, Object...)} will potentially yield
	 * incorrect results if the doubles and constants are given as strings, each
	 * containing only a double. If a unit is specified in the string, they are safe
	 * to use.
	 * 
	 * @param n    the index of the Constant requested, starting with 0
	 * @param args the arguments to search for the nth constant
	 * @return the Constant
	 */
	public static final Optional<Constant> constantFromArgs(int n, Object... args) {
		int counter = 0;
		for (Object o : args) {
			if (o instanceof Constant c && counter++ == n) {
				return Optional.of(c);
			} else if (o instanceof String s && Constant.isConstant(s) && counter++ == n) {
				return Optional.of(Constant.of(s));
			} else if (isConvertibleToDouble(o) && counter++ == n) {
				return Optional.of(Constant.of(toDouble(o)));
			}
		}
		return Optional.empty();
	}

	/**
	 * Prints the metric tensor of the given coordinate system at the specified
	 * position.
	 * 
	 * @param sys      the coordinate system
	 * @param position the position
	 */
	public static final void printMetricTensor(CoordinateSystem sys, double[] position) {
		double[][] metricTensor = sys.metricTensor(position, TensorIndexType.COVARIANT);
		printMatrix("g", metricTensor);
	}

	/**
	 * Prints the Christoffel symbols of the second kind for the given index.
	 * 
	 * @param sys      the coordinate system
	 * @param position the position
	 * @param m        the index of the first tangent direction (that is not
	 *                 derived)
	 */
	public static final void printChristoffelSymbols2ndKind(CoordinateSystem sys, double[] position, int m) {
		double[][] christoffels = new double[position.length][position.length];
		for (int i = 0; i < position.length; i++) {
			for (int j = 0; j < position.length; j++) {
				christoffels[i][j] = sys.christoffelSymbol2ndKind(position, m, i, j);
			}
		}
		printMatrix("Γ^%d_ij".formatted(m), christoffels);
	}

	/**
	 * Prints the Christoffel symbols of the first kind for the given index.
	 * 
	 * @param sys      the coordinate system
	 * @param position the position
	 * @param k        the <i>k</i> index in [<i>ij</i>,<i>k</i>]
	 */
	public static final void printChristoffelSymbols1stKind(CoordinateSystem sys, double[] position, int k) {
		double[][] christoffels = new double[position.length][position.length];
		for (int i = 0; i < position.length; i++) {
			for (int j = 0; j < position.length; j++) {
				christoffels[i][j] = sys.christoffelSymbol1stKind(position, k, i, j);
			}
		}
		printMatrix("[ij,%d]".formatted(k), christoffels);
	}

	/**
	 * Prints the given matrix.
	 * 
	 * @param leftHandSide the left hand side, i.e. the symbol that represents the
	 *                     given values
	 * @param values       the matrix
	 */
	private static final void printMatrix(String leftHandSide, double[][] values) {
		boolean isOdd = Maths.isOdd(values.length);
		for (int i = 0; i < values.length; i++) {
			StringBuilder sb = new StringBuilder();
			
			if ((isOdd && i == (values.length - 1.0) / 2) || (!isOdd && i == values.length / 2)) {
				sb.append(leftHandSide);
				sb.append(" = ");
			} else {
				sb.append(" ".repeat(leftHandSide.length()));
				sb.append("   ");
			}
			
			if (i == 0) {
				sb.append("⎛");
			} else if (i == values.length - 1) {
				sb.append("⎝");
			} else {
				sb.append("⎜");
			}

			for (int j = 0; j < values[i].length; j++) {
				sb.append(" %10.5f ");
			}

			if (i == 0) {
				sb.append("⎞");
			} else if (i == values.length - 1) {
				sb.append("⎠");
			} else {
				sb.append("⎟");
			}

			System.out.println(String.format(Locale.ENGLISH, sb.toString(), (Object[]) Types.box(values[i])));
		}
	}
}
