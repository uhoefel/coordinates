package eu.hoefel.coordinates.axes;

import java.util.Collections;
import java.util.Comparator;
import java.util.NavigableSet;
import java.util.Objects;
import java.util.TreeSet;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.logging.Logger;
import java.util.stream.Stream;

import eu.hoefel.coordinates.CoordinateSystem;
import eu.hoefel.unit.Unit;

/**
 * Helper class for operations that include multiple {@link Axis axes}. This is
 * mainly intended for use in implementations of {@link CoordinateSystem}.
 * 
 * @author Udo Hoefel
 */
public final class Axes {

	/** The units of the dimensions. */
	private final Unit[] units;

	/** The indicator whether the default dimensions should also be overridden. */
	private final boolean overrideDefault;

	/**
	 * The default target dimension for an {@link Axis}. Usually, this should put
	 * the axis at the first entry in a {@link NavigableSet} using the
	 * {@link #AXIS_DIM_COMPARATOR}.
	 */
	public static final int DEFAULT_DIMENSION = -1;

	/**
	 * A default argument checker. Mainly intended for use together with
	 * {@link #fromArgs(NavigableSet, Consumer, Object...)}, such that the user only
	 * needs to provide a name (e.g. "Cartesian") to obtain the necessary consumer.
	 * Checks if all given arguments are either {@link Axis} oder {@link Axes}.
	 */
	public static final BiConsumer<String, Object[]> DEFAULT_ARG_CHECK = (name, args) -> {
		for (Object arg : args) {
			if (!(arg instanceof Axis) && !(arg instanceof Axes)) {
				throw new IllegalArgumentException(("%s coordinates can handle i) no arguments (default) "
						+ "and ii) objects that are either Axes or Axis objects. However, you passed on a %s object.")
								.formatted(name, arg.getClass().getSimpleName()));
			}
		}
	};

	/** Compares two axis by comparing their target dimension. */
	public static final Comparator<Axis> AXIS_DIM_COMPARATOR = (axis1, axis2) -> {
		if (axis1.dimension() > axis2.dimension()) return 1;
		if (axis1.dimension() < axis2.dimension()) return -1;
		return 0;
	};

	private static final Logger logger = Logger.getLogger(Axes.class.getName());

	/**
	 * Constructor.
	 * 
	 * @param units           the units of the axes, e.g. "m"
	 * @param extraUnits      additional units required to parse the given units
	 * @param overrideDefault if true, indicates that the default axes should also
	 *                        be changed
	 */
	private Axes(String[] units, Unit[][] extraUnits, boolean overrideDefault) {
		this(Stream.of(units).map(unit -> Unit.of(unit, extraUnits)).toArray(Unit[]::new), overrideDefault);
	}

	/**
	 * Constructor.
	 * 
	 * @param units           the units of the axes, e.g.
	 *                        {@link eu.hoefel.unit.si.SiBaseUnit#METER}
	 * @param overrideDefault if true, indicates that the default axes should also
	 *                        be changed
	 */
	private Axes(Unit[] units, boolean overrideDefault) {
		this.units = units;
		this.overrideDefault = overrideDefault;
	}

	/**
	 * Constructs an axes object that contains the units for all axes, except for
	 * the default axes.
	 * 
	 * @param units the units for the axes
	 * @return the axes with the specified units
	 */
	public static final Axes withUnits(Unit... units) {
		return new Axes(units, false);
	}

	/**
	 * Constructs an axes object that contains the units for all axes, except for
	 * the default axes.
	 * 
	 * @param units      the units for the axes
	 * @param extraUnits additional units required for the parsing of {@code units}
	 * @return the axes with the specified units
	 */
	public static final Axes withUnits(String[] units, Unit[]... extraUnits) {
		return new Axes(units, extraUnits, false);
	}

	/**
	 * Constructs an axes object that contains the units for all axes, including the
	 * default axes.
	 * 
	 * @param units the units for the axes
	 * @return the axes with the specified units
	 */
	public static final Axes allWithUnits(Unit... units) {
		return new Axes(units, true);
	}

	/**
	 * Constructs an axes object that contains the units for all axes, including the
	 * default axes.
	 * 
	 * @param units      the units for the axes
	 * @param extraUnits additional units required for the parsing of {@code units}
	 * @return the axes with the specified units
	 */
	public static final Axes allWithUnits(String[] units, Unit[]... extraUnits) {
		return new Axes(units, extraUnits, true);
	}

	/**
	 * Parses the given arguments and constructs a set of axes from them.
	 * 
	 * @param defaultAxes the default axes to use
	 * @param args        the arguments from which to construct the axes
	 * @return the sorted axes
	 */
	private static final NavigableSet<Axis> parseArgs(NavigableSet<Axis> defaultAxes, Object... args) {
		Objects.requireNonNull(defaultAxes);
		Objects.requireNonNull(args);

		var axes = new TreeSet<>(AXIS_DIM_COMPARATOR);
		axes.addAll(defaultAxes);

		for (Object arg : args) {
			if (arg instanceof Axis axis) {
				axes.removeIf(val -> val.dimension() == axis.dimension());
				axes.add(axis);
			}
		}

		for (Object arg : args) {
			if (arg instanceof Axes a) {
				if (a.overrideDefault) {
					// we override everything in order
					int counter = 0;
					for (Axis axis : axes) {
						axes.removeIf(val -> val.dimension() == axis.dimension());
						axes.add(new Axis(axis.dimension(), a.units[counter++]));
					}
				} else {
					Axis startAxis = axes.higher(Axis.forDim(Axes.DEFAULT_DIMENSION));
					if (startAxis != null) {
						var axesToChange = new TreeSet<>(axes.tailSet(startAxis));
						int counter = 0;
						for (Axis axis : axesToChange) {
							axes.removeIf(val -> val.dimension() == axis.dimension());
							axes.add(new Axis(axis.dimension(), a.units[counter++]));
						}
					}
				}
			}
		}
		
		return axes;
	}

	/**
	 * Creates a set of axes from the given arguments. Intended as a helper method
	 * in the constructor when implementing new coordinate systems. Does not check
	 * the inputs, in contrast to {@link #fromArgs(NavigableSet, Consumer,
	 * Object...)}.
	 * 
	 * @param defaultAxes the default axes to use
	 * @param args        the (input) arguments from which to construct the axes
	 * @return the sorted axes (usually of a coordinate system)
	 */
	public static final NavigableSet<Axis> fromArgs(NavigableSet<Axis> defaultAxes, Object... args) {
		return parseArgs(defaultAxes, args);
	}

	/**
	 * Creates a set of axes from the given arguments. Intended as a helper method
	 * in the constructor when implementing new coordinate systems.
	 * 
	 * @param defaultAxes the default axes to use
	 * @param argCheck    performs a check if the arguments match the expectations
	 * @param args        the (input) arguments from which to construct the axes
	 * @return the sorted axes (usually of a coordinate system)
	 */
	public static final NavigableSet<Axis> fromArgs(NavigableSet<Axis> defaultAxes, Consumer<Object[]> argCheck, Object... args) {
		Objects.requireNonNull(defaultAxes);
		Objects.requireNonNull(argCheck);
		Objects.requireNonNull(args);

		argCheck.accept(args);
		return parseArgs(defaultAxes, args);
	}

	/**
	 * Creates an unmodifiable set of the given axes. If an axis is added with a
	 * target dimension identical to an already added axis a warning is issued. This
	 * is mainly intended to be used within constructors of {@link CoordinateSystem}
	 * implementations.
	 * 
	 * @param axes the axes to put together
	 * @return the set of axes
	 */
	public static final NavigableSet<Axis> of(Axis... axes) {
		Objects.requireNonNull(axes);

		NavigableSet<Axis> ret = new TreeSet<>(AXIS_DIM_COMPARATOR);
		for (Axis axis : axes) {
			if (!ret.add(axis)) {
				logger.warning(() -> "Overriding already existing axis for dimension %d!".formatted(axis.dimension()));
			}
		}
		return Collections.unmodifiableNavigableSet(ret);
	}
}
