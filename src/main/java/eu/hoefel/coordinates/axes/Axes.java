package eu.hoefel.coordinates.axes;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.NavigableSet;
import java.util.Objects;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.logging.Logger;
import java.util.stream.Stream;

import eu.hoefel.coordinates.CoordinateSystem;
import eu.hoefel.unit.Unit;
import eu.hoefel.unit.UnitContext;
import eu.hoefel.unit.UnitContextMatch;
import eu.hoefel.unit.Units;

/**
 * Helper class for operations that include multiple {@link Axis axes}. This is
 * mainly intended for use in implementations of {@link CoordinateSystem}.
 * 
 * @author Udo Hoefel
 */
public final class Axes {

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

    /**
     * Compares two axis by comparing first their target dimension, then, if
     * necessary, their first symbol and finally, if necessary, their axis name.
     */
    public static final Comparator<Axis> AXIS_DIM_COMPARATOR = Comparator.comparing(Axis::dimension)
            .thenComparing(Axis::unit, (u1, u2) -> Units.COMPARATOR_FOR_UNIT_ORDERING.compare(u1.symbols().get(0), u2.symbols().get(0)))
            .thenComparing(Axis::name);

    /** The units of the dimensions. */
    private final Unit[] units;

    /** The context in which to use the units to try to infer the axis name. */
    private final UnitContext[] contexts;

    /** The axis names. */
    private final String[] names;

    /** The indicator whether the default dimensions should also be overridden. */
    private final boolean overrideDefault;

    /** An empty set representing no axis names. */
    private static final NavigableSet<String> EMPTY_AXIS_NAMES = Collections.unmodifiableNavigableSet(new TreeSet<>());

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
        this(Stream.of(units).map(unit -> Unit.of(unit, extraUnits)).toArray(Unit[]::new), null, null, overrideDefault);
    }

    /**
     * Constructor.
     * 
     * @param units           the units of the axes, e.g.
     *                        {@link eu.hoefel.unit.si.SiBaseUnit#METER}
     * @param names           the names of the axes (e.g. "length")
     * @param contexts        the contexts that may be used to allo axis name
     *                        inference in combination with {@code units}
     * @param overrideDefault if true, indicates that the default axes should also
     *                        be changed
     */
    private Axes(Unit[] units, String[] names, UnitContext[] contexts, boolean overrideDefault) {
        this.units = units;
        this.names = names;
        this.contexts = contexts;
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
        return new Axes(units, null, null, false);
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
        return new Axes(units, null, null, true);
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
     * Constructs an axes object that contains the contexts for all axes, except for
     * the default axes.
     * 
     * @param contexts the contexts for the axes
     * @return the axes with the specified contexts
     */
    public static final Axes inContexts(UnitContext... contexts) {
        return new Axes(null, null, contexts, false);
    }

    /**
     * Constructs an axes object that contains the contexts for all axes, including
     * the default axes.
     * 
     * @param contexts the contexts for the axes
     * @return the axes with the specified contexts
     */
    public static final Axes allInContexts(UnitContext... contexts) {
        return new Axes(null, null, contexts, true);
    }

    /**
     * Constructs an axes object that contains the names for all axes, except for
     * the default axes.
     * 
     * @param names the names for the axes
     * @return the axes with the specified names
     */
    public static final Axes withNames(String... names) {
        return new Axes(null, names, null, false);
    }

    /**
     * Constructs an axes object that contains the names for all axes, including the
     * default axes.
     * 
     * @param names the names for the axes
     * @return the axes with the specified names
     */
    public static final Axes allWithNames(String... names) {
        return new Axes(null, names, null, true);
    }

    /**
     * Checks whether the current Axes contain information about the units to use on
     * the axes.
     * 
     * @return true if unit information is available
     */
    private boolean hasUnitInfo() {
        return units != null;
    }

    /**
     * Checks whether the current Axes contain information about the axis names.
     * 
     * @return true if axis name information is available
     */
    private boolean hasNameInfo() {
        return names != null;
    }

    /**
     * Checks whether the current Axes contain information about the contexts.
     * 
     * @return true if context information is available
     */
    private boolean hasContextInfo() {
        return contexts != null;
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

        var unitInfos = new ArrayList<Axes>();
        var nameInfos = new ArrayList<Axes>();
        var contextInfos = new ArrayList<Axes>();
        for (Object arg : args) {
            if (arg instanceof Axes a) {
                if (a.hasUnitInfo()) unitInfos.add(a);
                if (a.hasNameInfo()) nameInfos.add(a);
                if (a.hasContextInfo()) contextInfos.add(a);
            }
        }

        unitInfos.forEach(a -> updateUnits(a, axes));
        contextInfos.forEach(a -> updateContexts(a, axes));
        nameInfos.forEach(a -> updateNames(a, axes));

        // directly provided Axis should take precedence
        for (Object arg : args) {
            if (arg instanceof Axis axis) {
                axes.removeIf(val -> val.dimension() == axis.dimension());
                axes.add(axis);
            }
        }

        return axes;
    }

    /**
     * Updates the {@code axisSet} with the unit information from the {@code axes}.
     * 
     * @param axes    the axes containing unit infos
     * @param axisSet the axis set to update
     */
    private static void updateUnits(Axes axes, TreeSet<Axis> axisSet) {
        AtomicInteger counter = new AtomicInteger();
        BiFunction<Axis, Axes, Axis> nac = (axis, laxes) -> {
            int i = counter.getAndIncrement();
            if (i > laxes.units.length - 1) {
                i = 0;
            }
            return new Axis(axis.dimension(), laxes.units[i], axis.name());
        };
        update(axes, axisSet, nac);
    }

    /**
     * Updates the {@code axisSet} with the context information from the {@code axes}.
     * 
     * @param axes    the axes containing context infos
     * @param axisSet the axis set to update
     */
    private static void updateContexts(Axes axes, TreeSet<Axis> axisSet) {
        AtomicInteger counter = new AtomicInteger();
        BiFunction<Axis, Axes, Axis> nac = (axis, laxes) -> {
            int i = counter.getAndIncrement();
            if (i > laxes.contexts.length - 1) {
                i = 0;
            }
            return new Axis(axis.dimension(), axis.unit(), tryToDetermineName(axis.unit(), laxes.contexts[i]));
        };

        update(axes, axisSet, nac);
    }

    /**
     * Updates the {@code axisSet} with the names information from the {@code axes}.
     * 
     * @param axes    the axes containing names infos
     * @param axisSet the axis set to update
     */
    private static void updateNames(Axes axes, TreeSet<Axis> axisSet) {
        AtomicInteger counter = new AtomicInteger();
        BiFunction<Axis, Axes, Axis> nac = (axis, laxes) -> {
            int i = counter.getAndIncrement();
            if (i > laxes.names.length - 1) {
                i = 0;
            }
            return new Axis(axis.dimension(), axis.unit(), laxes.names[i]);
        };
        update(axes, axisSet, nac);
    }

    /**
     * Tries to infer the name of an axis by checking the given unit in the
     * specified contexts, by checking all values of {@link UnitContextMatch}, going
     * from strict to less strict matches.
     * 
     * @param unit     the unit
     * @param contexts the contexts
     * @return the inferred axis name, or an empty String if no unique matching name
     *         was found
     */
    private static String tryToDetermineName(Unit unit, UnitContext... contexts) {
        // try to determine the name of the unit
        NavigableSet<String> names = EMPTY_AXIS_NAMES;
        for (var matcher : UnitContextMatch.values()) {
            names = Units.inContext(unit, matcher, contexts);
            if (names.size() == 1) break;
        }

        if (names.size() == 1) {
            return names.first();
        }

        // so we might have no value or multiple values, implying that we cannot
        // determine a name uniquely -> no name used at all
        return "";
    }

    /**
     * Updates the axis set via the {@code axisUpdater} and the {@code axes}.
     * 
     * @param axes        the axes containing info for the update of the axis set
     * @param axisSet     the axis set to update
     * @param axisUpdater the function used to update an axis
     */
    private static void update(Axes axes, TreeSet<Axis> axisSet, BiFunction<Axis, Axes, Axis> axisUpdater) {
        if (axes.overrideDefault) {
            // we override everything in order
            var axesToChange = new TreeSet<>(axisSet);
            for (Axis axis : axesToChange) {
                axisSet.removeIf(val -> val.dimension() == axis.dimension());
                axisSet.add(axisUpdater.apply(axis, axes));
            }
        } else {
            // try to get the first non-default axis 
            Axis startAxis = Axis.fromSet(axisSet, 0);

            // ensure that it is not the default axis that was used as a fallback
            if (startAxis.dimension() != Axes.DEFAULT_DIMENSION) {
                var axesToChange = new TreeSet<>(axisSet.tailSet(startAxis, true));
                for (Axis axis : axesToChange) {
                    axisSet.removeIf(val -> val.dimension() == axis.dimension());
                    axisSet.add(axisUpdater.apply(axis, axes));
                }
            }
        }
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
     * @throws NullPointerException if {@code defaultAxes}, {@code argCheck} or
     *                              {@code args} is null
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
     * @throws NullPointerException if {@code axes} is null
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
