package eu.hoefel.coordinates.axes;

import java.util.Iterator;
import java.util.NavigableSet;
import java.util.Objects;

import eu.hoefel.coordinates.CoordinateSystem;
import eu.hoefel.unit.Unit;

/**
 * Describes an axis within a {@link CoordinateSystem}.
 * 
 * @param dimension the dimension to which this axis belongs. Note that the
 *                  first dimension should be at 0, not 1. Has to be &ge;0,
 *                  except for {@link Axes#DEFAULT_DIMENSION}.
 * @param unit      the unit, e.g. the unit of "m" or "s m^-2". May not be null.
 * @param name      the name describing this axis. This may be "x", or "B_z" or
 *                  whatever else describes the axis. <em>It should not contain
 *                  unit information!</em> Hence, names like "B_z in T" should
 *                  be avoided. Potential plotting routines should build up the
 *                  "in T" part from the unit information provided by
 *                  {@code unit}. May be empty, but not null.
 * @author Udo Hoefel
 */
public final record Axis(int dimension, Unit unit, String name) {

    /**
     * Creates a new Axis.
     * 
     * @param dimension the dimension to which this axis belongs. Note that the
     *                  first dimension should be at 0, not 1. Has to be &ge;0,
     *                  except for {@link Axes#DEFAULT_DIMENSION}.
     * @param unit      the unit, e.g. the unit of "m" or "s m^-2". May not be null.
     * @param name      the name describing this axis. This may be "x", or "B_z" or
     *                  whatever else describes the axis. <em>It should not contain
     *                  unit information!</em> Hence, names like "B_z in T" should
     *                  be avoided. Potential plotting routines should build up the
     *                  "in T" part from the unit information provided by
     *                  {@code unit}. May be empty, but not null.
     * @throws NullPointerException     if {@code unit} or {@code name} are null
     * @throws IllegalArgumentException if {@code dimension} is negative and not the
     *                                  {@link Axes#DEFAULT_DIMENSION default
     *                                  dimension}.
     */
    public Axis {
        Objects.requireNonNull(unit);
        Objects.requireNonNull(name);

        if (dimension < 0 && dimension != Axes.DEFAULT_DIMENSION) {
            throw new IllegalArgumentException("Negative dimensions are not supported!");
        }
    }

    /**
     * Constructs a new axis with the specified properties. Note that no operation
     * will be performed on the axis.
     * 
     * @param dimension  the dimension to which this axis belongs. Note that the
     *                   first dimension should be at 0, not 1.
     * @param units      the units, e.g. "m" or "s m^-2"
     * @param extraUnits additional units required for the parsing of {@code unit}
     */
    public Axis(int dimension, String units, Unit... extraUnits) {
        this(dimension, Unit.of(units, extraUnits), "");
    }

    /**
     * Gets the axis of the specified dimension from the given set. Returns the
     * default axis if a default axis is available but no specific axis.
     * 
     * @param axes      the axes to search in
     * @param dimension the dimension
     * @return the axis corresponding to the specified dimension
     * @throws IllegalArgumentException if neither an axis with the exact dimension
     *                                  was found nor a default dimension that could
     *                                  be used as a fallback was available
     */
    public static Axis fromSet(NavigableSet<Axis> axes, int dimension) {
        Objects.requireNonNull(axes);

        int size = axes.size();

        // let's see if we have a default dimension
        boolean hasDefault = axes.first().dimension() == Axes.DEFAULT_DIMENSION;
        int index = hasDefault ? dimension + 1 : dimension;

        // iterate either form the top down or from the bottom up
        boolean descend = index >= size / 2;
        Iterator<Axis> itr = descend ? axes.descendingIterator() : axes.iterator();
        int maxSteps = descend ? size - index : index + 1;

        for (int i = 0; i < maxSteps && itr.hasNext(); i++) {
            Axis axis = itr.next();
            if (axis.dimension() == dimension) return axis;
        }

        // So we don't have an exact match for the desired dimension. Is there a default
        // dimension to fall back to?
        if (hasDefault) return axes.first();

        throw new IllegalArgumentException(
                "No axis with dimension " + dimension + " known. Also, no default dimension is available to use.");
    }
}
