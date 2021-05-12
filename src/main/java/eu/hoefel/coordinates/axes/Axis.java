package eu.hoefel.coordinates.axes;

import java.util.NavigableSet;

import eu.hoefel.coordinates.CoordinateSystem;
import eu.hoefel.unit.Unit;
import eu.hoefel.unit.Units;
import eu.hoefel.unit.si.SiBaseUnit;

/**
 * Describes an axis within a {@link CoordinateSystem}.
 * 
 * @param dimension  the dimension to which this axis belongs. Note that the
 *                   first dimension should be at 0, not 1.
 * @param unit       the unit, e.g. the unit of "m" or "s m^-2"
 * 
 * @author Udo Hoefel
 */
public final record Axis(int dimension, Unit unit) {

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
		this(dimension, Unit.of(units, extraUnits));
	}

	/**
	 * Constructs a new axis with the specified properties. Note that no operation
	 * will be performed on the axis and the dimension this axis is assigned to is
	 * the {@link Axes#DEFAULT_DIMENSION}.
	 * 
	 * @param units      the units, e.g. "m" or "s m^-2"
	 * @param extraUnits additional units required for the parsing of {@code unit}
	 */
	public Axis(String units, Unit... extraUnits) {
		this(Axes.DEFAULT_DIMENSION, units, extraUnits);
	}

	/**
	 * Constructs a new axis with the specified properties. Note that no operation
	 * will be performed on the axis and the dimension this axis is assigned to is
	 * the {@link Axes#DEFAULT_DIMENSION}.
	 * 
	 * @param unit the unit, e.g. {@link SiBaseUnit#METER}
	 */
	public Axis(Unit unit) {
		this(Axes.DEFAULT_DIMENSION, unit);
	}

	/**
	 * Gets an axis of which only the dimension is guaranteed to be useful, as it
	 * corresponds to the given value. The main intention for this method is to have
	 * a clean way to get an axis from the {@link CoordinateSystem#axes()}.
	 * 
	 * @param dimension the dimension the axis should correspond to
	 * @return an axis placeholder with the specified dimension set correctly
	 */
	public static final Axis forDim(int dimension) {
		return new Axis(dimension, Units.EMPTY_UNIT);
	}

	/**
	 * Gets the axis of the specified dimension from the given set. Returns the
	 * default axis if a default axis is available but no specific axis.
	 * 
	 * @param axes      the axes to search in
	 * @param dimension the dimension
	 * @return the axis corresponding to the specified dimension
	 */
	public static final Axis fromSet(NavigableSet<Axis> axes, int dimension) {
		Axis axis = axes.floor(Axis.forDim(dimension));
		if (axis == null || axis.dimension() != dimension) {
			axis = axes.floor(Axis.forDim(Axes.DEFAULT_DIMENSION));
			if (axis == null || axis.dimension() != Axes.DEFAULT_DIMENSION) {
				throw new IllegalArgumentException("No axis with dimension " + dimension + " known. Also, no default dimension is available to use.");
			}
		}
		return axis;
	}
}
