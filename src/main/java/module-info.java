/**
 * This module provides support for programmatic (generic) coordinate systems,
 * such as {@link CartesianCoordinateSystem}.
 * <p>
 * The {@link eu.hoefel.coordinates} is the package containing the base classes
 * needed for the handling of coordinate systems, as well as some basic
 * implementations, e.g. for Cartesian, cylindrical and spherical coordinate
 * systems. It requires both the {@link eu.hoefel.coordinates.axes} and
 * {@link eu.hoefel.coordinates.tensors} packages, as well as the external
 * dependency {@link eu.hoefel.units} to manage the units on the axes.
 * <p>
 * The {@link eu.hoefel.coordinates.axes} is the package containing axes related
 * classes necessary for each coordinate system. It requires the external
 * dependency {@link eu.hoefel.units} to manage the units on the axes.
 * 
 * The {@link eu.hoefel.coordinates.tensors} is the package containing classes
 * for basic tensor transformations.
 * 
 * @author Udo Hoefel
 */
module eu.hoefel.coordinates {
	exports eu.hoefel.coordinates;
	exports eu.hoefel.coordinates.axes;
	exports eu.hoefel.coordinates.tensors;

	// ugh...
	opens eu.hoefel.coordinates to org.junit.platform.commons;
	opens eu.hoefel.coordinates.axes to org.junit.platform.commons;
	opens eu.hoefel.coordinates.tensors to org.junit.platform.commons;

	requires java.logging;

	requires org.junit.jupiter.api;
	
	requires eu.hoefel.utils;
	requires transitive eu.hoefel.unit;
}