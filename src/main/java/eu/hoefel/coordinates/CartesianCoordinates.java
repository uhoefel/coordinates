package eu.hoefel.coordinates;

import java.util.Collections;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Objects;
import java.util.TreeMap;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import eu.hoefel.coordinates.axes.Axes;
import eu.hoefel.coordinates.axes.Axis;
import eu.hoefel.coordinates.tensors.TensorTransformation;
import eu.hoefel.unit.Unit;
import eu.hoefel.unit.si.SiBaseUnit;
import eu.hoefel.utils.Maths;

/**
 * The Cartesian coordinate system, named in honor of René Descartes (with his
 * surname being Cartesius in its latinized form), connected Euclidean geometry
 * and algebra for the first time in the 17<sup>th</sup> century. It is an
 * orthogonal coordinate system.
 * <p>
 * An example for a 2D Cartesian coordinate system can be seen here:<br>
 * <img src="doc-files/cartesian.svg" alt="cartesian coordinate system"><br>
 * The two axes are denoted <i>x</i> and <i>y</i>. The origin (0,0) is shown in
 * purple. Further example points, and how the values of their components relate
 * to the axes are shown in red, green and blue.
 * 
 * <p>
 * The Cartesian coordinate system is used widely throughout engineering,
 * physics, mathematics and many more branches of basic and applied science.
 * <p>
 * To get a new instance of a default Cartesian coordinate system with
 * {@link SiBaseUnit#METER} as the unit on its axes you can just do:
 * <p>
 * <code>var cartesian = new CartesianCoordinates();</code>
 * <p>
 * Note that it can in principle handle a dimensionality up to the maximum array
 * size.
 * 
 * @param axes the axes defining the coordinate system
 *
 * @author Udo Hoefel
 */
@CoordinateSystemSymbols({"cartesian", "cart"})
public final record CartesianCoordinates(NavigableSet<Axis> axes) implements CoordinateSystem {

	// Note that we do not use validatePosition here, as it is extremely unlikely
	// that users want an array with such a humongous size (and if they do they will
	// run into other issues anyway) and the check is probably not worth the cost.

	/** The default axes. */
	public static final NavigableSet<Axis> DEFAULT_AXES = Axes.of(new Axis(SiBaseUnit.METER));

	/**
	 * The consumer useful for checking the arguments in
 	 * {@link #CartesianCoordinates(Object...)}.
	 */
	private static final Consumer<Object[]> ARG_CHECK = args -> Axes.DEFAULT_ARG_CHECK.accept("Cartesian", args);

	/** Constructs a new Cartesian coordinate system. */
	public CartesianCoordinates {
		Objects.requireNonNull(axes, "Axes may not be null. "
				+ "Use the DEFAULT_AXES or the empty constructor to get a reasonable default.");
	}

	/**
	 * Constructs a new Cartesian coordinate system.
	 * 
	 * @param args the arguments, must be either {@link Axes} or {@link Axis}, which
	 *			   will take precedence over the {@link #DEFAULT_AXES} if given. If
	 *			   no arguments are given, the default axes will be used.
	 */
	public CartesianCoordinates(Object... args) {
		this(Axes.fromArgs(DEFAULT_AXES, ARG_CHECK, args));
	}

	@Override
	public int dimension() {
		return Maths.MAX_ARRAY_SIZE;
	}

	@Override
	public boolean isOrthogonal() {
		return true;
	}
	
	@Override
	public boolean isBasic() {
		return true;
	}

	@Override
	public NavigableMap<Integer, Unit> toBaseUnits() {
		var map = axes.stream().collect(Collectors.toMap(Axis::dimension, Axis::unit));
		return Collections.unmodifiableNavigableMap(new TreeMap<>(map));
	}

	@Override
	public double[] toBasePoint(double[] position) {
		return position.clone();
	}

	@Override
	public double[] fromBasePoint(double[] position) {
		return position.clone();
	}

	@Override
	public double[] toBaseVector(double[] position, double[] vector) {
		return vector.clone();
	}

	@Override
	public double[] fromBaseVector(double[] position, double[] vector) {
		return vector.clone();
	}

	@Override
	public double metricCoefficient(double[] position, TensorTransformation behavior, int i, int j) {
		if (i < 0 || j < 0) {
			throw new IllegalArgumentException("Metric coefficient not available for i=%d, j=%d (too low dimension)"
					.formatted(i, j));
		}
		return i == j ? 1 : 0;
	}

	@Override
	public double[][] metricTensor(double[] position, TensorTransformation behavior) {
		return Maths.eye(position.length);
	}

	@Override
	public double jacobianDeterminant(double[] position) {
		return 1;
	}

	@Override
	public double christoffelSymbol1stKind(double[] position, int i, int j, int k) {
		int dim = position.length;
		if (i < 0 || j < 0 || k < 0 || i >= dim || j >= dim || k >= dim) {
			throw new IllegalArgumentException(
					"i, j and k may not be <0 or exceed %d, but they were i=%d, j=%d and k=%d"
							.formatted(dim, i, j, k));
		}

		return 0;
	}

	@Override
	public double christoffelSymbol2ndKind(double[] position, int m, int i, int j) {
		int dim = position.length;
		if (m < 0 || i < 0 || j < 0 || m >= dim || i >= dim || j >= dim) {
			throw new IllegalArgumentException(
					"m, i and j may not be <0 or exceed %d, but they were m=%d, i=%d and j=%d"
							.formatted(dim, m, i, j));
		}

		return 0;
	}

	@Override
	public double riemannTensor(double[] position, int mu, int nu, int rho, int sigma) {
		int dim = position.length;
		if (mu < 0 || nu < 0 || rho < 0 || sigma < 0 || mu >= dim || nu >= dim || rho >= dim || sigma >= dim) {
			throw new IllegalArgumentException(
					"mu, nu, rho and sigma may not be <0 or exceed %d, but they were mu=%d, nu=%d, rho=%d and sigma=%d"
							.formatted(dim, mu, nu, rho, sigma));
		}

		return 0;
	}

	@Override
	public boolean isFlat() {
		return true;
	}

	@Override
	public double ricciTensor(double[] position, int mu, int nu) {
		int dim = position.length;
		if (mu < 0 || nu < 0 || mu >= dim || nu >= dim) {
			throw new IllegalArgumentException(
					"mu and nu may not be <0 or exceed %d, but they were mu=%d and nu=%d".formatted(dim, mu, nu));
		}

		return 0;
	}

	@Override
	public double ricciScalar(double[] position) {
		return 0;
	}
}
