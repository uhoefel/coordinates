package eu.hoefel.coordinates;

import java.util.Collections;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Objects;
import java.util.TreeMap;
import java.util.function.Consumer;
import java.util.function.Function;

import eu.hoefel.coordinates.axes.Axes;
import eu.hoefel.coordinates.axes.Axis;
import eu.hoefel.coordinates.tensors.TensorIndexType;
import eu.hoefel.coordinates.tensors.TensorTransformation;
import eu.hoefel.unit.Unit;
import eu.hoefel.unit.Units;

/**
 * 6-sphere coordinates are a 3D coordinate system. To derive them consider 3D
 * {@link CartesianCoordinates} inverted across the 2-sphere (i.e. the standard
 * 2 dimensional sphere). Their name stems from the fact that keeping one
 * dimension constant forms a lobe pointing from the origin in the direction of
 * an axis. As this means that for each axis one lobe points toward the positive
 * direction and one in the negative direction, there are six lobes in total.
 * <p>
 * Note that they are not related to the 6-sphere (i.e. the 6 dimensional
 * sphere).
 * 
 * @param axes the axes defining the coordinate system, not null
 * 
 * @author Udo Hoefel
 */
@CoordinateSystemSymbols({"6sphere", "6sph"})
public final record SixSphereCoordinates(NavigableSet<Axis> axes) implements CoordinateSystem {

	/** The default axes. */
    public static final NavigableSet<Axis> DEFAULT_AXES = Axes.of(
            new Axis(0, Unit.of("m^-1"), "u"),
            new Axis(1, Unit.of("m^-1"), "v"), 
            new Axis(2, Unit.of("m^-1"), "w"));

	/**
	 * The consumer useful for checking the arguments in
 	 * {@link #SixSphereCoordinates(Object...)}.
	 */
	private static final Consumer<Object[]> ARG_CHECK = args -> Axes.DEFAULT_ARG_CHECK.accept("6-Sphere", args);

	/** Constructs a new 6-sphere coordinate system. */
	public SixSphereCoordinates {
		Objects.requireNonNull(axes, "Axes may not be null. "
				+ "Use the DEFAULT_AXES or the empty constructor to get a reasonable default.");
	}

	/**
	 * Constructs a new 6-sphere coordinate system.
	 * 
	 * @param args the arguments, must be either {@link Axes} or {@link Axis}, which
	 *             will take precedence over the {@link #DEFAULT_AXES} if given. If
	 *             no arguments are given, the default axes will be used.
	 */
	public SixSphereCoordinates(Object... args) {
		this(Axes.fromArgs(DEFAULT_AXES, ARG_CHECK, args));
	}

	/**
	 * Validates the position, i.e. it throws an exception if a dimension of the 
	 * position is out of range.
	 * 
	 * @param position the position to validate
	 * @throw IllegalArgumentException if the assumptions about the dimensionality 
	 *		  or the valid range of any dimension of the input are violated.
	 */
	private void validatePosition(double[] position) {
		Objects.requireNonNull(position);

		if (position.length > dimension()) {
			throw new IllegalArgumentException(
					"The given dimensionality exceeds the maximum supported dimensionality (%d vs %d)"
							.formatted(position.length, dimension()));
		}
	}

	@Override
	public int dimension() {
		return 3;
	}

	@Override
	public boolean isOrthogonal() {
		return true;
	}

	@Override
	public NavigableMap<Integer, Unit> toBaseUnits() {
		NavigableMap<Integer, Unit> map = new TreeMap<>();
		map.put(0, Unit.of(Units.simplify(axis(0).unit().symbols().get(0) + "^-1")));
		map.put(1, Unit.of(Units.simplify(axis(0).unit().symbols().get(0) + "^-2 " + axis(1).unit().symbols().get(0))));
		map.put(2, Unit.of(Units.simplify(axis(0).unit().symbols().get(0) + "^-2 " + axis(2).unit().symbols().get(0))));
		return Collections.unmodifiableNavigableMap(map);
	}

	@Override
	public double[] toBasePoint(double[] position) {
		validatePosition(position);
		return Transformations.sixSphereToCartesian(position);
	}

	@Override
	public double[] fromBasePoint(double[] position) {
		return Transformations.cartesianToSixSphere(position);
	}

	@Override
	public double[] toBaseVector(double[] position, double[] vector) {
		validatePosition(position);

		double[] vectorInBaseSys = new double[vector.length];
		vectorInBaseSys[0] = Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*(Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(position[0],2))*vector[0]+-2*position[0]*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[1]+-2*position[0]*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[2];
		vectorInBaseSys[1] = -2*position[0]*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[0]+Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*(Math.pow(position[0],2)+Math.pow(position[2],2)-1*Math.pow(position[1],2))*vector[1]+-2*position[1]*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[2];
		vectorInBaseSys[2] = -2*position[0]*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[0]+-2*position[1]*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[1]+Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*(Math.pow(position[0],2)+Math.pow(position[1],2)-1*Math.pow(position[2],2))*vector[2];
		return vectorInBaseSys;
	}

	@Override
	public double[] fromBaseVector(double[] position, double[] vector) {
		double[] vectorInCurrentSys = new double[vector.length];
		vectorInCurrentSys[0] = Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*(Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(position[0],2))*vector[0]-2*position[0]*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[1]-2*position[0]*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[2];
		vectorInCurrentSys[1] = -2*position[0]*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[0]+Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*(Math.pow(position[0],2)+Math.pow(position[2],2)-1*Math.pow(position[1],2))*vector[1]-2*position[1]*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[2];
		vectorInCurrentSys[2] = -2*position[0]*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[0]-2*position[1]*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*vector[1]+Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-2)*(Math.pow(position[0],2)+Math.pow(position[1],2)-1*Math.pow(position[2],2))*vector[2];
		return vectorInCurrentSys;
	}

	@Override
	public double metricCoefficient(double[] position, TensorTransformation behavior, int i, int j) {
		validatePosition(position);

		if (i < 0 || j < 0) {
			throw new IllegalArgumentException(("Metric coefficient not available for i=%d, j=%d "
					+ "(too low dimension, only 3 dimensions [0,1,2] are supported for sixsphere coordinates)")
					.formatted(i, j));
		} else if (i > 2 || j > 2) {
			throw new IllegalArgumentException(("Metric coefficient not available for i=%d, j=%d "
					+ "(too high dimension, only 3 dimensions [0,1,2] are supported for sixsphere coordinates)")
					.formatted(i, j));
		}

		if (behavior instanceof TensorIndexType tit) {
			return switch (tit) {
				case COVARIANT -> {
					if (i == 0 && j == 0) yield Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-4)*(Math.pow((Math.pow(position[0],2)-1*Math.pow(position[1],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[0],2)*(Math.pow(position[1],2)+Math.pow(position[2],2)));
					if (i == 1 && j == 1) yield Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-4)*(Math.pow((Math.pow(position[1],2)-1*Math.pow(position[0],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[1],2)*(Math.pow(position[0],2)+Math.pow(position[2],2)));
					if (i == 2 && j == 2) yield Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-4)*(Math.pow((Math.pow(position[2],2)-1*Math.pow(position[0],2)-1*Math.pow(position[1],2)),2)+4*Math.pow(position[2],2)*(Math.pow(position[0],2)+Math.pow(position[1],2)));
					yield 0;
				}
				case CONTRAVARIANT -> {
					if (i == 0 && j == 0) yield Math.pow((Math.pow((Math.pow(position[0],2)-1*Math.pow(position[1],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[0],2)*(Math.pow(position[1],2)+Math.pow(position[2],2))),-1)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),4);
					if (i == 1 && j == 1) yield Math.pow((Math.pow((Math.pow(position[1],2)-1*Math.pow(position[0],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[1],2)*(Math.pow(position[0],2)+Math.pow(position[2],2))),-1)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),4);
					if (i == 2 && j == 2) yield Math.pow((Math.pow((Math.pow(position[2],2)-1*Math.pow(position[0],2)-1*Math.pow(position[1],2)),2)+4*Math.pow(position[2],2)*(Math.pow(position[0],2)+Math.pow(position[1],2))),-1)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),4);
					yield 0;
				}
			};
		}

		// Fall back to metric tensor (which might fall back to this method) for mixed behavior
		Function<double[], double[][]> metricTensor = pos -> metricTensor(pos, TensorIndexType.COVARIANT);
		return TensorIndexType.COVARIANT.transform(this, metricTensor, behavior).apply(position)[i][j];
	}

	@Override
	public double[][] metricTensor(double[] position, TensorTransformation behavior) {
		validatePosition(position);

		int dim = 3;
		double[][] g = new double[dim][dim];

		// Note that we skip all elements that are zero anyways
		g[0][0] = metricCoefficient(position, behavior, 0, 0);
		g[1][1] = metricCoefficient(position, behavior, 1, 1);
		g[2][2] = metricCoefficient(position, behavior, 2, 2);
		return g;
	}

	@Override
	public double jacobianDeterminant(double[] position) {
		validatePosition(position);
		return Math.pow((Math.pow(position[0],12)+Math.pow(position[1],12)+Math.pow(position[2],12)+6*Math.pow(position[0],2)*Math.pow(position[1],10)+6*Math.pow(position[0],2)*Math.pow(position[2],10)+6*Math.pow(position[0],10)*Math.pow(position[1],2)+6*Math.pow(position[0],10)*Math.pow(position[2],2)+6*Math.pow(position[1],2)*Math.pow(position[2],10)+6*Math.pow(position[1],10)*Math.pow(position[2],2)+15*Math.pow(position[0],4)*Math.pow(position[1],8)+15*Math.pow(position[0],4)*Math.pow(position[2],8)+15*Math.pow(position[0],8)*Math.pow(position[1],4)+15*Math.pow(position[0],8)*Math.pow(position[2],4)+15*Math.pow(position[1],4)*Math.pow(position[2],8)+15*Math.pow(position[1],8)*Math.pow(position[2],4)+20*Math.pow(position[0],6)*Math.pow(position[1],6)+20*Math.pow(position[0],6)*Math.pow(position[2],6)+20*Math.pow(position[1],6)*Math.pow(position[2],6)+30*Math.pow(position[0],2)*Math.pow(position[1],2)*Math.pow(position[2],8)+30*Math.pow(position[0],2)*Math.pow(position[1],8)*Math.pow(position[2],2)+30*Math.pow(position[0],8)*Math.pow(position[1],2)*Math.pow(position[2],2)+60*Math.pow(position[0],2)*Math.pow(position[1],4)*Math.pow(position[2],6)+60*Math.pow(position[0],2)*Math.pow(position[1],6)*Math.pow(position[2],4)+60*Math.pow(position[0],4)*Math.pow(position[1],2)*Math.pow(position[2],6)+60*Math.pow(position[0],4)*Math.pow(position[1],6)*Math.pow(position[2],2)+60*Math.pow(position[0],6)*Math.pow(position[1],2)*Math.pow(position[2],4)+60*Math.pow(position[0],6)*Math.pow(position[1],4)*Math.pow(position[2],2)+90*Math.pow(position[0],4)*Math.pow(position[1],4)*Math.pow(position[2],4)),-1.0/2);
	}

	@Override
	public double christoffelSymbol1stKind(double[] position, int i, int j, int k) {
		validatePosition(position);

		int dim = position.length;
		if (i < 0 || j < 0 || k < 0 || i >= dim || j >= dim || k >= dim) {
			throw new IllegalArgumentException(
					"i, j and k may not be <0 or exceed %d, but they were i=%d, j=%d and k=%d"
							.formatted(dim, i, j, k));
		}

		if (i == 0 && j == 0 && k == 0) return -2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(position[0],2)),2)+4*Math.pow(position[0],2)*(Math.pow(position[1],2)+Math.pow(position[2],2)));
		if (i == 0 && j == 0 && k == 1) return 2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[2],2)-1*Math.pow(position[1],2)),2)+4*Math.pow(position[1],2)*(Math.pow(position[0],2)+Math.pow(position[2],2)));
		if (i == 0 && j == 0 && k == 2) return 2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[2],2)*(Math.pow(position[0],2)+Math.pow(position[1],2)));
		if (i == 0 && j == 1 && k == 0) return -2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(position[0],2)),2)+4*Math.pow(position[0],2)*(Math.pow(position[1],2)+Math.pow(position[2],2)));
		if (i == 0 && j == 1 && k == 1) return -2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[2],2)-1*Math.pow(position[1],2)),2)+4*Math.pow(position[1],2)*(Math.pow(position[0],2)+Math.pow(position[2],2)));
		if (i == 0 && j == 2 && k == 0) return -2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(position[0],2)),2)+4*Math.pow(position[0],2)*(Math.pow(position[1],2)+Math.pow(position[2],2)));
		if (i == 0 && j == 2 && k == 2) return -2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[2],2)*(Math.pow(position[0],2)+Math.pow(position[1],2)));
		if (i == 1 && j == 0 && k == 0) return -2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(position[0],2)),2)+4*Math.pow(position[0],2)*(Math.pow(position[1],2)+Math.pow(position[2],2)));
		if (i == 1 && j == 0 && k == 1) return -2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[2],2)-1*Math.pow(position[1],2)),2)+4*Math.pow(position[1],2)*(Math.pow(position[0],2)+Math.pow(position[2],2)));
		if (i == 1 && j == 1 && k == 0) return 2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(position[0],2)),2)+4*Math.pow(position[0],2)*(Math.pow(position[1],2)+Math.pow(position[2],2)));
		if (i == 1 && j == 1 && k == 1) return -2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[2],2)-1*Math.pow(position[1],2)),2)+4*Math.pow(position[1],2)*(Math.pow(position[0],2)+Math.pow(position[2],2)));
		if (i == 1 && j == 1 && k == 2) return 2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[2],2)*(Math.pow(position[0],2)+Math.pow(position[1],2)));
		if (i == 1 && j == 2 && k == 1) return -2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[2],2)-1*Math.pow(position[1],2)),2)+4*Math.pow(position[1],2)*(Math.pow(position[0],2)+Math.pow(position[2],2)));
		if (i == 1 && j == 2 && k == 2) return -2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[2],2)*(Math.pow(position[0],2)+Math.pow(position[1],2)));
		if (i == 2 && j == 0 && k == 0) return -2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(position[0],2)),2)+4*Math.pow(position[0],2)*(Math.pow(position[1],2)+Math.pow(position[2],2)));
		if (i == 2 && j == 0 && k == 2) return -2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[2],2)*(Math.pow(position[0],2)+Math.pow(position[1],2)));
		if (i == 2 && j == 1 && k == 1) return -2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[2],2)-1*Math.pow(position[1],2)),2)+4*Math.pow(position[1],2)*(Math.pow(position[0],2)+Math.pow(position[2],2)));
		if (i == 2 && j == 1 && k == 2) return -2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[2],2)*(Math.pow(position[0],2)+Math.pow(position[1],2)));
		if (i == 2 && j == 2 && k == 0) return 2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(position[0],2)),2)+4*Math.pow(position[0],2)*(Math.pow(position[1],2)+Math.pow(position[2],2)));
		if (i == 2 && j == 2 && k == 1) return 2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[2],2)-1*Math.pow(position[1],2)),2)+4*Math.pow(position[1],2)*(Math.pow(position[0],2)+Math.pow(position[2],2)));
		if (i == 2 && j == 2 && k == 2) return -2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-5)*(Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)-1*Math.pow(position[2],2)),2)+4*Math.pow(position[2],2)*(Math.pow(position[0],2)+Math.pow(position[1],2)));
		return 0;
	}

	@Override
	public double christoffelSymbol2ndKind(double[] position, int m, int i, int j) {
		validatePosition(position);

		int dim = position.length;
		if (m < 0 || i < 0 || j < 0 || m >= dim || i >= dim || j >= dim) {
			throw new IllegalArgumentException(
					"m, i and j may not be <0 or exceed %d, but they were m=%d, i=%d and j=%d"
							.formatted(dim, m, i, j));
		}

		if (m == 0 && i == 0 && j == 0) return -2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 0 && i == 0 && j == 1) return -2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 0 && i == 0 && j == 2) return -2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 0 && i == 1 && j == 0) return -2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 0 && i == 1 && j == 1) return 2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 0 && i == 2 && j == 0) return -2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 0 && i == 2 && j == 2) return 2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 1 && i == 0 && j == 0) return 2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 1 && i == 0 && j == 1) return -2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 1 && i == 1 && j == 0) return -2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 1 && i == 1 && j == 1) return -2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 1 && i == 1 && j == 2) return -2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 1 && i == 2 && j == 1) return -2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 1 && i == 2 && j == 2) return 2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 2 && i == 0 && j == 0) return 2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 2 && i == 0 && j == 2) return -2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 2 && i == 1 && j == 1) return 2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 2 && i == 1 && j == 2) return -2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 2 && i == 2 && j == 0) return -2*position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 2 && i == 2 && j == 1) return -2*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		if (m == 2 && i == 2 && j == 2) return -2*position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1);
		return 0;
	}

	@Override
	public double riemannTensor(double[] position, int mu, int nu, int rho, int sigma) {
		validatePosition(position);

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
		validatePosition(position);

		int dim = position.length;
		if (mu < 0 || nu < 0 || mu >= dim || nu >= dim) {
			throw new IllegalArgumentException(
					"mu and nu may not be <0 or exceed %d, but they were mu=%d and nu=%d".formatted(dim, mu, nu));
		}

		return 0;
	}

	@Override
	public double ricciScalar(double[] position) {
		validatePosition(position);
		return 0;
	}

	@Override
	public <T> double magnitude(double[] position, TensorTransformation transformation, Function<double[], T> tensorfield) {
		validatePosition(position);
		return CoordinateSystem.super.magnitude(position, transformation, tensorfield);
	}

	@Override
	public double ds(double[] position, int i, double dui) {
		validatePosition(position);
		return CoordinateSystem.super.ds(position, i, dui);
	}

	@Override
	public double dA(double[] position, int i, int j, double dui, double duj) {
		validatePosition(position);
		return CoordinateSystem.super.dA(position, i, j, dui, duj);
	}

	@Override
	public double dV(double[] position, double... du) {
		validatePosition(position);
		return CoordinateSystem.super.dV(position, du);
	}

	@Override
	public double dot(double[] position, TensorIndexType behavior, double[] v1, double[] v2) {
		validatePosition(position);
		return CoordinateSystem.super.dot(position, behavior, v1, v2);
	}

	@Override
	public double[] cross(double[] position, TensorIndexType behavior, double[] v1, double[] v2, double[]... vn) {
		validatePosition(position);
		return CoordinateSystem.super.cross(position, behavior, v1, v2, vn);
	}

	@Override
	public <T> T div(double[] position, TensorTransformation componentBehavior, Function<double[], T[]> field) {
		validatePosition(position);
		return CoordinateSystem.super.div(position, componentBehavior, field);
	}

	@Override
	public <T> T[] grad(double[] position, TensorTransformation componentBehavior, Function<double[], T> field) {
		validatePosition(position);
		return CoordinateSystem.super.grad(position, componentBehavior, field);
	}

	@Override
	public double[] curl(double[] position, TensorTransformation componentBehavior, Function<double[], double[]> vectorfield) {
		validatePosition(position);
		return CoordinateSystem.super.curl(position, componentBehavior, vectorfield);
	}
}
