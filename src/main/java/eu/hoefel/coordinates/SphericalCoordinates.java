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
import eu.hoefel.unit.si.SiBaseUnit;
import eu.hoefel.unit.si.SiDerivedUnit;
import eu.hoefel.utils.Maths;

/**
 * The spherical coordinate system is an extension of the {@link PolarCoordinates} to 3D (and, in fact, <i>n</i>-D).
 * The radius <i>ρ</i> and multiple angles <i>φ</i><sub>1</sub>…<i>φ</i><sub><i>n</i></sub> are used to determine the position of a point.
 * <p>
 * They are defined here as follows:<br>
 * {@code x_1 = ρ cos(φ_1)}<br>
 * {@code x_2 = ρ sin(φ_1)cos(φ_2)}<br>
 * {@code x_3 = ρ sin(φ_1)sin(φ_2)cos(φ_3)}<br>
 * <pre>
 *     ⋮
 * </pre>
 * {@code x_n = ρ sin(φ_1) … sin(φ_(n-1))}
 * <p>
 * to allow a straightforward extension to <i>n</i>-D.
 * For the use with this class they should be at least 3D.
 * <p>
 * Spherical coordinates are widely used, e.g. in astronomy.
 * 
 * @param axes the axes defining the coordinate system, not null. Note that the
 *             units of all dimensions, except the first, need to be effectively
 *             dimensionless.
 * 
 * @author Udo Hoefel
 * 
 * @see <a href="https://doi.org/10.2307/2308932">A Derivation of n-Dimensional Spherical Coordinates</a> 
 */
@CoordinateSystemSymbols({"spherical", "sph"})
public final record SphericalCoordinates(NavigableSet<Axis> axes) implements CoordinateSystem {

	/** The default axes. */
	public static final NavigableSet<Axis> DEFAULT_AXES = Axes.of(
			new Axis(SiDerivedUnit.RADIAN),
			new Axis(0, SiBaseUnit.METER));

	/**
	 * The consumer useful for checking the arguments in
	 * {@link #SphericalCoordinates(Object...)}.
	 */
	private static final Consumer<Object[]> ARG_CHECK = args -> Axes.DEFAULT_ARG_CHECK.accept("Spherical", args);

	/** Constructs a new spherical coordinate system. */
	public SphericalCoordinates {
		Objects.requireNonNull(axes, "Axes may not be null. "
				+ "Use the DEFAULT_AXES or the empty constructor to get a reasonable default.");

		for (Axis axis : axes) {
			if (axis.dimension() == 0) continue;
			
			if (!Units.convertible(axis.unit(), Units.EMPTY_UNIT)) {
				throw new IllegalArgumentException("The unit of dimension %d (%s) needs to be effectively dimensionless."
						.formatted(axis.dimension(), axis.unit()));
			}
		}
	}

	/**
	 * Constructs a new spherical coordinate system.
	 * 
	 * @param args the arguments, must be either {@link Axes} or {@link Axis}, which
	 *			   will take precedence over the {@link #DEFAULT_AXES} if given. If
	 *			   no arguments are given, the default axes will be used.
	 */
	public SphericalCoordinates(Object... args) {
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
	public NavigableMap<Integer, Unit> toBaseUnits() {
		NavigableMap<Integer, Unit> map = new TreeMap<>();
		// the radius alone determines the axes in the corresponding base coordinate system
		map.put(Axes.DEFAULT_DIMENSION, axis(0).unit());
		return Collections.unmodifiableNavigableMap(new TreeMap<>(map));
	}

	@Override
	public double[] toBasePoint(double[] position) {
		if (position.length == 3) {
			double[] pointInBase = new double[3];
			pointInBase[0] = position[0]*Math.cos(position[1]);
			pointInBase[1] = position[0]*Math.cos(position[2])*Math.sin(position[1]);
			pointInBase[2] = position[0]*Math.sin(position[1])*Math.sin(position[2]);
			return pointInBase;
		} else if (position.length == 4) {
			double[] pointInBase = new double[4];
			pointInBase[0] = position[0]*Math.cos(position[1]);
			pointInBase[1] = position[0]*Math.cos(position[2])*Math.sin(position[1]);
			pointInBase[2] = position[0]*Math.cos(position[3])*Math.sin(position[1])*Math.sin(position[2]);
			pointInBase[3] = position[0]*Math.sin(position[1])*Math.sin(position[2])*Math.sin(position[3]);
			return pointInBase;
		}
		
		return Transformations.sphericalToCartesian(position);
	}

	@Override
	public double[] fromBasePoint(double[] position) {
		if (position.length == 3) {
			double[] pointInCurrentSystem = new double[3];
			pointInCurrentSystem[0] = Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),0.5);
			pointInCurrentSystem[1] = Math.atan2(Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)),0.5),position[0]);
			pointInCurrentSystem[2] = Math.atan2(position[2],position[1]);
			return pointInCurrentSystem;
		} else if (position.length == 4) {
			double[] pointInCurrentSystem = new double[4];
			pointInCurrentSystem[0] = Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),0.5);
			pointInCurrentSystem[1] = Math.atan2(Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),0.5),position[0]);
			pointInCurrentSystem[2] = Math.atan2(Math.pow((Math.pow(position[2],2)+Math.pow(position[3],2)),0.5),position[1]);
			pointInCurrentSystem[3] = Math.atan2(position[3],position[2]);
			return pointInCurrentSystem;
		}
		
		return Transformations.cartesianToSpherical(position);
	}

	@Override
	public double[] toBaseVector(double[] position, double[] vector) {
		if (position.length == 3) {
			double[] vectorInBaseSys = new double[vector.length];
			vectorInBaseSys[0] = Math.cos(position[1])*vector[0]+Math.cos(position[2])*Math.sin(position[1])*vector[1]+Math.sin(position[1])*Math.sin(position[2])*vector[2];
			vectorInBaseSys[1] = -1*position[0]*Math.sin(position[1])*vector[0]+position[0]*Math.cos(position[1])*Math.cos(position[2])*vector[1]+position[0]*Math.cos(position[1])*Math.sin(position[2])*vector[2];
			vectorInBaseSys[2] = -1*position[0]*Math.sin(position[1])*Math.sin(position[2])*vector[1]+position[0]*Math.cos(position[2])*Math.sin(position[1])*vector[2];
			return vectorInBaseSys;
		} else if (position.length == 4) {
			double[] vectorInBaseSys = new double[vector.length];
			vectorInBaseSys[0] = Math.cos(position[1])*vector[0]+Math.cos(position[2])*Math.sin(position[1])*vector[1]+Math.cos(position[3])*Math.sin(position[1])*Math.sin(position[2])*vector[2]+Math.sin(position[1])*Math.sin(position[2])*Math.sin(position[3])*vector[3];
			vectorInBaseSys[1] = -1*position[0]*Math.sin(position[1])*vector[0]+position[0]*Math.cos(position[1])*Math.cos(position[2])*vector[1]+position[0]*Math.cos(position[1])*Math.cos(position[3])*Math.sin(position[2])*vector[2]+position[0]*Math.cos(position[1])*Math.sin(position[2])*Math.sin(position[3])*vector[3];
			vectorInBaseSys[2] = -1*position[0]*Math.sin(position[1])*Math.sin(position[2])*vector[1]+position[0]*Math.cos(position[2])*Math.cos(position[3])*Math.sin(position[1])*vector[2]+position[0]*Math.cos(position[2])*Math.sin(position[1])*Math.sin(position[3])*vector[3];
			vectorInBaseSys[3] = -1*position[0]*Math.sin(position[1])*Math.sin(position[2])*Math.sin(position[3])*vector[2]+position[0]*Math.cos(position[3])*Math.sin(position[1])*Math.sin(position[2])*vector[3];
			return vectorInBaseSys;
		}
		return CoordinateSystem.super.toBaseVector(position, vector);
	}

	@Override
	public double[] fromBaseVector(double[] position, double[] vector) {
		if (position.length == 3) {
			double[] vectorInCurrentSys = new double[vector.length];
			vectorInCurrentSys[0] = position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1.0/2)*vector[0]-1*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)),0.5)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1)*vector[1];
			vectorInCurrentSys[1] = position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1.0/2)*vector[0]+position[0]*position[1]*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)),-1.0/2)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1)*vector[1]-1*position[2]*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)),-1)*vector[2];
			vectorInCurrentSys[2] = position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1.0/2)*vector[0]+position[0]*position[2]*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)),-1.0/2)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)),-1)*vector[1]+position[1]*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)),-1)*vector[2];
			return vectorInCurrentSys;
		} else if (position.length == 4) {
			double[] vectorInCurrentSys = new double[vector.length];
			vectorInCurrentSys[0] = position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1.0/2)*vector[0]-1*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),0.5)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1)*vector[1];
			vectorInCurrentSys[1] = position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1.0/2)*vector[0]+position[0]*position[1]*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1.0/2)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1)*vector[1]-1*Math.pow((Math.pow(position[2],2)+Math.pow(position[3],2)),0.5)*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1)*vector[2];
			vectorInCurrentSys[2] = position[2]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1.0/2)*vector[0]+position[0]*position[2]*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1.0/2)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1)*vector[1]+position[1]*position[2]*Math.pow((Math.pow(position[2],2)+Math.pow(position[3],2)),-1.0/2)*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1)*vector[2]-1*position[3]*Math.pow((Math.pow(position[2],2)+Math.pow(position[3],2)),-1)*vector[3];
			vectorInCurrentSys[3] = position[3]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1.0/2)*vector[0]+position[0]*position[3]*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1.0/2)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1)*vector[1]+position[1]*position[3]*Math.pow((Math.pow(position[2],2)+Math.pow(position[3],2)),-1.0/2)*Math.pow((Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(position[3],2)),-1)*vector[2]+position[2]*Math.pow((Math.pow(position[2],2)+Math.pow(position[3],2)),-1)*vector[3];
			return vectorInCurrentSys;
		}
		return CoordinateSystem.super.fromBaseVector(position, vector);
	}

	@Override
	public double metricCoefficient(double[] position, TensorTransformation behavior, int i, int j) {
		if (i < 0 || j < 0) {
			throw new IllegalArgumentException("Metric coefficient not available for i=%d, j=%d (too low dimension)"
					.formatted(i, j));
		}
		
		if (behavior instanceof TensorIndexType tit) {
			return switch (tit) {
				case COVARIANT -> {
					if (i == 0 && j == 0) yield 1;
					if (i == 1 && j == 1) yield Math.pow(position[0],2);
					if (i == 2 && j == 2) yield Math.pow(position[0],2)*Math.pow(Math.sin(position[1]),2);
					if (i == 3 && j == 3) yield Math.pow(position[0],2)*Math.pow(Math.sin(position[1]),2)*Math.pow(Math.sin(position[2]),2);
					if (i > 3 && i == j) yield CoordinateSystem.super.metricCoefficient(position, behavior, i, j);
					yield 0;
				}
				case CONTRAVARIANT -> {
					if (i == 0 && j == 0) yield 1;
					if (i == 1 && j == 1) yield Math.pow(position[0],-2);
					if (i == 2 && j == 2) yield Math.pow(position[0],-2)*Math.pow(Math.sin(position[1]),-2);
					if (i == 3 && j == 3) yield Math.pow(position[0],-2)*Math.pow(Math.sin(position[1]),-2)*Math.pow(Math.sin(position[2]),-2);
					if (i > 3 && i == j) yield CoordinateSystem.super.metricCoefficient(position, behavior, i, j);
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
		int dim = position.length;
		double[][] g = new double[dim][dim];
		for (int i = 0; i < dim; i++) {
			g[i][i] = metricCoefficient(position, behavior, i, i);
		}
		return g;
	}

	@Override
	public double jacobianDeterminant(double[] position) {
		int dim = position.length;
		return switch (dim) {
			case 3 -> Math.pow(position[0],2)*Math.abs(Math.sin(position[1]));
			case 4 -> Math.pow(position[0],3)*Math.pow(Math.sin(position[1]),2)*Math.abs(Math.sin(position[2]));
			default -> CoordinateSystem.super.jacobianDeterminant(position);
		};
	}

	@Override
	public double christoffelSymbol1stKind(double[] position, int i, int j, int k) {
		int dim = position.length;
		if (i < 0 || j < 0 || k < 0 || i >= dim || j >= dim || k >= dim) {
			throw new IllegalArgumentException(
					"i, j and k may not be <0 or exceed %d, but they were i=%d, j=%d and k=%d"
							.formatted(dim, i, j, k));
		}

		if (i == 0 && j == 1 && k == 1) return position[0];
		if (i == 0 && j == 2 && k == 2) return position[0]*Math.pow(Math.sin(position[1]),2);
		if (i == 1 && j == 0 && k == 1) return position[0];
		if (i == 1 && j == 1 && k == 0) return -1*position[0];
		if (i == 1 && j == 2 && k == 2) return 0.5*Math.pow(position[0],2)*Math.sin(2*position[1]);
		if (i == 2 && j == 0 && k == 2) return position[0]*Math.pow(Math.sin(position[1]),2);
		if (i == 2 && j == 1 && k == 2) return 0.5*Math.pow(position[0],2)*Math.sin(2*position[1]);
		if (i == 2 && j == 2 && k == 0) return -1*position[0]*Math.pow(Math.sin(position[1]),2);
		if (i == 2 && j == 2 && k == 1) return -1.0/2*Math.pow(position[0],2)*Math.sin(2*position[1]);
		
		if (dim > 4) {
			// fallback to numerics for higher dim
			return CoordinateSystem.super.christoffelSymbol1stKind(position,i, j, k);
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
		
		if (m == 0 && i == 1 && j == 1) return -1*position[0];
		if (m == 0 && i == 2 && j == 2) return -1*position[0]*Math.pow(Math.sin(position[1]),2);
		if (m == 0 && i == 3 && j == 3) return -1*position[0]*Math.pow(Math.sin(position[1]),2)*Math.pow(Math.sin(position[2]),2);
		if (m == 1 && i == 0 && j == 1) return Math.pow(position[0],-1);
		if (m == 1 && i == 1 && j == 0) return Math.pow(position[0],-1);
		if (m == 1 && i == 2 && j == 2) return -1.0/2*Math.sin(2*position[1]);
		if (m == 1 && i == 3 && j == 3) return -1*Math.pow(Math.sin(position[2]),2)*Math.cos(position[1])*Math.sin(position[1]);
		if (m == 2 && i == 0 && j == 2) return Math.pow(position[0],-1);
		if (m == 2 && i == 1 && j == 2) return Math.pow(Math.tan(position[1]),-1);
		if (m == 2 && i == 2 && j == 0) return Math.pow(position[0],-1);
		if (m == 2 && i == 2 && j == 1) return Math.pow(Math.tan(position[1]),-1);
		if (m == 2 && i == 3 && j == 3) return -1.0/2*Math.sin(2*position[2]);
		if (m == 3 && i == 0 && j == 3) return Math.pow(position[0],-1);
		if (m == 3 && i == 1 && j == 3) return Math.pow(Math.tan(position[1]),-1);
		if (m == 3 && i == 2 && j == 3) return Math.pow(Math.tan(position[2]),-1);
		if (m == 3 && i == 3 && j == 0) return Math.pow(position[0],-1);
		if (m == 3 && i == 3 && j == 1) return Math.pow(Math.tan(position[1]),-1);
		if (m == 3 && i == 3 && j == 2) return Math.pow(Math.tan(position[2]),-1);
		
		if (dim > 4) {
			// fallback to numerics for higher dim
			return CoordinateSystem.super.christoffelSymbol2ndKind(position, m, dim, j);
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
