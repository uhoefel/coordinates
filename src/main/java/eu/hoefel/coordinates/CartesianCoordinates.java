package eu.hoefel.coordinates;

import java.util.Collections;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Objects;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import eu.hoefel.coordinates.axes.Axes;
import eu.hoefel.coordinates.axes.Axis;
import eu.hoefel.coordinates.tensors.TensorIndexType;
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
 * <img src="doc-files/cartesian.svg" alt="Cartesian coordinate system"><br>
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
 * @param dimension the dimension, needs to be &gt;0. By default 1.
 *
 * @author Udo Hoefel
 */
@CoordinateSystemSymbols({"cartesian", "cart"})
public final record CartesianCoordinates(NavigableSet<Axis> axes, int dimension) implements CoordinateSystem {

    /** The default axes in 1D. */
    public static final NavigableSet<Axis> DEFAULT_AXES_1D = Axes.of(
            new Axis(Axes.DEFAULT_DIMENSION, SiBaseUnit.METER, ""),
            // special implementation for x as it is often used
            new Axis(0, SiBaseUnit.METER, "x"));

    /** The default axes in 2D. */
    public static final NavigableSet<Axis> DEFAULT_AXES_2D = Axes.of(
            new Axis(Axes.DEFAULT_DIMENSION, SiBaseUnit.METER, ""),
            // special implementation for x,y as they are often used
            new Axis(0, SiBaseUnit.METER, "x"),
            new Axis(1, SiBaseUnit.METER, "y"));

    /** The default axes in 3D. */
    public static final NavigableSet<Axis> DEFAULT_AXES_3D = Axes.of(
            new Axis(Axes.DEFAULT_DIMENSION, SiBaseUnit.METER, ""),
            // special implementation for x,y,z as they are used most frequently
            new Axis(0, SiBaseUnit.METER, "x"),
            new Axis(1, SiBaseUnit.METER, "y"),
            new Axis(2, SiBaseUnit.METER, "z"));

    /** The default axes in <i>n</i>D. */
    public static final NavigableSet<Axis> DEFAULT_AXES_ND = Axes.of(
            new Axis(Axes.DEFAULT_DIMENSION, SiBaseUnit.METER, ""));

    /**
     * Constructs a new Cartesian coordinate system.
     * 
     * @param axes      the axes defining the coordinate system
     * @param dimension the dimension, needs to be &gt;0. By default 1.
     * @throws NullPointerException     if {@code axes} are null
     * @throws IllegalArgumentException if {@code dimension} is &lt;1
     */
    public CartesianCoordinates {
        Objects.requireNonNull(axes, "Axes may not be null. "
                + "Use the DEFAULT_AXES_<your dimension>D or the empty constructor to get a reasonable default.");

        if (dimension < 1) {
            throw new IllegalArgumentException("Dimension must be >0!");
        }
    }

    /**
     * Constructs a new Cartesian coordinate system.
     * 
     * @param args the arguments, must be either {@link Axes} or {@link Axis}, which
     *             will take precedence over the default axes (which are determined
     *             dynamically depending on the chosen dimension) if given. The
     *             arguments also need to contain an integer to specify the
     *             dimension. If no {@link Axes} or {@link Axis} arguments are
     *             given, the default axes appropriate for the given dimension will
     *             be used (with special cases for 1D, 2D and 3D).
     */
    public CartesianCoordinates(Object... args) {
        // not particularly nice that we have to walk twice through the args to get the
        // dimension, but I cannot think of something nicer right now
        this(Axes.fromArgs(selectDefaultAxesFromDimension(CoordinateSystems.intFromArgs(0, args).orElse(1)), args),
                CoordinateSystems.intFromArgs(0, args).orElse(1));
    }

    /**
     * Selects the appropriate default axes for the given dimension.
     * 
     * @param dimension the dimension of the Cartesian coordinates, determining the
     *                  corresponding default axes
     * @return the default axes, never null
     */
    private static final NavigableSet<Axis> selectDefaultAxesFromDimension(int dimension) {
        return switch (dimension) {
            case 1  -> DEFAULT_AXES_1D;
            case 2  -> DEFAULT_AXES_2D;
            case 3  -> DEFAULT_AXES_3D;
            default -> DEFAULT_AXES_ND;
        };
    }

    /**
     * Validates the position, i.e. it throws an exception if a dimension of the 
     * position is out of range.
     * 
     * @param position the position to validate
     * @throw IllegalArgumentException if the assumptions about the dimensionality 
     *          or the valid range of any dimension of the input are violated.
     */
    private void validatePosition(double[] position) {
        Objects.requireNonNull(position);

        if (position.length > dimension) {
            throw new IllegalArgumentException(
                    "The given dimensionality exceeds the specified dimensionality (%d vs %d)"
                            .formatted(position.length, dimension));
        }
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
        validatePosition(position);
        return position.clone();
    }

    @Override
    public double[] fromBasePoint(double[] position) {
        return position.clone();
    }

    @Override
    public double[] toBaseVector(double[] position, double[] vector) {
        validatePosition(position);
        return vector.clone();
    }

    @Override
    public double[] fromBaseVector(double[] position, double[] vector) {
        return vector.clone();
    }

    @Override
    public double metricCoefficient(double[] position, TensorTransformation behavior, int i, int j) {
        validatePosition(position);

        if (i < 0 || j < 0) {
            throw new IllegalArgumentException(("Metric coefficient not available for i=%d, j=%d "
                    + "(too low dimension, only 3 dimensions [0,1,2] are supported for Cartesian coordinates)")
                    .formatted(i, j));
        } else if (i > dimension || j > dimension) {
            throw new IllegalArgumentException(("Metric coefficient not available for i=%d, j=%d "
                    + "(too high dimension, only %d dimensions are supported for Cartesian coordinates)")
                    .formatted(i, j, dimension));
        }

        return i == j ? 1 : 0;
    }

    @Override
    public double[][] metricTensor(double[] position, TensorTransformation behavior) {
        validatePosition(position);
        return Maths.eye(position.length);
    }

    @Override
    public double jacobianDeterminant(double[] position) {
        validatePosition(position);
        return 1;
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
    public double[] curl(double[] position, TensorTransformation componentBehavior, UnaryOperator<double[]> vectorfield) {
        validatePosition(position);
        return CoordinateSystem.super.curl(position, componentBehavior, vectorfield);
    }
}
