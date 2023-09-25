package eu.hoefel.coordinates;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Objects;
import java.util.Set;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.function.UnaryOperator;
import java.util.stream.IntStream;

import eu.hoefel.coordinates.axes.Axes;
import eu.hoefel.coordinates.axes.Axis;
import eu.hoefel.coordinates.tensors.TensorIndexType;
import eu.hoefel.coordinates.tensors.TensorTransformation;
import eu.hoefel.unit.Unit;
import eu.hoefel.utils.Maths;
import eu.hoefel.utils.Reflections;
import eu.hoefel.utils.Types;

/**
 * A coordinate system is a system allowing coordinates (represented by one or
 * more ordered numbers) to mark a position on a manifold (e.g. Euclidean space,
 * here represented by the {@link CartesianCoordinates}). This allows to treat
 * geometrical problems analytically.
 * <p>
 * Note further that a coordinate system may or may not require state. For
 * example, {@link BipolarCoordinates} require state (one needs to define where
 * the focii are positioned), while e.g. {@link CylindricalCoordinates} can have
 * state (e.g. by having non-default {@link Unit units} on the axes), but do not
 * require it.
 * <p>
 * Note that for some abstract systems complex numbers may be necessary. This is
 * currently not supported. It may be worth thinking about a change there once
 * project <a href="https://openjdk.java.net/projects/valhalla/">Valhalla</a> is
 * implemented.
 * 
 * @author Udo Hoefel
 * @see <a href=
 *      "https://csm.mech.utah.edu/content/wp-content/uploads/2011/03/GoBagCurvilinear.pdf">Curvilinear
 *      Analysis in a Euclidean Space by Rebecca Moss Brannon</a>
 * @see <a href="http://faculty.ce.berkeley.edu/sanjay/ce231mse211/curvi.pdf">A
 *      Quick Overview of Curvilinear Coordinates</a>
 * @apiNote Implementations have to either implement
 *          {@link #toBasePoint(double[])} and {@link #fromBasePoint(double[])}.
 *          A limited subset of features can be provided by providing the
 *          {@link #metricTensor(double[], TensorTransformation)} alone. It is
 *          strongly recommended to use a {@link java.lang.Record record} for
 *          the implementation, otherwise it may not work directly with
 *          {@code eu.hoefel.quantities}.
 */
public interface CoordinateSystem {

    /**
     * Gets the (defined) axes of this coordinate system. Note that it may contain
     * an axis with the {@link Axes#DEFAULT_DIMENSION}.
     * 
     * @return the axes
     */
    public NavigableSet<Axis> axes();

    /**
     * Gets the axis of the specified dimension.
     * 
     * @param dimension the dimension
     * @return the axis corresponding to the specified dimension
     */
    default Axis axis(int dimension) {
        return Axis.fromSet(axes(), dimension);
    }

    /**
     * Gets the dimensionality that the coordinate system handles. It can handle
     * <em>only</em> that dimensionality.
     * <p>
     * For example:<br>
     * {@link CartesianCoordinates} can handle <em>whatever</em> dimensionality
     * you specify, while {@link CylindricalCoordinates} can <em>only</em>
     * handle 3D coordinates.
     * 
     * @return the dimensionality that can be handled by the coordinate system
     */
    public int dimension();

    /**
     * Gets the class that represents the base coordinate system implementation to
     * the current coordinate system.
     * 
     * @return the base coordinates
     */
    default Class<? extends CoordinateSystem> baseCoordinates() {
        return CartesianCoordinates.class;
    }

    /**
     * Transforms the given position into the base coordinate system of the current
     * coordinate system.
     * 
     * @param position the position in the current coordinate system
     * @return the position in the base coordinate system
     */
    public double[] toBasePoint(double[] position);

    /**
     * Gets the units of the axes of the base coordinate system corresponding to the
     * current coordinate system.
     * <p>
     * An example code for cylindrical coordinates could look like this:
     * <p>
     * {@code var map = Map.of(0, axis(0).unit(), 1, axis(0).unit(), 2, axis(2).unit());}<br>
     * {@code return Collections.unmodifiableNavigableMap(new TreeMap<>(map));}
     * <p>
     * Note that default axes may be put in via {@link Axes#DEFAULT_DIMENSION}.
     * <p>
     * Due to Javas lack of operator overloading it has to be specified separately
     * from {@link #toBasePoint(double[])}.
     * 
     * @return the map containing the axis dimension and the corresponding unit
     */
    public NavigableMap<Integer, Unit> toBaseUnits();

    /**
     * Transforms the given position from the base coordinate system into the
     * current coordinate system.
     * 
     * @param position the position in the base coordinate system
     * @return the position in the current coordinate system
     */
    public double[] fromBasePoint(double[] position);

    /**
     * Transforms the given vector at the specified position into the base
     * coordinate system.
     * 
     * @param position the position (in the current coordinate system) at which to
     *                 evaluate the derivatives necessary for the transformation.
     *                 May not be null.
     * @param vector   the vector (with respect to the natural basis) to transform.
     *                 May not be null.
     * @return the vector in the base coordinate system (using natural basis
     *         vectors, i.e. they are not necessarily normalized)
     */
    default double[] toBaseVector(double[] position, double[] vector) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(vector);

        return vectorTransformation(position, vector, this::toBasePoint);
    }

    /**
     * Transforms the given vector at the specified position from the base
     * coordinate system to this coordinate system.
     * 
     * @param position the position (in the base coordinate system) at which to
     *                 evaluate the derivatives necessary for the transformation.
     *                 May not be null.
     * @param vector   the vector (with respect to the natural basis) to transform.
     *                 May not be null
     * @return the vector in this coordinate system (using natural basis vectors,
     *         i.e. they are not necessarily normalized)
     */
    default double[] fromBaseVector(double[] position, double[] vector) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(vector);

        return vectorTransformation(position, vector, this::fromBasePoint);
    }

    /**
     * Transforms the vector at the given position, either to or from the base
     * coordinate system.
     * 
     * @param position       the position at which to evaluate the derivatives
     *                       necessary for the transformation
     * @param vector         the vector (with respect to the natural basis) to
     *                       transform
     * @param transformation the transformation to use for the calculation of the
     *                       partial derivatives
     * @return the transformed vector (using natural basis vectors, i.e. they are
     *         not necessarily normalized)
     */
    private double[] vectorTransformation(double[] position, double[] vector, UnaryOperator<double[]> transformation) {
        double[][] m = new double[vector.length][vector.length];
        for (int i = 0; i < m.length; i++) {
            double[] d1 = Maths.partialDerivatives(transformation, 1, position, i);
            for (int j = 0; j < d1.length; j++) {
                m[j][i] = d1[j];
            }
        }

        return Maths.matrixVectorMul(m, vector);
    }

    /**
     * Gets the symbols representing this coordinate system, e.g. "cartesian" and
     * "cart" for {@link CartesianCoordinates}.
     * 
     * @return the list of symbols, starting with the preferred symbol
     * 
     * @implNote The default implementation takes the symbols as obtained from the
     *           {@link CoordinateSystemSymbols} annotation that is recommended to be
     *           present for the coordinate system implementation. The default
     *           implementation will use an empty list if this annotation is not
     *           present, which means that string based methods like e.g.
     *           {@link CoordinateSystems#transform(double[], String, String, Set)}
     *           will not recognize the coordinate system if the given arguments do
     *           not allow a successful instantiation.
     */
    default List<String> symbols() {
        return CoordinateSystems.symbolsFromClass(this.getClass());
    }

    /**
     * Checks whether the current coordinate system is a base coordinate system
     * (i.e., a coordinate system that will be used in coordinate transformations to
     * avoid the combinatorial explosion that otherwise would arise if every
     * coordinate system would have to implement conversions to all other coordinate
     * systems).
     * 
     * @return true if the coordinate system is a base coordinate system
     */
    default boolean isBasic() {
        return false;
    }

    /**
     * Gets a new coordinate system from the specified representative name and the
     * given arguments. Note that an empty name will always return
     * {@link CoordinateSystems#IDENTITY_COORDINATE_SYSTEM}.
     * <p>
     * Example usage:<br>
     * <code>
     * CoordinateSystem c = CoordinateSystem.from("cart"); // standard "stateless" cartesian coordinates<br>
     * c = CoordinateSystem.from("bipol", 2); // stateful bipolar coordinates (a=2)<br>
     * </code>
     * 
     * @param name the name representing the coordinate system, e.g. "cart" for the
     *             {@link CartesianCoordinates}. May not be null.
     * @param args the arguments to pass on to the constructor that fits most
     *             closely the signature of the arguments. Can handle varargs. May
     *             not be null.
     * @return a new instance of the specified coordinate system
     */
    public static CoordinateSystem from(String name, Object... args) {
        Objects.requireNonNull(name);
        Objects.requireNonNull(args);

        return from(name, CoordinateSystems.DEFAULT_COORDINATE_SYSTEMS, args);
    }

    /**
     * Gets a new coordinate system from the specified representative name and the
     * given arguments. Note that an empty name will always return
     * {@link CoordinateSystems#IDENTITY_COORDINATE_SYSTEM}.
     * <p>
     * Example usage:<br>
     * <code>
     * Set&lt;Class&lt;? extends CoordinateSystem&gt;&gt; s = Set.of(CartesianCoordinates.class, ToroidalCoordinates.class);<br>
     * CoordinateSystem c = CoordinateSystem.from("cart", s);<br>
     * c = CoordinateSystem.from("tor", s, 2); // stateful toroidal coordinates (a=2)<br>
     * </code>
     * 
     * @param name             the name representing the coordinate system, e.g.
     *                         "cart" for the {@link CartesianCoordinates}. May not
     *                         be null.
     * @param extraCoordinates the coordinate systems of which one contains the
     *                         specified name as an allowed representation. May not
     *                         be null.
     * @param args             the arguments to pass on to the constructor that fits
     *                         most closely the signature of the arguments. Can
     *                         handle varargs. May not be null.
     * @return a new instance of the specified coordinate system
     */
    public static CoordinateSystem from(String name, Set<? extends Class<? extends CoordinateSystem>> extraCoordinates, Object... args) {
        Objects.requireNonNull(name);
        Objects.requireNonNull(extraCoordinates);
        Objects.requireNonNull(args);

        if (name.isEmpty()) return CoordinateSystems.IDENTITY_COORDINATE_SYSTEM;

        List<String> knownSymbols = new ArrayList<>();
        for (var clazz : extraCoordinates) {
            List<String> symbols = CoordinateSystems.symbolsFromClass(clazz);
            if (symbols.contains(name)) {
                return Reflections.newInstance(clazz, args);
            }
            knownSymbols.addAll(symbols);
        }

        throw new IllegalArgumentException(
                "Unknown coordinate system (%s). The known symbols are %s".formatted(name, knownSymbols));
    }

    /**
     * The Jacobian determinant (also known as functional determinant) contains
     * important information about the behavior of the metric near the specified
     * coordinates. For example, if the determinant is positive the orientation of
     * the metric is preserved at the specified point, while it reverses orientation
     * if the determinant is negative.
     * 
     * @param position the position, not null
     * @return the Jacobian determinant at the specified coordinates
     * @see <a href=
     *      "https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant#Jacobian_determinant">Wikipedia</a>
     */
    default double jacobianDeterminant(double[] position) {
        Objects.requireNonNull(position);
        return Math.sqrt(Maths.determinant(metricTensor(position, TensorIndexType.COVARIANT)));
    }

    /**
     * Checks whether the coordinate system is orthogonal, i.e. all non-diagonal
     * metric coefficients are zero.
     * 
     * @return true if the coordinate system is orthogonal everywhere
     */
    public boolean isOrthogonal();

    /**
     * Gets the natural basis vectors. Note that the components are expressed in the
     * base coordinate system.
     * 
     * @param position the position at which the natural basis vectors are requested
     * @param type     the desired component index type of the basis vectors
     * @return the natural basis vectors in the format [i][j] with i the index for
     *         the basis vector, i.e. <i>g</i><sub><i>i</i></sub> (or
     *         <i>g</i><sup><i>i</i></sup>) and j the index in the <i>i</i>th base
     *         vector
     */
    private double[][] naturalBasisVectors(double[] position, TensorIndexType type) {
        // natural base vectors g_i
        double[][] gi = new double[position.length][position.length];
        for (int i = 0; i < gi.length; i++) {
            gi[i] = Maths.partialDerivatives(this::toBasePoint, 1, position, i);
        }

        if (type == TensorIndexType.COVARIANT) {
            return gi;
        }

        var indices = TensorTransformation.with(TensorIndexType.COVARIANT, TensorIndexType.CONTRAVARIANT);
        return indices.transform(this, pos -> gi, TensorIndexType.CONTRAVARIANT).apply(position);
    }

    /**
     * Calculates the magnitude of a tensor field.
     * 
     * @param <T>            specifies the rank of the tensor, i.e. for a 0th order
     *                       tensor (scalar) field it should be Double, for a first
     *                       order tensor (vector) field it should be double[] and
     *                       so on
     * @param position       the position at which to evaluate the tensor field, not
     *                       null
     * @param transformation the transformation information (i.e. the index
     *                       positions) of the given tensor field, not null
     * @param tensorfield    the tensor field of arbitrary order, not null
     * @return the magnitude of the tensor field
     */
    @SuppressWarnings("unchecked")
    default <T> double magnitude(double[] position, TensorTransformation transformation, Function<double[], T> tensorfield) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(transformation);
        Objects.requireNonNull(tensorfield);

        T tensor = tensorfield.apply(position);
        if (tensor instanceof Double) {
            // tensor of 0th order
            return Math.abs((double) tensor);
        }

        // tensor of nth order
        var covariantField = transformation.transform(this, tensorfield, TensorIndexType.COVARIANT);
        var contravariantField = transformation.transform(this, tensorfield, TensorIndexType.CONTRAVARIANT);

        T coTensor = covariantField.apply(position);
        T contraTensor = contravariantField.apply(position);

        int dim = Array.getLength(tensor);
        int rank = Types.dimension(tensor.getClass());
        double[] summands = new double[(int) Math.pow(dim, rank)];

        Object coTensorFlattening = coTensor;
        Object contraTensorFlattening = contraTensor;

        for (int i = 0; i < rank - 1; i++) {
            // these casts are safe
            coTensorFlattening = Maths.flatten((T[]) coTensorFlattening);
            contraTensorFlattening = Maths.flatten((T[]) contraTensorFlattening);
        }

        // we necessarily have a double[] now
        double[] flatCoTensor = (double[]) coTensorFlattening;
        double[] flatContraTensor = (double[]) contraTensorFlattening;

        for (int i = 0; i < summands.length; i++) {
            summands[i] = flatCoTensor[i] * flatContraTensor[i];
        }

        return Math.sqrt(Maths.compensatedSum(summands));
    }

    /**
     * Gets the metric tensor at the given coordinates for the specified
     * transformation behavior (i.e., either {@link TensorIndexType#COVARIANT
     * covariant}, {@link TensorIndexType#CONTRAVARIANT contravariant} or mixed).
     * 
     * @param position the position at which to get the metric tensor, not null
     * @param behavior the tensor transformation behavior of the metric
     *                 coefficients, not null
     * @return the metric tensor
     */
    default double[][] metricTensor(double[] position, TensorTransformation behavior) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(behavior);

        int dim = position.length;
        double[][] partialDerivatives = new double[dim][dim];

        for (int i = 0; i < dim; i++) {
            partialDerivatives[i] = Maths.partialDerivatives(this::toBasePoint, 1, position, i);
        }

        double[][] g = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                for (int k = 0; k < dim; k++) {
                    g[i][j] += partialDerivatives[i][k] * partialDerivatives[j][k];
                }
                g[j][i] = g[i][j];
            }
        }

        // In principle we could use something along the lines of
        // TensorIndexType.COVARIANT.transform(this, pos -> g, behavior).apply(position)
        // but this transformation relies on the metricTensor method we are in here, so 
        // we would get a stackoverflow error. Consequently, we just invert the metric.

        if (behavior == TensorIndexType.CONTRAVARIANT) {
            g = Maths.inverse(g);
        }

        return g;
    }
    
    /**
     * Gets a single metric coefficient from the metric tensor for the specified
     * tensor transformation behavior (i.e., either {@link TensorIndexType#COVARIANT
     * covariant}, {@link TensorIndexType#CONTRAVARIANT contravariant} or mixed).
     * 
     * @param position the position at which to get the metric tensor coefficient,
     *                 not null
     * @param behavior the tensor transformation behavior of the metric coefficient,
     *                 not null
     * @param i        the row index in the metric tensor
     * @param j        the column index in the metric tensor
     * @return the metric tensor coefficient
     * 
     * @implNote The default implementation gets the full metric tensor (which might
     *           be an expensive operation) and returns just the single coefficient
     *           the user asked for. If multiple coefficients are required it is
     *           advisable to use {@link #metricTensor(double[], TensorTransformation)}
     *           and access the coefficients from the locally stored full tensor.
     */
    default double metricCoefficient(double[] position, TensorTransformation behavior, int i, int j) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(behavior);

        return metricTensor(position, behavior)[i][j];
    }

    /**
     * Gets the line element for a specific dimension.
     * 
     * @param position the position (with contravariant components) at which to get
     *                 the line element, not null
     * @param i        the dimension for which to get the line element (dimension
     *                 indices start at 0)
     * @param dui      the (small) change with respect to the <i>i</i>th coordinate
     * @return the line element for the specified dimension, i.e.
     *         &radic;<span style="text-decoration:
     *         overline"><i>g</i><sub><i>ii</i></sub></span>
     *         d<i>u</i><sub><i>i</i></sub>
     * 
     * @see <a href=
     *      "https://personalpages.manchester.ac.uk/staff/Andrew.Hazel/MATH45061/MATH45061_Ch1.pdf">Describing
     *      the Physical World: Vectors &amp; Tensors, page 18</a>
     * 
     * @throws IllegalArgumentException if i is greater than the supported
     *                                  dimensionality
     */
    default double ds(double[] position, int i, double dui) {
        Objects.requireNonNull(position);

        if (i >= dimension()) {
            throw new IllegalArgumentException(
                    "Specified dimension greater than supported dimensionality (given: %d, supported: %d)"
                    .formatted(i + 1, dimension()));
        }

        double gii = metricCoefficient(position, TensorIndexType.COVARIANT, i, i);
        return Math.sqrt(gii) * dui;
    }

    /**
     * Gets the surface element for the specified dimensions.
     * 
     * @param position the position (with contravariant components) at which to get
     *                 the surface element, not null
     * @param i        the first of the two dimensions for which to get the surface
     *                 element (dimension indices start at 0)
     * @param j        the second of the two dimensions for which to get the surface
     *                 element (dimension indices start at 0)
     * @param dui      the (small) change with respect to the <i>i</i>th coordinate
     * @param duj      the (small) change with respect to the <i>j</i>th coordinate
     * @return the surface element for the specified dimensions, i.e.
     *         &radic;<span style="text-decoration:
     *         overline"><i>g</i><sub><i>ii</i></sub><i>g</i><sub><i>jj</i></sub> -
     *         <i>g</i><sub><i>ij</i></sub><i>g</i><sub><i>ji</i></sub></span>
     *         d<i>u</i><sub><i>i</i></sub>d<i>u</i><sub><i>j</i></sub>
     * 
     * @see <a href=
     *      "https://personalpages.manchester.ac.uk/staff/Andrew.Hazel/MATH45061/MATH45061_Ch1.pdf">Describing
     *      the Physical World: Vectors &amp; Tensors, page 18</a>
     * 
     * @throws IllegalArgumentException if i or j are greater than the supported
     *                                  dimensionality, or i==j as this does not
     *                                  create a surface
     */
    default double dA(double[] position, int i, int j, double dui, double duj) {
        Objects.requireNonNull(position);

        if (i == j) {
            throw new IllegalArgumentException(
                    "To get a surface element the given dimensions i and j must be different, but you specified i=j=%d"
                    .formatted(i));
        } else if (i >= dimension()) {
            throw new IllegalArgumentException(
                    "Specified dimension for i greater than supported dimensionality (given: %d, supported: %d)"
                    .formatted(i + 1, dimension()));
        } else if (j >= dimension()) {
            throw new IllegalArgumentException(
                    "Specified dimension for j greater than supported dimensionality (given: %d, supported: %d)"
                    .formatted(j + 1, dimension()));
        }

        double[][] g = metricTensor(position, TensorIndexType.COVARIANT);
        return Math.sqrt(g[i][i] * g[j][j] - g[i][j] * g[j][i]) * dui * duj;
    }

    /**
     * Gets the volume element at the given position.
     * <p>
     * It measures the volume of a small (in principle, infinitesimally small, in
     * practice it depends on the given du) parallelepiped. In general it is defined
     * to be
     * <p>
     * d<i>V</i> = |det(<i>J</i>)|
     * d<i>u</i><sub>1</sub>d<i>u</i><sub>2</sub>...d<i>u</i><sub><i>n</i></sub>,
     * <p>
     * with <i>J</i> being the {@link #jacobianDeterminant(double[]) Jacobian
     * determinant}.
     * <p>
     * From a mathematical perspective, a volume form on a <i>n</i>-dimensional
     * oriented and differentiable manifold is a differential form of degree
     * <i>n</i> that vanishes nowhere.
     * 
     * @param position the position (with contravariant components) at which to get
     *                 the volume element. May not be null.
     * @param du       the (small) changes with respect to the 1,2,...,<i>n</i>th
     *                 coordinate, not null
     * @return the volume element
     * 
     * @see <a href=
     *      "https://personalpages.manchester.ac.uk/staff/Andrew.Hazel/MATH45061/MATH45061_Ch1.pdf">Describing
     *      the Physical World: Vectors &amp; Tensors, page 19</a>
     * 
     * @throws IllegalArgumentException if du is null or the length of the given du
     *                                  is not the same as the length of the
     *                                  position
     */
    default double dV(double[] position, double... du) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(du);

        if (du.length != position.length) {
            throw new IllegalArgumentException(
                    "A volume form needs to be of the same dimensionality as the given coordinates. "
                            + "You specified a %d dimensional position, but du was of dimensionality %d"
                            .formatted(position.length, du.length));
        }

        double ret = Math.abs(jacobianDeterminant(position));
        for (int i = 0; i < du.length; i++) {
            ret *= du[i];
        }
        return ret;
    }

    /**
     * Calculates the dot (also known as vector inner or scalar) product.
     * 
     * @param position the position at which to evaluate the dot product, not null
     * @param behavior the tensor index type for the components of v1 and v2, not
     *                 null
     * @param v1       the first vector. Needs to be of the same length as v2. May
     *                 not be null.
     * @param v2       the second vector. Needs to be of the same length as v1. May
     *                 not be null.
     * @return the dot product
     */
    default double dot(double[] position, TensorIndexType behavior, double[] v1, double[] v2) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(behavior);
        Objects.requireNonNull(v1);
        Objects.requireNonNull(v2);

        int dims = position.length;
        double[][] g = metricTensor(position, behavior.flip());
        double[] summands = new double[dims * dims];
        for (int i = 0; i < dims; i++) {
            for (int j = 0; j < dims; j++) {
                summands[i * dims + j] = g[i][j] * v1[i] * v2[j];
            }
        }

        return Maths.compensatedSum(summands);
    }

    /**
     * Calculates the cross product of v1 and v2 at the given position. Note that
     * for dimensions &gt;3 this is generalized to be the Hodge dual of the exterior
     * product and <i>n</i>-1 vectors, also called external product.
     * 
     * @param position the position at which to evaluate the cross/external product,
     *                 not null
     * @param behavior the transformation property of v1 and v2 (and potentially
     *                 vn), either {@link TensorIndexType#COVARIANT covariant} or
     *                 {@link TensorIndexType#CONTRAVARIANT contravariant}, not null
     * @param v1       the first vector, not null
     * @param v2       the second vector, not null
     * @param vn       additional vectors if calculated for dimensions &gt;3, not
     *                 null
     * @return the cross/external product with the basis change behavior flipped as
     *         compared to the given vectors
     */
    default double[] cross(double[] position, TensorIndexType behavior, double[] v1, double[] v2, double[]... vn) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(behavior);
        Objects.requireNonNull(v1);
        Objects.requireNonNull(v2);
        Objects.requireNonNull(vn);

        int dim = position.length;

        if (vn.length + 2 > dim - 1) {
            throw new IllegalArgumentException("Wrong number of vectors (%d) given for the requested dimension (%d)"
                    .formatted(vn.length + 2, dim));
        }

        double[] returnVector = new double[dim];

        // collect all dim-1 vectors
        double[][] vs = new double[vn.length + 2][dim];
        vs[0] = v1;
        vs[1] = v2;
        System.arraycopy(vn, 0, vs, 2, vn.length);

        double sqrtg = jacobianDeterminant(position);

        int[][] permutations = Maths.permutations(IntStream.rangeClosed(1, dim).toArray());
        for (int[] sumIndices : permutations) {
            int leviCivita = Maths.leviCivita(sumIndices);
            if (leviCivita == 0) continue;

            double factor = 1;
            for (int j = 0; j < vs.length; j++) {
                factor *= vs[j][sumIndices[j + dim - vs.length] - 1];
            }

            returnVector[sumIndices[0] - 1] += leviCivita * factor / sqrtg;
        }

        return behavior.transform(this, pos -> returnVector, behavior.flip()).apply(position);
    }

    /**
     * Calculates the divergence of the given tensor field.
     * <p>
     * Note that due to current limitations with primitives/generics one cannot use
     * a {@code Function<double[],double[]}, but one will have to fall back using a
     * {@code Function<double[],Double[]}.
     * 
     * @param <T>               the type of the return of the tensor field. If the
     *                          tensor field is of first order this must be
     *                          {@code double[]}, if it is of second order it needs
     *                          to be {@code double[][]}.
     * @param position          the position at which the divergence should be
     *                          calculated, not null
     * @param componentBehavior the tensor field component behavior, not null
     * @param field             the tensor field of first order (i.e. a vector
     *                          field) or second order (i.e. a matrix field), not
     *                          null
     * @return the divergence of the tensor field at the specified position
     * @throws ClassCastException if any of the conditions for T are violated
     */
    default <T> T div(double[] position, TensorTransformation componentBehavior, Function<double[], T[]> field) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(componentBehavior);
        Objects.requireNonNull(field);

        Object t = field.apply(position);
        if (t instanceof double[]) {
            @SuppressWarnings("unchecked")
            T div = (T) divergenceOfVectorField(position, componentBehavior, pos -> (double[]) (Object) field.apply(pos));
            return div;
        } else if (t instanceof Double[]) { // ugh, will be good when this will be replaceable with double[]. Valhalla probably...
            @SuppressWarnings("unchecked")
            T div = (T) divergenceOfVectorField(position, componentBehavior, pos -> {
                var applied = field.apply(pos);
                return Types.unbox(applied);
            });
            return div;
        } else if (t instanceof double[][] tensor) {
            @SuppressWarnings("unchecked")
            T div = (T) divergenceOfMatrixField(position, tensor, componentBehavior, pos -> (double[][]) field.apply(pos));
            return div;
        }

        throw new IllegalArgumentException("f needs to return either a double[] (i.e. a vector) or a double[][] (i.e. a tensor)");
    }

    /**
     * Calculates the divergence of the given vector field.
     * 
     * @param position          the position at which the divergence should be
     *                          calculated
     * @param componentBehavior the vector field component behavior
     * @param vectorfield       the vector field
     * @return the vector field divergence at the given position
     */
    private Double divergenceOfVectorField(double[] position, TensorTransformation componentBehavior, UnaryOperator<double[]> vectorfield) {
        double sqrtg = jacobianDeterminant(position);

        UnaryOperator<double[]> detTimesVector = pos -> {
            double jd = jacobianDeterminant(pos);
            double[] vec = componentBehavior.transform(this, vectorfield, TensorIndexType.CONTRAVARIANT).apply(pos);
            for (int i = 0; i < vec.length; i++) {
                vec[i] *= jd;
            }
            return vec;
        };

        double[] div = new double[position.length];
        for (int i = 0; i < position.length; i++) {
            div[i] = Maths.partialDerivatives(detTimesVector, 1, position, i)[i];
        }

        return Maths.compensatedSum(div) / sqrtg;
    }

    /**
     * Calculates the divergence of the given tensor field. The tensor is of second
     * rank, i.e. it is a matrix field.
     * 
     * @param position          the position at which the divergence should be
     *                          calculated
     * @param tensor            the tensor of second rank at the given position
     * @param componentBehavior the tensor field transformation behavior
     * @param field             the tensor (of rank 2) field
     * @return the tensor field divergence at the given position
     */
    private double[] divergenceOfMatrixField(double[] position, double[][] tensor,
            TensorTransformation componentBehavior, Function<double[], double[][]> field) {

        int dimension = position.length;
        double[][] summands = new double[tensor.length][(int) (Math.pow(tensor.length, 2) + 2*Math.pow(tensor.length, 3))];

        // we need σ^ij further down
        var tensorfield = componentBehavior.transform(this, field, TensorIndexType.CONTRAVARIANT);

        // this gives the natural base vectors g_i (note that g is _not_ the metric)
        // that are used for the summations below
        double[][] gi = naturalBasisVectors(position, TensorIndexType.COVARIANT);

        int[] counters = new int[tensor.length];
        for (int vecIndex = 0; vecIndex < tensor.length; vecIndex++) {
            for (int i = 0; i < tensor.length; i++) {
                // ∂σ^ij / ∂z^j
                final int dimi = i; // final for lambda
                for (int j = 0; j < dimension; j++) {
                    final int dimj = j; // final for lambda

                    DoubleUnaryOperator func = pos -> tensorfield.apply(Maths.updatePosition(position, dimj, pos))[dimi][dimj];
                    summands[vecIndex][counters[vecIndex]++] = gi[i][vecIndex] * Maths.derivative(func, 1, position[j]);
                }

                // σ^lj Γ^i_lj
                for (int l = 0; l < dimension; l++) {
                    for (int j = 0; j < dimension; j++) {
                        summands[vecIndex][counters[vecIndex]++] = gi[i][vecIndex] * tensor[l][j] * christoffelSymbol2ndKind(position, i, l, j);
                    }
                }

                // σ^il Γ^j_lj
                for (int l = 0; l < dimension; l++) {
                    for (int j = 0; j < dimension; j++) {
                        summands[vecIndex][counters[vecIndex]++] = gi[i][vecIndex] * tensor[i][l] * christoffelSymbol2ndKind(position, j, l, j);
                    }
                }
            }
        }

        double[] div = new double[tensor.length];
        for (int i = 0; i < summands.length; i++) {
            div[i] = Maths.compensatedSum(summands[i]);
        }

        return div;
    }

    /**
     * Calculates the physical gradient of the given tensor field. This method can
     * handle scalar (tensor of zeroth order), vector (tensor of first order) and
     * matrix (tensor of second order) fields.
     * 
     * @param <T>               the type of the return of the tensor field. If the
     *                          tensor field is of zeroth order this must be
     *                          {@code Double}, if it is of first order this must be
     *                          {@code double[]} and if it is of second order it
     *                          needs to be {@code double[][]}.
     * @param position          the position at which the gradient should be
     *                          calculated, not null
     * @param componentBehavior the tensor field component behavior, not null
     * @param field             the field. May be a scalar, vector or tensor/matrix
     *                          field, not null
     * @return the physical gradient (i.e. the gradient in terms of the physical
     *         basis)
     * @throws IllegalArgumentException if the conditions for T are violated
     */
    default <T> T[] grad(double[] position, TensorTransformation componentBehavior, Function<double[], T> field) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(componentBehavior);
        Objects.requireNonNull(field);

        T t = field.apply(position);
        if (t instanceof Double) {
            // gradient of scalar field
            @SuppressWarnings("unchecked")
            T[] grad = (T[]) gradientOfScalarField(position, (Function<double[], Double>) field);
            return grad;
        } else if (t instanceof double[]) {
            // gradient of vector field
            @SuppressWarnings("unchecked")
            T[] grad = (T[]) gradientOfVectorField(position, componentBehavior, (Function<double[], double[]>) field);
            return grad;
        } else if (t instanceof double[][] tensor) {
             // gradient of matrix field
            throw new UnsupportedOperationException("Currently not implemented. "
                    + "If you need this please create an issue at https://github.com/uhoefel/coordinates");
//            @SuppressWarnings("unchecked")
//            T[] grad = (T[]) gradientOfMatrixField(position, tensor, componentBehavior, (Function<double[], double[][]>) field);
//            return grad;
        } else {
            throw new IllegalArgumentException(
                    "f needs to return either a Double (i.e. it is a scalar field), a double[] (i.e. i is a vector field) "
                    + "or a double[][] (i.e. it is a tensor field)");
        }
    }

    /**
     * Calculates the physical gradient of a scalar field.
     * 
     * @param position the position at which the gradient should be calculated
     * @param field    the scalar field
     * @return the physical gradient in terms of the physical basis
     */
    private Double[] gradientOfScalarField(double[] position, Function<double[],Double> field) {
        // calculate derivatives only once
        double[] derivatives = new double[position.length];
        for (int i = 0; i < position.length; i++) {
            // ∂f / ∂z^i
            final int dimi = i;
            DoubleUnaryOperator func = pos -> field.apply(Maths.updatePosition(position, dimi, pos));
            derivatives[i] = Maths.derivative(func, 1, position[i]);
        }

        Double[] gradient; // boxing is not so nice, will be gone with Valhalla
        if (isOrthogonal()) {
            gradient = Types.box(derivatives.clone());
            for (int i = 0; i < gradient.length; i++) {
                gradient[i] /= Math.sqrt(metricCoefficient(position, TensorIndexType.COVARIANT, i, i));
            }
        } else {
            gradient = new Double[derivatives.length];
            double[][] gContra = metricTensor(position, TensorIndexType.CONTRAVARIANT);
            for (int i = 0; i < gradient.length; i++) {
                // (∇u)_i = h_i * g^ij * (∂f / ∂z^j)
                double[] gijAj = new double[derivatives.length];
                for (int j = 0; j < derivatives.length; j++) {
                    gijAj[j] = gContra[i][j] * derivatives[j];
                }
                double hi = Math.sqrt(metricCoefficient(position, TensorIndexType.COVARIANT, i, i));
                gradient[i] = hi * Maths.compensatedSum(gijAj);
            }
        }
        return gradient;
    }

    /**
     * Gets the physical gradient of a vector field.
     * 
     * @param position          the position at which the gradient should be
     *                          calculated
     * @param componentBehavior the vector field component behavior
     * @param field             the vector field
     * @return the physical gradient in terms of the physical basis
     */
    private double[][] gradientOfVectorField(double[] position, TensorTransformation componentBehavior,
            Function<double[],double[]> field) {
        var fieldWithCovariantComponents = componentBehavior.transform(this, field, TensorIndexType.COVARIANT);
        double[] vector = fieldWithCovariantComponents.apply(position);
        double[][] gradient = new double[vector.length][vector.length];
        double[][] g = metricTensor(position, TensorIndexType.COVARIANT);

        // see https://www.iith.ac.in/~ashok/Maths_Lectures/TutorialB/tprimer.pdf equation (107)
        if (isOrthogonal()) {
            for (int i = 0; i < gradient.length; i++) {
                for (int j = 0; j < gradient.length; j++) {
                    // ∂_i(h_j A_j)
                    final int dimi = i; // for lambda
                    final int dimj = j; // for lambda
                    DoubleUnaryOperator func = pos -> {
                        double[] newPosition = Maths.updatePosition(position, dimi, pos);
                        double hj = Math.sqrt(metricCoefficient(newPosition, TensorIndexType.COVARIANT, dimj, dimj));
                        double Aj = fieldWithCovariantComponents.apply(newPosition)[dimj];
                        return hj * Aj;
                    };

                    double dihjAj = Maths.derivative(func, 1, position[j]);

                    double[] GkijhkAk = new double[g.length];
                    for (int k = 0; k < g.length; k++) {
                        GkijhkAk[k] = christoffelSymbol2ndKind(position, k, i, j) * Math.sqrt(g[k][k]) * vector[k];
                    }
                    gradient[i][j] = (dihjAj - Maths.compensatedSum(GkijhkAk)) / Math.sqrt(g[i][i] * g[j][j]);
                }
            }
        } else {
            double[][] gContra = metricTensor(position, TensorIndexType.CONTRAVARIANT);

            // via eq. (32) in https://www.iith.ac.in/~ashok/Maths_Lectures/TutorialB/tprimer.pdf
            for (int p = 0; p < gradient.length; p++) {
                for (int q = 0; q < gradient.length; q++) {
                    double gpigqjTij = 0;
                    for (int i = 0; i < gradient.length; i++) {
                        for (int j = 0; j < gradient.length; j++) {
                            // ∂_i A_j =  ∂_i(∑_n (g_nj * A_n / h_n))
                            final int dimi = i; // for lambda
                            final int dimj = j; // for lambda
                            DoubleUnaryOperator func = pos -> {
                                double[] newPosition = Maths.updatePosition(position, dimi, pos);
                                double sum = 0;
                                for (int n = 0; n < newPosition.length; n++) {
                                    double An = fieldWithCovariantComponents.apply(newPosition)[n];
                                    double gnj = metricCoefficient(newPosition, TensorIndexType.COVARIANT, n, dimj);
                                    double hn = Math.sqrt(metricCoefficient(newPosition, TensorIndexType.COVARIANT, n, n));
                                    sum += gnj * An / hn;
                                }
                                return sum;
                            };
                            double diSngnjAnhn = Maths.derivative(func, 1, position[i]);

                            // ∑_k Γ^k_ij A_k = ∑_k Γ^k_ij (∑_n (g_nk * A_n / h_n))
                            double GkijAk = 0;
                            for (int k = 0; k < g.length; k++) {
                                double Gkij = christoffelSymbol2ndKind(position, k, i, j);
                                double sum = 0;
                                for (int n = 0; n < g.length; n++) {
                                    sum += g[n][k] * vector[n] / Math.sqrt(g[n][n]);
                                }
                                GkijAk += Gkij * sum;
                            }
                            gpigqjTij += gContra[p][i] * gContra[q][j] * (diSngnjAnhn - GkijAk);
                        }
                    }
                    double hphq = Math.sqrt(g[p][p]*g[q][q]);
                    gradient[p][q] = hphq * gpigqjTij;
                }
            }
        }

        return gradient;
    }

    /**
     * Gets the component of the Riemann curvature tensor with the specified
     * indices.
     * 
     * @param position the position at which the Riemann tensor should be
     *                 calculated, not null
     * @param mu       the upper index
     * @param nu       the first lower index
     * @param rho      the second lower index
     * @param sigma    the third lower index
     * @return the Riemann tensor <i>R</i><sup><i>µ</i></sup><sub><i>νρσ</i></sub>
     */
    default double riemannTensor(double[] position, int mu, int nu, int rho, int sigma) {
        Objects.requireNonNull(position);

        int dim = position.length;
        if (mu < 0 || nu < 0 || rho < 0 || sigma < 0 || mu >= dim || nu >= dim || rho >= dim || sigma >= dim) {
            throw new IllegalArgumentException(
                    "mu, nu, rho and sigma may not be <0 or exceed %d, but they were mu=%d, nu=%d, rho=%d and sigma=%d"
                            .formatted(dim, mu, nu, rho, sigma));
        }

        // Γ^µ_νσ
        DoubleUnaryOperator Gmunus = pos -> {
            double[] newPosition = Maths.updatePosition(position, rho, pos);
            return christoffelSymbol2ndKind(newPosition, mu, nu, sigma);
        };

        // ∂_ρ Γ^µ_νσ
        double drGmunus = Maths.derivative(Gmunus, 1, position[rho]);

        // Γ^µ_νρ
        DoubleUnaryOperator Gmunur = pos -> {
            double[] newPosition = Maths.updatePosition(position, sigma, pos);
            return christoffelSymbol2ndKind(newPosition, mu, nu, rho);
        };

        // ∂_σ Γ^µ_νρ
        double dsGmunur = Maths.derivative(Gmunur, 1, position[sigma]);

        double[] summands = new double[2 + 2 * position.length];
        summands[0] = drGmunus;
        summands[1] = -dsGmunur;

        int counter = 2;
        for (int i = 0; i < position.length; i++) {
            // Γ^µ_ρi Γ^i_νσ
            summands[counter++] = christoffelSymbol2ndKind(position, mu, rho, i) * christoffelSymbol2ndKind(position, i, nu, sigma);

            // -Γ^µ_σi Γ^i_νρ
            summands[counter++] = -christoffelSymbol2ndKind(position, mu, sigma, i) * christoffelSymbol2ndKind(position, i, nu, rho);
        }

        return Maths.compensatedSum(summands);
    }

    /**
     * Gets whether the coordinate system is locally flat, i.e. the
     * {@link #riemannTensor(double[], int, int, int, int)} vanishes everywhere.
     * 
     * @return true if the coordinate system is locally flat
     */
    public boolean isFlat();

    /**
     * Gets the component of the Ricci tensor with the specified indices. It roughly
     * can be seen as a measure of how much the current coordinate system deviates
     * (locally) from Euclidean space.
     * 
     * @param position the position at which the Ricci tensor should be calculated,
     *                 not null
     * @param mu       the first index
     * @param nu       the second index
     * @return the Ricci tensor <i>R</i><sub><i>µν</i></sub>
     */
    default double ricciTensor(double[] position, int mu, int nu) {
        Objects.requireNonNull(position);

        double[] summands = new double[position.length];
        for (int rho = 0; rho < summands.length; rho++) {
            summands[rho] = riemannTensor(position, rho, mu, rho, nu);
        }
        return Maths.compensatedSum(summands);
    }

    /**
     * Gets the Ricci scalar, also known as scalar curvature, representing the
     * deviation that a small geodesic ball at the given position has from a
     * corresponding ball in Euclidean space.
     * 
     * @param position the position at which the Ricci scalar should be calculated,
     *                 not null
     * @return the Ricci scalar <i>R</i>
     */
    default double ricciScalar(double[] position) {
        Objects.requireNonNull(position);

        double[] summands = new double[position.length * position.length];
        for (int mu = 0; mu < summands.length; mu++) {
            for (int nu = 0; nu < summands.length; nu++) {
                double gmunu = metricCoefficient(position, TensorIndexType.CONTRAVARIANT, mu, nu);
                summands[(mu+1)*nu] = gmunu * ricciTensor(position, mu, nu);
            }
        }
        return Maths.compensatedSum(summands);
    }

    /**
     * Calculates the curl of the given vector field at the given position.
     * 
     * @param position          the position at which the curl should be calculated,
     *                          not null
     * @param componentBehavior the transformation behavior of the vector field
     *                          components, not null
     * @param vectorfield       the vector field, not null
     * @return the curl with contravariant components
     */
    // It would be good to provide a method that also handles matrix fields
    default double[] curl(double[] position, TensorTransformation componentBehavior, UnaryOperator<double[]> vectorfield) {
        Objects.requireNonNull(position);
        Objects.requireNonNull(componentBehavior);
        Objects.requireNonNull(vectorfield);

        if (position.length != 3) {
            throw new UnsupportedOperationException("Currently only implemented for 3D. "
                    + "See https://math.stackexchange.com/questions/337971/can-the-curl-operator-be-generalized-to-non-3d "
                    + "for a potential generalization.");
        }

        var covariantField = componentBehavior.transform(this, vectorfield, TensorIndexType.COVARIANT);
        // C^k = ε^ijk / sqrt(g) * (∂_i A_j)
        double[] curl = new double[position.length];
        for (int i = 0; i < curl.length; i++) {
            for (int j = 0; j < curl.length; j++) {
                for (int k = 0; k < curl.length; k++) {
                    // ∂_i A_j =  ∂_i(∑_n (g_nj * A_n / h_n))
                    final int dimi = i; // for lambda
                    final int dimj = j; // for lambda
                    DoubleUnaryOperator func = pos -> {
                        double[] newPosition = Maths.updatePosition(position, dimi, pos);
                        double sum = 0;
                        for (int n = 0; n < newPosition.length; n++) {
                            double An = covariantField.apply(newPosition)[n];
                            double gnj = metricCoefficient(newPosition, TensorIndexType.COVARIANT, n, dimj);
                            double hn = Math.sqrt(metricCoefficient(newPosition, TensorIndexType.COVARIANT, n, n));
                            sum += gnj * An / hn;
                        }
                        return sum;
                    };
                    double diAj = Maths.derivative(func, 1, position[i]);

                    curl[k] += Maths.leviCivita(i+1,j+1,k+1) * diAj;
                }
            }
        }

        double sqrtg = jacobianDeterminant(position);
        for (int i = 0; i < curl.length; i++) {
            curl[i] /= sqrtg;
        }

        return curl;
    }

    /**
     * Gets the Christoffel symbol (or coefficient of connection)
     * <i>&Gamma;</i><sup><i>m</i></sup><sub><i>ij</i></sub>=<i>e</i><sup><i>m</i></sup>
     * ∂<i>e</i><sub><i>i</i></sub>/∂<i>x</i><sup><i>j</i></sup> of the second kind
     * at the specified position (cf. Arfken, 2005, p. 155, eq. (2.138)). Note that
     * they are <em>not</em> components of a third-rank tensor. They characterize
     * how the base vectors vary in space.
     * 
     * @param position the position at which to get the Christoffel symbol, not null
     * @param m        the index of the first tangent direction (that is not
     *                 derived)
     * @param i        the index of the second tangent direction (that is partially
     *                 derived with respect to <i>j</i>)
     * @param j        the "derivative" index
     * @return the Christoffel symbol of the second kind
     */
    default double christoffelSymbol2ndKind(double[] position, int m, int i, int j) {
        Objects.requireNonNull(position);

        int dimension = position.length;
        double ret = 0;
        for (int k = 0; k < dimension; k++) {
            ret += metricCoefficient(position, TensorIndexType.CONTRAVARIANT, m, k)
                    * (metricCoefficientPartialDerivative(position, TensorIndexType.COVARIANT, i, k, j)
                            + metricCoefficientPartialDerivative(position, TensorIndexType.COVARIANT, j, k, i)
                            - metricCoefficientPartialDerivative(position, TensorIndexType.COVARIANT, i, j, k));
        }
        return 0.5 * ret;
    }

    /**
     * Gets the Christoffel symbol
     * <i>&Gamma;</i><sub><i>ijk</i></sub>=[<i>ij</i>,<i>k</i>]=<i>g</i><sub><i>mk</i></sub><i>&Gamma;</i><sup><i>m</i></sup><sub><i>ij</i></sub>
     * of the first kind at the specified position. As Christoffel symbols of the
     * second kind, Christoffel symbols of the first kind characterize how the base
     * vectors vary in space, but in metric-corrected coordinates.
     * 
     * @param position the position at which to get the Christoffel symbol, not null
     * @param i        the first index
     * @param j        the second index
     * @param k        the third index
     * @return the Christoffel symbol of the first kind
     */
    default double christoffelSymbol1stKind(double[] position, int i, int j, int k) {
        Objects.requireNonNull(position);

        double ret = 0;
        for (int m = 0; m < position.length; m++) {
            ret += metricCoefficient(position, TensorIndexType.COVARIANT, m, k) * christoffelSymbol2ndKind(position, m, i, j);
        }
        return ret;
    }

    /**
     * Gets the partial derivative of the metric coefficient, i.e.
     * ∂<i>g</i><sub><i>ij</i></sub>/∂<i>q</i><sup><i>k</i></sup> or
     * ∂<i>g</i><sup><i>ij</i></sup>/∂<i>q</i><sub><i>k</i></sub>.
     * 
     * @param position the position at which to get the derivative of the metric
     *                 coefficient
     * @param behavior the behavior of the metric coefficient with respect to a
     *                 basis change behavior, i.e. either
     *                 {@link TensorIndexType#COVARIANT covariant} or
     *                 {@link TensorIndexType#CONTRAVARIANT contravariant}
     * @param i        the row index in the metric tensor
     * @param j        the column index in the metric tensor
     * @param k        the dimension (i.e. the index of the position) to partially
     *                 derive
     * @return the partial derivative of the metric coefficient
     *         <i>g</i><sub><i>ij</i></sub> (or <i>g</i><sup><i>ij</i></sup>) with
     *         respect to the <i>k</i>th dimension
     */
    private double metricCoefficientPartialDerivative(double[] position, TensorIndexType behavior, int i, int j, int k) {
        return Maths.derivative(x -> metricCoefficient(Maths.updatePosition(position, k, x), behavior, i, j), 1, position[k]);
    }
}
