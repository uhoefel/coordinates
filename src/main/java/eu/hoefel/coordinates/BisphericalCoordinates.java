package eu.hoefel.coordinates;

import java.util.Collections;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Objects;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.function.UnaryOperator;

import eu.hoefel.coordinates.axes.Axes;
import eu.hoefel.coordinates.axes.Axis;
import eu.hoefel.coordinates.tensors.TensorIndexType;
import eu.hoefel.coordinates.tensors.TensorTransformation;
import eu.hoefel.unit.Unit;
import eu.hoefel.unit.Units;
import eu.hoefel.unit.constant.Constant;
import eu.hoefel.unit.si.SiDerivedUnit;
import eu.hoefel.utils.Maths;

/**
 * Bispherical coordinates are obtained from {@link BipolarCoordinates} by
 * rotating about the axis connecting its foci. Consequently, the foci separated
 * by 2<i>a</i> in bipolar coordinates are still points in bispherical
 * coordinates.
 * <p>
 * They are defined as follows:<br>
 * {@code x = a sin(σ)cos(φ)/(cosh(τ)-cos(σ))}<br>
 * {@code y = a sin(σ)sin(φ)/(cosh(τ)-cos(σ))}<br>
 * {@code z = a sinh(τ)/(cosh(τ)-cos(σ))}<br>
 * <p>
 * Note the similarity to {@link ToroidalCoordinates}, which are created by
 * rotating about the foci-separating axis.
 * 
 * @param axes the axes defining the coordinate system, not null
 * @param a    the distance of the foci to the origin, needs to be larger than 0
 * 
 * @author Udo Hoefel
 */
@CoordinateSystemSymbols({"bispherical", "bisph"})
public final record BisphericalCoordinates(NavigableSet<Axis> axes, Constant a) implements CoordinateSystem {

    /** The default axes. */
    public static final NavigableSet<Axis> DEFAULT_AXES = Axes.of(
            new Axis(0, SiDerivedUnit.RADIAN, "σ"),
            new Axis(1, SiDerivedUnit.RADIAN, "τ"),
            new Axis(2, SiDerivedUnit.RADIAN, "φ"));

    /** Constructs a new bispherical coordinate system. */
    public BisphericalCoordinates {
        Objects.requireNonNull(axes, "Axes may not be null. "
                + "Use the DEFAULT_AXES or the constructor that just requires ['a'] to get a reasonable default.");

        if (!Units.convertible(Axis.fromSet(axes, 0).unit(), Units.EMPTY_UNIT)) {
            throw new IllegalArgumentException("The unit of dimension 0 (%s) needs to be effectively dimensionless."
                    .formatted(Axis.fromSet(axes, 0).unit()));
        }

        if (!Units.convertible(Axis.fromSet(axes, 1).unit(), Units.EMPTY_UNIT)) {
            throw new IllegalArgumentException("The unit of dimension 1 (%s) needs to be effectively dimensionless."
                    .formatted(Axis.fromSet(axes, 1).unit()));
        }

        if (!Units.convertible(Axis.fromSet(axes, 2).unit(), Units.EMPTY_UNIT)) {
            throw new IllegalArgumentException("The unit of dimension 2 (%s) needs to be effectively dimensionless."
                    .formatted(Axis.fromSet(axes, 2).unit()));
        }
    }

    /**
     * Constructs a new bispherical coordinate system.
     * 
     * @param args the arguments, must be either {@link Axes} or {@link Axis}, which
     *               will take precedence over the {@link #DEFAULT_AXES} if given. If
     *               no arguments are given, the default axes will be used.
     *               Further required arguments may also be passed in here, but (if 
     *               they are doubles or {@link Constant}) they have to be in the 
     *               same order in which they are specified in the record declaration.
     */
    public BisphericalCoordinates(Object... args) {
        this(Axes.fromArgs(DEFAULT_AXES, args),
                CoordinateSystems.constantFromArgs(0, args)
                    .orElseThrow(() -> new IllegalArgumentException("No Constant found!")));
    }

    /**
     * Constructs a new bispherical coordinate system.
     * 
     * @param a the distance of the foci to the origin, needs to be larger than 0
     */
    public BisphericalCoordinates(Constant a) {
        this(DEFAULT_AXES, a);
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
        map.put(0, a.unit());
        map.put(1, a.unit());
        map.put(2, a.unit());
        return Collections.unmodifiableNavigableMap(map);
    }

    @Override
    public double[] toBasePoint(double[] position) {
        validatePosition(position);
        return Transformations.bisphericalToCartesian(position, a.value());
    }

    @Override
    public double[] fromBasePoint(double[] position) {
        double[] pointInCurrentSystem = new double[3];
        pointInCurrentSystem[0] = Math.acos(Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1.0/2)*(Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(a.value(),2)));
        pointInCurrentSystem[1] = Maths.arsinh(2*position[2]*a.value()*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1.0/2));
        pointInCurrentSystem[2] = Math.atan2(position[1],position[0]);
        return pointInCurrentSystem;
    }

    @Override
    public double[] toBaseVector(double[] position, double[] vector) {
        validatePosition(position);

        double[] vectorInBaseSys = new double[vector.length];
        vectorInBaseSys[0] = a.value()*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-2)*(-1+Math.cos(position[0])*Math.cosh(position[1]))*Math.cos(position[2])*vector[0]+a.value()*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-2)*(-1+Math.cos(position[0])*Math.cosh(position[1]))*Math.sin(position[2])*vector[1]+-1*a.value()*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-2)*Math.sin(position[0])*Math.sinh(position[1])*vector[2];
        vectorInBaseSys[1] = -1*a.value()*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-2)*Math.cos(position[2])*Math.sin(position[0])*Math.sinh(position[1])*vector[0]+-1*a.value()*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-2)*Math.sin(position[2])*Math.sin(position[0])*Math.sinh(position[1])*vector[1]+-1*a.value()*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-2)*(-1+Math.cos(position[0])*Math.cosh(position[1]))*vector[2];
        vectorInBaseSys[2] = a.value()*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*Math.sin(position[2])*Math.sin(position[0])*vector[0]+-1*a.value()*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*Math.cos(position[2])*Math.sin(position[0])*vector[1];
        return vectorInBaseSys;
    }

    @Override
    public double[] fromBaseVector(double[] position, double[] vector) {
        double[] vectorInCurrentSys = new double[vector.length];
        vectorInCurrentSys[0] = -1*Math.pow((1-1*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(a.value(),2)),2)),-1.0/2)*(2*position[0]*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1.0/2)-2*position[0]*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-3.0/2)*(Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2))*(Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(a.value(),2)))*vector[0]-4*position[0]*position[2]*a.value()*Math.pow((1+4*Math.pow(position[2],2)*Math.pow(a.value(),2)*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1)),-1.0/2)*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-3.0/2)*(Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2))*vector[1]-1*position[1]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)),-1)*vector[2];
        vectorInCurrentSys[1] = -1*Math.pow((1-1*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(a.value(),2)),2)),-1.0/2)*(2*position[1]*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1.0/2)-2*position[1]*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-3.0/2)*(Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2))*(Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(a.value(),2)))*vector[0]-4*position[1]*position[2]*a.value()*Math.pow((1+4*Math.pow(position[2],2)*Math.pow(a.value(),2)*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1)),-1.0/2)*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-3.0/2)*(Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2))*vector[1]+position[0]*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)),-1)*vector[2];
        vectorInCurrentSys[2] = -1*Math.pow((1-1*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1)*Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(a.value(),2)),2)),-1.0/2)*(2*position[2]*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1.0/2)+Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-3.0/2)*(-2*position[2]*(Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2))+4*position[2]*Math.pow(a.value(),2))*(Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)-1*Math.pow(a.value(),2)))*vector[0]+Math.pow((1+4*Math.pow(position[2],2)*Math.pow(a.value(),2)*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1)),-1.0/2)*(2*a.value()*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-1.0/2)+2*position[2]*a.value()*Math.pow((Math.pow((Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2)),2)-4*Math.pow(position[2],2)*Math.pow(a.value(),2)),-3.0/2)*(-2*position[2]*(Math.pow(position[0],2)+Math.pow(position[1],2)+Math.pow(position[2],2)+Math.pow(a.value(),2))+4*position[2]*Math.pow(a.value(),2)))*vector[1];
        return vectorInCurrentSys;
    }

    @Override
    public double metricCoefficient(double[] position, TensorTransformation behavior, int i, int j) {
        validatePosition(position);

        if (i < 0 || j < 0) {
            throw new IllegalArgumentException(("Metric coefficient not available for i=%d, j=%d "
                    + "(too low dimension, only 3 dimensions [0,1,2] are supported for bispherical coordinates)")
                    .formatted(i, j));
        } else if (i > 2 || j > 2) {
            throw new IllegalArgumentException(("Metric coefficient not available for i=%d, j=%d "
                    + "(too high dimension, only 3 dimensions [0,1,2] are supported for bispherical coordinates)")
                    .formatted(i, j));
        }

        if (behavior instanceof TensorIndexType tit) {
            return switch (tit) {
                case COVARIANT -> {
                    if (i == 0 && j == 0) yield Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-4)*(Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2));
                    if (i == 1 && j == 1) yield Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-4)*(Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2));
                    if (i == 2 && j == 2) yield Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-2)*Math.pow(Math.sin(position[0]),2);
                    yield 0;
                }
                case CONTRAVARIANT -> {
                    if (i == 0 && j == 0) yield Math.pow(a.value(),-2)*Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),4);
                    if (i == 1 && j == 1) yield Math.pow(a.value(),-2)*Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),4);
                    if (i == 2 && j == 2) yield Math.pow(a.value(),-2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)*Math.pow(Math.sin(position[0]),-2);
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
        return Math.pow(a.value(),2)*Math.pow(Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-10)*(1+Math.pow(Math.cos(position[0]),4)*Math.pow(Math.cosh(position[1]),4)+Math.pow(Math.sin(position[0]),4)*Math.pow(Math.sinh(position[1]),4)-4*Math.pow(Math.cos(position[0]),3)*Math.pow(Math.cosh(position[1]),3)-4*Math.cos(position[0])*Math.cosh(position[1])+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)+6*Math.pow(Math.cos(position[0]),2)*Math.pow(Math.cosh(position[1]),2)-1.0/32*(-1+Math.cos(4*position[0]))*(-1+Math.cosh(4*position[1]))-4*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)*Math.cos(position[0])*Math.cosh(position[1])),0.5)*Math.abs(a.value()*Math.sin(position[0]));
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

        if (i == 0 && j == 0 && k == 0) return Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-5)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sin(position[0]);
        if (i == 0 && j == 0 && k == 1) return -1*Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-5)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sinh(position[1]);
        if (i == 0 && j == 1 && k == 0) return Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-5)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sinh(position[1]);
        if (i == 0 && j == 1 && k == 1) return Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-5)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sin(position[0]);
        if (i == 0 && j == 2 && k == 2) return -1*Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-3)*(-1+Math.cos(position[0])*Math.cosh(position[1]))*Math.sin(position[0]);
        if (i == 1 && j == 0 && k == 0) return Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-5)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sinh(position[1]);
        if (i == 1 && j == 0 && k == 1) return Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-5)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sin(position[0]);
        if (i == 1 && j == 1 && k == 0) return -1*Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-5)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sin(position[0]);
        if (i == 1 && j == 1 && k == 1) return Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-5)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sinh(position[1]);
        if (i == 1 && j == 2 && k == 2) return Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-3)*Math.pow(Math.sin(position[0]),2)*Math.sinh(position[1]);
        if (i == 2 && j == 0 && k == 2) return -1*Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-3)*(-1+Math.cos(position[0])*Math.cosh(position[1]))*Math.sin(position[0]);
        if (i == 2 && j == 1 && k == 2) return Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-3)*Math.pow(Math.sin(position[0]),2)*Math.sinh(position[1]);
        if (i == 2 && j == 2 && k == 0) return Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-3)*(-1+Math.cos(position[0])*Math.cosh(position[1]))*Math.sin(position[0]);
        if (i == 2 && j == 2 && k == 1) return -1*Math.pow(a.value(),2)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-3)*Math.pow(Math.sin(position[0]),2)*Math.sinh(position[1]);
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

        if (m == 0 && i == 0 && j == 0) return Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sin(position[0]);
        if (m == 0 && i == 0 && j == 1) return Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sinh(position[1]);
        if (m == 0 && i == 1 && j == 0) return Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sinh(position[1]);
        if (m == 0 && i == 1 && j == 1) return -1*Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sin(position[0]);
        if (m == 0 && i == 2 && j == 2) return Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*(-1+Math.cos(position[0])*Math.cosh(position[1]))*(-1*Math.cosh(position[1])+Math.cos(position[0]))*Math.sin(position[0]);
        if (m == 1 && i == 0 && j == 0) return -1*Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sinh(position[1]);
        if (m == 1 && i == 0 && j == 1) return Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sin(position[0]);
        if (m == 1 && i == 1 && j == 0) return Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sin(position[0]);
        if (m == 1 && i == 1 && j == 1) return Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*(-1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),2)+2*Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+2*Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2))*Math.sinh(position[1]);
        if (m == 1 && i == 2 && j == 2) return -1*Math.pow((Math.pow((-1+Math.cos(position[0])*Math.cosh(position[1])),2)+Math.pow(Math.sin(position[0]),2)*Math.pow(Math.sinh(position[1]),2)),-1)*Math.pow(Math.sin(position[0]),2)*(-1*Math.cosh(position[1])+Math.cos(position[0]))*Math.sinh(position[1]);
        if (m == 2 && i == 0 && j == 2) return -1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*Math.pow(Math.sin(position[0]),-1)*(-1+Math.cos(position[0])*Math.cosh(position[1]));
        if (m == 2 && i == 1 && j == 2) return Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*Math.sinh(position[1]);
        if (m == 2 && i == 2 && j == 0) return -1*Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*Math.pow(Math.sin(position[0]),-1)*(-1+Math.cos(position[0])*Math.cosh(position[1]));
        if (m == 2 && i == 2 && j == 1) return Math.pow((-1*Math.cosh(position[1])+Math.cos(position[0])),-1)*Math.sinh(position[1]);
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
