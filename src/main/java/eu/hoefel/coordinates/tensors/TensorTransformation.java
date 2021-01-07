package eu.hoefel.coordinates.tensors;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import eu.hoefel.coordinates.CoordinateSystem;
import eu.hoefel.utils.Maths;
import eu.hoefel.utils.Types;

/**
 * Interface for handling tensor index transformations, i.e.
 * {@link TensorIndexType#COVARIANT covariant},
 * {@link TensorIndexType#CONTRAVARIANT contravariant} and mixed tensors
 * can be transformed to any other index representation.
 * 
 * @author Udo Hoefel
 */
public sealed interface TensorTransformation permits TensorIndexType, MixedTensorTransformation {

	/**
	 * Creates a new tensor transformation with the specified index transformation
	 * properties. The order matters.
	 * 
	 * @param properties the index transformation properties
	 * @return a new immutable TensorTransformation
	 */
	public static TensorTransformation with(TensorIndexType... properties) {
		if (properties.length == 1) return properties[0];

		return new MixedTensorTransformation(List.of(properties));
	}

	/**
	 * Gets the tensor index transformation properties.
	 * 
	 * @return the transformation properties of the indices. For the special cases
	 *         with respect to the length of the list see above.
	 */
	public List<TensorIndexType> indexTypes();

	/**
	 * Changes the tensor index transformation behavior of the given tensor field
	 * from the current tensor transformation to the specified one, i.e. it lowers
	 * and raises all the necessary indices. This works for tensors of any rank.
	 * Note that it is valid to have either target transformation properties of
	 * length one even if the current tensor transformation has more than one index
	 * or the current tensor transformation properties may have only one
	 * transformation property but the target may have more than one index. In both
	 * these cases the single value tensor index transformation property will be
	 * extended to as many indices as necessary.
	 * 
	 * @param <T>            the return type of the tensor field, i.e. it needs to
	 *                       return either a {@code double[]} (corresponding to a
	 *                       tensor field of rank 1) or any higher dimensional array
	 *                       thereof.
	 * @param sys            the coordinate system
	 * @param tensorfield    the tensor field of any rank
	 * @param targetBehavior the desired tensor index transformation behavior
	 * @return the tensor field with the appropriately raised and lowered indices
	 */
	@SuppressWarnings("unchecked")
	default <T> Function<double[], T> transform(CoordinateSystem sys, Function<double[], T> tensorfield,
			TensorTransformation targetBehavior) {
		return pos -> {
			T tensor = tensorfield.apply(pos);
			int actualRank = Types.dimension(tensor.getClass());

			List<TensorIndexType> transformationProperties = indexTypes();
			int rank = indexTypes().size();

			List<TensorIndexType> targetTransformationProperties = targetBehavior.indexTypes();
			int targetRank = targetTransformationProperties.size();
			
			// handle special cases
			if (rank == 1 && targetRank > 1 && actualRank == targetRank) {
				TensorIndexType[] tps = new TensorIndexType[targetRank];
				Arrays.fill(tps, transformationProperties.get(0));
				transformationProperties = List.of(tps);
				rank = targetRank;
			} else if (rank > 1 && targetRank == 1 && actualRank == rank) {
				TensorIndexType[] tps = new TensorIndexType[rank];
				Arrays.fill(tps, targetTransformationProperties.get(0));
				targetTransformationProperties = List.of(tps);
				targetRank = rank;
			} else if (rank == 1 && targetRank == 1 && actualRank > 1) {
				TensorIndexType[]  tps = new TensorIndexType[actualRank];
				Arrays.fill(tps, transformationProperties.get(0));
				transformationProperties = List.of(tps);
				rank = actualRank;
				
				tps = new TensorIndexType[actualRank];
				Arrays.fill(tps, targetTransformationProperties.get(0));
				targetTransformationProperties = List.of(tps);
				targetRank = actualRank;
			} else if (rank == actualRank && actualRank != targetRank) {
				throw new IllegalArgumentException(
						"Cannot transform tensor T%s to tensor T%s ('d' represents a covariant, 'u' a contravariant component)"
								.formatted(indexPositionRepresentation(transformationProperties.stream()),
										indexPositionRepresentation(targetTransformationProperties.stream())));
			} else if (rank != actualRank && actualRank != targetRank) {
				TensorIndexType[]  tps = new TensorIndexType[actualRank];
				Arrays.fill(tps, transformationProperties.get(0));
				transformationProperties = List.of(tps);
				throw new IllegalArgumentException(
						"Cannot transform tensor T%s to tensor T%s ('d' represents a covariant, 'u' a contravariant component)"
								.formatted(indexPositionRepresentation(transformationProperties.stream()),
										indexPositionRepresentation(targetTransformationProperties.stream())));
			}
			
			// special case for performance of vector field calculations (probably by far
			// most used)
			if (tensor instanceof double[] vector) {
				if (this == targetBehavior) return (T) vector;
				return (T) switchVectorIndexPosition(sys, pos, vector, transformationProperties.get(0).flip());
			}

			int dim = Array.getLength(tensor);
			double[][] metricTensorCovariant = sys.metricTensor(pos, TensorIndexType.COVARIANT);
			double[][] metricTensorContravariant = sys.metricTensor(pos, TensorIndexType.CONTRAVARIANT);
			
			for (int i = 0; i < targetRank; i++) {
				TensorIndexType currentTransformation = transformationProperties.get(i);
				if (targetTransformationProperties.get(i) == currentTransformation) continue;

				double[][] g;
				if (currentTransformation == TensorIndexType.COVARIANT) {
					g = metricTensorContravariant;
				} else {
					g = metricTensorCovariant;
				}
				
				T preCurrentTransformationTensor = Maths.deepCopyPrimitiveArray(Types.unbox(tensor));
				
				int numEvals = (int) Math.pow(dim, targetRank - 1);
				Supplier<int[]> indices = tensorIndices(targetRank, dim, i);
				for (int dummy = 0; dummy < numEvals; dummy++) {
					int[] idx = indices.get(); // this updates each loop iteration
					double[] tensorValues = getTensorValues(preCurrentTransformationTensor, targetRank, dim, i, idx);
					double[] result = new double[dim];
					for (int j = 0; j < dim; j++) {
						for (int k = 0; k < dim; k++) {
							result[j] += g[j][k] * tensorValues[k];
						}
					}
					setTensorValues(tensor, result, i, idx);
				}
			}
			
			return tensor;
		};
	}

	/**
	 * Switches the index position of the vector (corresponding to a tensor of rank
	 * 1), i.e. if the vector currently has {@link TensorIndexType#COVARIANT
	 * covariant} components they will be returned as
	 * {@link TensorIndexType#CONTRAVARIANT contravariant} components and
	 * vice versa.
	 * 
	 * @param sys                    the coordinate system
	 * @param position               the position
	 * @param vector                 the vector whose index should switch position
	 * @param transformationProperty the target index position, i.e. covariant for a
	 *                               lower index and contravariant for an upper
	 *                               index
	 * @return the vector with covariant components if they were given as
	 *         contravariant components and vice versa
	 */
	private static double[] switchVectorIndexPosition(CoordinateSystem sys, double[] position, double[] vector, TensorIndexType transformationProperty) {
		int dim = vector.length;
        double[][] g = sys.metricTensor(position, transformationProperty);
		double[] ret = new double[dim];
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				ret[i] += g[i][j] * vector[j];
			}
		}
        return ret;
	}

	/**
	 * Creates a simple string representation of the index positions, "d" represents
	 * a covariant (i.e. lower) and "u" a contravariant (i.e. upper) index.
	 * 
	 * @param s the stream of transformation properties to represent
	 * @return the representation, e.g. "ddu" for a tensor whose first two indices
	 *         behave covariantly and the last index behaves contravariantly
	 */
	private static String indexPositionRepresentation(Stream<TensorIndexType> s) {
		return s.map(tp -> tp == TensorIndexType.COVARIANT ? "d" : "u").collect(Collectors.joining(""));
	}

	/**
	 * Provides an index supplier. Note that the supplier varies the returned
	 * indices whenever {@link Supplier#get()} is called. The indices are varied,
	 * one at a time, up to the given {@code dimension} (exclusive). The only index
	 * in the returned array that is not varied is the one specified in
	 * {@code stableIndex}.
	 * 
	 * @param rank        the rank of tensor for which the indices shall be created
	 * @param dimension   the dimension of the tensor
	 * @param stableIndex the index of the indices that should not vary
	 * @return an index supplier for a tensor of the specified rank and dimension.
	 *         The supplier will iterate subsequently through all possible indices
	 *         save the given {@code stableIndex}.
	 */
	private static Supplier<int[]> tensorIndices(int rank, int dimension, int stableIndex) {
		int[] indices = new int[rank];
		boolean[] zeroAlreadyUsed = new boolean[rank];
		return () -> {
			for (int i = indices.length - 1; i >= 0; i--) {
				if (indices[i] == dimension - 1 || i == stableIndex) {
					continue;
				} else if (!zeroAlreadyUsed[i]) {
					zeroAlreadyUsed[i] = true;
					continue;
				}
				indices[i]++;
				break;
			}
			return indices;
		};
	}

	/**
	 * Gets the values of the arbitrarily ranked tensor in the following way:<br>
	 * For example, <code>indices={1,2,3}</code>, <code>dimension=3</code> and
	 * <code>index=1</code> would return the tensor values at the indices
	 * <code>[1][0,1,2][3]</code>.
	 * 
	 * @param <R>       the double array type, i.e. it needs to be a double[] or
	 *                  higher dimensional array
	 * @param array     the tensor, i.e. a double[] array for a tensor of rank 1,
	 *                  double[][] for a tensor of order 2 and so on
	 * @param rank      the rank (i.e., the dimensionality of the array, e.g. rank 2
	 *                  would correspond to a double[][] array)
	 * @param dimension the length of the array
	 * @param index     the index of which the values are desired
	 * @param indices   the indices to fetch (except that the specified index in the
	 *                  indices will be varied)
	 * @return the values of the tensor for the fixed, specified indices and the
	 *         given, variable index in the indices
	 */
	private static <R> double[] getTensorValues(R array, int rank, int dimension, int index, int[] indices) {
		List<Double> values = new ArrayList<>();
		int[] summandIndices = indices.clone();
		for (int i = 0; i < dimension; i++) {
			summandIndices[index] = i;
			values.add(getTensorValue(array, summandIndices));
		}
		return values.stream().mapToDouble(f->f).toArray();
	}

	/**
	 * Gets the value of the array at the given indices.
	 * 
	 * @param <R>     the double array type, i.e. it needs to be a double[] or
	 *                higher dimensional array
	 * @param array   the tensor, i.e. a double[] array for a tensor of rank 1,
	 *                double[][] for a tensor of order 2 and so on
	 * @param indices the indices to fetch
	 * @return the value at the given indices
	 */
	private static <R> double getTensorValue(R array, int... indices) {
		Object o = array;
		for (int i = 0; i < indices.length - 1; i++) {
			o = Array.get(o, indices[i]);
		}
		return Array.getDouble(o, indices[indices.length - 1]);
	}

	/**
	 * Sets the values of the arbitrarily ranked tensor in the following way:<br>
	 * For example, <code>indices={1,2,3}</code>, <code>dimension=3</code> and
	 * <code>index=1</code> would set the values at the indices
	 * <code>[1][0,1,2][3]</code>.
	 * 
	 * @param <R>     the double array type, i.e. it needs to be a double[] or
	 *                higher dimensional array
	 * @param array   the tensor, i.e. a double[] array for a tensor of rank 1,
	 *                double[][] for a tensor of order 2 and so on
	 * @param values  the values to store
	 * @param index   the index in which the values are to be stored
	 * @param indices the indices to store the values to (the specified index in the
	 *                indices will be varied and hence a double[] can be stored)
	 */
	private static <R> void setTensorValues(R array, double[] values, int index, int[] indices) {
		int[] summandIndices = indices.clone();
		for (int i = 0; i < values.length; i++) {
			summandIndices[index] = i;
			setTensorValue(array, values[i], summandIndices);
		}
	}

	/**
	 * Sets the value of the array at the given indices.
	 * 
	 * @param <R>     the double array type, i.e. it needs to be a double[] or
	 *                higher dimensional array
	 * @param array   the tensor, i.e. a double[] array for a tensor of rank 1,
	 *                double[][] for a tensor of order 2 and so on
	 * @param value   the value to set
	 * @param indices the indices for which to set the value
	 */
	private static <R> void setTensorValue(R array, double value, int... indices) {
		Object o = array;
		for (int i = 0; i < indices.length - 1; i++) {
			o = Array.get(o, indices[i]);
		}
		if (Types.elementType(array.getClass()).isPrimitive()) {
			Array.setDouble(o, indices[indices.length - 1], value);
		} else {
			Array.set(o, indices[indices.length - 1], value);
		}
	}
}
