package eu.hoefel.coordinates.tensors;

import java.util.List;

/**
 * Contains the transformation properties of the indices of a tensor of
 * arbitrary rank.
 * 
 * @param indexTypes the transformation properties of the indices of a tensor
 * 
 * @author Udo Hoefel
 */
record MixedTensorTransformation(List<TensorIndexType> indexTypes) implements TensorTransformation {}
