package eu.hoefel.coordinates.tensors;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

import java.util.function.Function;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import eu.hoefel.coordinates.CoordinateSystem;

/**
 * Tests for everything related to tensor transformations.
 * 
 * @author Udo Hoefel
 */
@DisplayName("Testing tensor transformations")
class TensorTests {

	@DisplayName("Testing vectorfield transformation")
	@Test
	void testVectorfieldTransformation() {
		CoordinateSystem sys = CoordinateSystem.from("cart", 3);
		double[] position = {2,2,3};
		Function<double[],double[]> tensorfield = Function.identity();

		Function<double[],double[]> transformedTensorfield = TensorIndexType.COVARIANT.transform(sys, tensorfield, TensorIndexType.COVARIANT);
		assertArrayEquals(tensorfield.apply(position), transformedTensorfield.apply(position));
		
		transformedTensorfield = TensorIndexType.COVARIANT.transform(sys, tensorfield, TensorIndexType.CONTRAVARIANT);
		assertArrayEquals(tensorfield.apply(position), transformedTensorfield.apply(position));
		
		sys = CoordinateSystem.from("cyl");
		transformedTensorfield = TensorIndexType.COVARIANT.transform(sys, tensorfield, TensorIndexType.COVARIANT);
		assertArrayEquals(tensorfield.apply(position), transformedTensorfield.apply(position));

		sys = CoordinateSystem.from("cyl");
		transformedTensorfield = TensorIndexType.COVARIANT.transform(sys, tensorfield, TensorIndexType.CONTRAVARIANT);
		assertArrayEquals(new double[] {2,0.5,3}, transformedTensorfield.apply(position));
		transformedTensorfield = TensorIndexType.CONTRAVARIANT.transform(sys, tensorfield, TensorIndexType.COVARIANT);
		assertArrayEquals(new double[] {2,8,3}, transformedTensorfield.apply(position));
	}

	@DisplayName("Testing tensorfield transformation")
	@Test
	void testTensorfieldTransformation() {
		CoordinateSystem sys = CoordinateSystem.from("cart", 3);
		double[] position = {2,3,4};
		Function<double[],double[][]> tensorfield = pos -> new double[][] {{pos[0],pos[1],pos[2]},{1.5*pos[0],1.5*pos[1],1.5*pos[2]},{2*pos[0],2*pos[1],2*pos[2]}};
		Function<double[],double[][]> transformedTensorfield = TensorIndexType.COVARIANT.transform(sys, tensorfield, TensorIndexType.COVARIANT);
		assertArrayEquals(tensorfield.apply(position), transformedTensorfield.apply(position));
		
		transformedTensorfield = TensorIndexType.COVARIANT.transform(sys, tensorfield, TensorIndexType.CONTRAVARIANT);
		assertArrayEquals(tensorfield.apply(position), transformedTensorfield.apply(position));
		
		sys = CoordinateSystem.from("cyl");
		transformedTensorfield = TensorIndexType.COVARIANT.transform(sys, tensorfield, TensorIndexType.COVARIANT);
		assertArrayEquals(tensorfield.apply(position), transformedTensorfield.apply(position));

		transformedTensorfield = TensorIndexType.COVARIANT.transform(sys, tensorfield, TensorTransformation.with(TensorIndexType.COVARIANT, TensorIndexType.CONTRAVARIANT));
		assertArrayEquals(new double[][] {{2,0.75,4},{3,1.125,6},{4,1.5,8}}, transformedTensorfield.apply(position));

		transformedTensorfield = TensorIndexType.COVARIANT.transform(sys, tensorfield, TensorIndexType.CONTRAVARIANT);
		assertArrayEquals(new double[][] {{2,0.75,4},{0.75,0.28125,1.5},{4,1.5,8}}, transformedTensorfield.apply(position));
	}
}
