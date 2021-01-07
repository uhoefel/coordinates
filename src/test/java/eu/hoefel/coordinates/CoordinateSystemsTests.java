package eu.hoefel.coordinates;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Random;
import java.util.function.Function;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ArgumentsSource;

import eu.hoefel.coordinates.axes.Axis;
import eu.hoefel.coordinates.tensors.TensorIndexType;
import eu.hoefel.utils.Types;

/**
 * Tests for everything related to the coordinate systems implementations.
 * 
 * @author Udo Hoefel
 */
@DisplayName("Testing coordinate systems")
class CoordinateSystemsTests {

	@DisplayName("Testing transformation")
	@ParameterizedTest
	@ArgumentsSource(DefaultCoordinateSystemProvider.class)
	void testTransformation(CoordinateSystem coordSys) {
		double[] coords;
		if (coordSys.dimension() > 2) {
			coords = new double[] { 1, 2, 3 };
		} else {
			coords = new double[] { 1, 2 };
		}
		assertArrayEquals(coords, coordSys.toBasePoint(coordSys.fromBasePoint(coords)), 1e-12);
	}

	@DisplayName("Testing transformation for identity coordinate system")
	@Test
	void testTransformation_IdentityCoordinateSystem() {
		double[] coords = {1,2,3,4,5};
		CoordinateSystem coordSys = CoordinateSystems.IDENTITY_COORDINATE_SYSTEM;
		assertArrayEquals(coords, coordSys.toBasePoint(coordSys.fromBasePoint(coords)));

		coords = new double[] {-3,12,52,856,1.14141414};
		assertArrayEquals(coords, CoordinateSystems.transform(coords, "", ""));
	}

	@DisplayName("Testing cross product")
	@Test
	void testCrossProduct() {
		CoordinateSystem cart = CoordinateSystem.from("cart");
		double[] w1 = {0,1,0};
		double[] w2 = {1,0,12};
		assertArrayEquals(new double[] {12,0,-1}, cart.cross(new double[] {0,0,0}, TensorIndexType.COVARIANT, w1, w2));
		
		double[] v1 = {0,0,1};
		double[] v2 = {0,1,0};
		assertArrayEquals(new double[] {-1,0,0}, cart.cross(new double[] {0,0,0}, TensorIndexType.COVARIANT, v1, v2));

		double[] u1 = {0,0,0,1};
		double[] u2 = {0,0,1,0};
		double[] u3 = {0,1,0,0};
		assertArrayEquals(new double[] {-1,0,0,0}, cart.cross(new double[] {0,0,0,0}, TensorIndexType.COVARIANT, u1, u2, u3));
	}

	@DisplayName("Testing dot product")
	@Test
	void testDotProduct() {
		CoordinateSystem cart = CoordinateSystem.from("cart");
		double[] u1 = {1,2,3,7};
		double[] u2 = {4,3,2,1};
		assertEquals(23.0, cart.dot(new double[] {0,0,0,0}, TensorIndexType.COVARIANT, u1, u2));
	}

	@DisplayName("Testing vector transformation")
	@ParameterizedTest
	@ArgumentsSource(DefaultCoordinateSystemProvider.class)
	void testVectorTransformation(CoordinateSystem coordSys) {
		double[] position;
		double[] vector;
		if (coordSys.dimension() > 2) {
			position = new double[] { 0.1, 0.2, 0.3 };
			vector = new double[] { 1, 2, 3 };
		} else {
			position = new double[] { 0.1, 0.2 };
			vector = new double[] { 1, 2 };
		}

		// quite inaccurate with the default implementation
		double[] vectorInCoordSys = coordSys.fromBaseVector(coordSys.toBasePoint(position), vector);
		assertArrayEquals(vector, coordSys.toBaseVector(position, vectorInCoordSys), 9e-4);
	}

	@DisplayName("Testing divergence")
	@Test
	void testDivergence() {
		// tensor rank 1
		CoordinateSystem csys = CoordinateSystem.from("cart");
		Function<double[], Double[]> field = xyz -> new Double[] { -xyz[1], xyz[0] * xyz[1], xyz[2] };

		int numPositionsToTest = 15;
		for (int i = 0; i < numPositionsToTest; i++) {
			double x = 5 * Math.random();
			double[] pos = new double[] {x, Math.random(), Math.random()};
			assertEquals(x+1, csys.div(pos, TensorIndexType.CONTRAVARIANT, field), 2.5e-6);
		}

		field = xyz -> new Double[] { xyz[0]*xyz[0]-xyz[1]*xyz[1], 2 * xyz[0] * xyz[1] };
		for (int i = 0; i < numPositionsToTest; i++) {
			double x = Math.random();
			double[] pos = new double[] {x, Math.random()};
			assertEquals(4*x, csys.div(pos, TensorIndexType.CONTRAVARIANT, field), 1e-5);
		}

		field = xyz -> new Double[] { xyz[0], xyz[1], xyz[2] };
		for (int i = 0; i < numPositionsToTest; i++) {
			double x = Math.random();
			double[] pos = new double[] {x, Math.random(), Math.random()};
			assertEquals(3, csys.div(pos, TensorIndexType.CONTRAVARIANT, field), 1e-5);
		}

		csys = CoordinateSystem.from("cyl");
		field = xyz -> new Double[] { xyz[0], 0.0, xyz[2] };
		for (int i = 0; i < numPositionsToTest; i++) {
			double x = Math.random();
			double[] pos = new double[] {x, Math.random(), Math.random()};
			assertEquals(3, csys.div(pos, TensorIndexType.CONTRAVARIANT, field), 1e-5);
		}

		// tensor rank 2
		csys = CoordinateSystem.from("cart");
		double[] pos = new double[] {5,2.5,3};
		Function<double[], double[][]> tensor = xyz -> new double[][] { {xyz[0], xyz[1], xyz[2]}, {1.5 * xyz[0], 1.5 * xyz[1], 1.5 * xyz[2]}, {2 * xyz[0], 2 * xyz[1], 2 * xyz[2]} };
		assertArrayEquals(new double[] {3,4.5,6}, csys.div(pos, TensorIndexType.COVARIANT, tensor), 2.5e-7);
	}

	@DisplayName("Testing magnitude")
	@Test
	void testMagnitude() {
		CoordinateSystem csys = CoordinateSystem.from("cart");
		int numPositionsToTest = 15;
		
		// test scalar field
		Function<double[], Double> scalarfield = xyz -> -xyz[1];
		for (int i = 0; i < numPositionsToTest; i++) {
			double y = (Math.random() > 0.5 ? 1 : -1) * 5 * Math.random();
			double[] pos = {Math.random(), y, Math.random()};
			assertEquals(Math.abs(y), csys.magnitude(pos, TensorIndexType.CONTRAVARIANT, scalarfield));
			assertEquals(Math.abs(y), csys.magnitude(pos, TensorIndexType.COVARIANT, scalarfield));
		}

		// test vector field
		Function<double[], double[]> vectorfield = xyz -> new double[] { -xyz[1], xyz[0] * xyz[1], xyz[2] };
		for (int i = 0; i < numPositionsToTest; i++) {
			double x = 5 * Math.random();
			double y = (Math.random() > 0.5 ? 1 : -1) * Math.random();
			double z = Math.random();
			double[] pos = {x, y, z};
			assertEquals(Math.sqrt(y*y + x*y*x*y + z*z), csys.magnitude(pos, TensorIndexType.CONTRAVARIANT, vectorfield), 2.5e-7);
			assertEquals(Math.sqrt(y*y + x*y*x*y + z*z), csys.magnitude(pos, TensorIndexType.COVARIANT, vectorfield), 2.5e-7);
		}

		csys = CoordinateSystem.from("cyl");
		vectorfield = xyz -> new double[] { xyz[0], xyz[1], xyz[2] };
		for (int i = 0; i < numPositionsToTest; i++) {
			double x = Math.random();
			double y = Math.random();
			double z = (Math.random() > 0.5 ? 1 : -1) * Math.random();
			double[] pos = new double[] {x, y, z};
			// the x^2 before y^2 is due to the metric tensor coefficient at g22=r^2
			assertEquals(Math.sqrt(x*x + x*x*y*y + z*z), csys.magnitude(pos, TensorIndexType.CONTRAVARIANT, vectorfield), 1e-7);
		}

		// tensor rank 2
		Function<double[], double[][]> tensorfield = xyz -> new double[][] {{ xyz[0], xyz[1], xyz[2] }, { 1.5*xyz[0], 1.5*xyz[1], 1.5*xyz[2] }, { 2*xyz[0], 2*xyz[1], 2*xyz[2] }};
		for (int i = 0; i < numPositionsToTest; i++) {
			double x = Math.random();
			double y = Math.random();
			double z = (Math.random() > 0.5 ? 1 : -1) * Math.random();
			double[] pos = new double[] {x, y, z};
			
			double correction = x*x;
			double expected = Math.sqrt(x*x + correction*y*y + z*z 
					+ correction*1.5*x*1.5*x + correction*correction*1.5*y*1.5*y + correction*1.5*z*1.5*z
					+ 2*x*2*x + correction*2*y*2*y + 2*z*2*z);
			assertEquals(expected, csys.magnitude(pos, TensorIndexType.CONTRAVARIANT, tensorfield), 1e-7);
			
			correction = 1/(x*x);
			expected = Math.sqrt(x*x + correction*y*y + z*z 
					+ correction*1.5*x*1.5*x + correction*correction*1.5*y*1.5*y + correction*1.5*z*1.5*z
					+ 2*x*2*x + correction*2*y*2*y + 2*z*2*z);
			assertEquals(expected, csys.magnitude(pos, TensorIndexType.COVARIANT, tensorfield), 1e-7);
		}
	}

	@DisplayName("Testing gradient")
	@Test
	void testGradient() {
		CoordinateSystem csys = CoordinateSystem.from("cart");
		int numPositionsToTest = 15;
		
		Random rnd = new Random(125L);
		
		// test scalar field
		Function<double[], Double> scalarfield = xyz -> -xyz[1];
		for (int i = 0; i < numPositionsToTest; i++) {
			double y = rnd.nextGaussian() + 1;
			double[] pos = {rnd.nextDouble(), y, rnd.nextDouble()};
			assertArrayEquals(new double[] {0, -1, 0}, Types.unbox(csys.grad(pos, TensorIndexType.COVARIANT, scalarfield)), 1e-6);
			assertArrayEquals(new double[] {0, -1, 0}, Types.unbox(csys.grad(pos, TensorIndexType.CONTRAVARIANT, scalarfield)), 1e-6);
		}
		
		csys = CoordinateSystem.from("cyl");
		scalarfield = position -> 2*position[0]+3*position[1]+4*position[2];
		for (int i = 0; i < numPositionsToTest; i++) {
			double x = rnd.nextDouble() + 1;
			double[] pos = {x, rnd.nextDouble(), rnd.nextDouble()};
			assertArrayEquals(new double[] {2,3/x,4}, Types.unbox(csys.grad(pos, TensorIndexType.COVARIANT, scalarfield)), 5e-5);
			assertArrayEquals(new double[] {2,3/x,4}, Types.unbox(csys.grad(pos, TensorIndexType.CONTRAVARIANT, scalarfield)), 5e-5);
		}

		// test vector field
		csys = CoordinateSystem.from("cart");
		Function<double[], double[]> vectorfield = position -> new double[] { position[0],position[1],position[0]*position[2] };
		double[] pos = {0.5,1,2};
		double[][] expected = {{1,0,2},{0,1,0},{0,0,0.5}};
		double[][] actual = csys.grad(pos, TensorIndexType.COVARIANT, vectorfield);
		for (int i = 0; i < actual.length; i++) {
			assertArrayEquals(expected[i], actual[i], 5e-5);
		}

		csys = CoordinateSystem.from("cyl");
		vectorfield = position -> new double[] { position[0],position[1],position[0]*position[2] };
		pos = new double[] {0.5,1,2.5};
		expected = new double[][] {{1,0,2.5},{-2,3,0},{0,0,0.5}};
		actual = csys.grad(pos, TensorIndexType.COVARIANT, vectorfield);
		for (int i = 0; i < actual.length; i++) {
			assertArrayEquals(expected[i], actual[i], 5e-5);
		}
	}

	@DisplayName("Testing curl")
	@Test
	void testCurl() {
		// test vector field
		CoordinateSystem csys = CoordinateSystem.from("cart");
		double[] pos = {0.5,1,2};
		Function<double[], double[]> vectorfield = position -> new double[] { position[1],-position[0],0 };
		assertArrayEquals(new double[] {0,0,-2}, csys.curl(pos, TensorIndexType.CONTRAVARIANT, vectorfield), 5e-5);

		vectorfield = position -> new double[] { 0,-position[0]*position[0],0 };
		assertArrayEquals(new double[] {0,0,-1}, csys.curl(pos, TensorIndexType.CONTRAVARIANT, vectorfield), 5e-5);

		csys = CoordinateSystem.from("cyl");
		vectorfield = position -> new double[] { position[1],-position[0],0 };
		assertArrayEquals(new double[] {0,0,-3}, csys.curl(pos, TensorIndexType.CONTRAVARIANT, vectorfield), 5e-5);
	}

	@DisplayName("Testing units on axes")
	@Test
	void testAxesUnits() {
		CoordinateSystem csys = CoordinateSystem.from("cart");
		CoordinateSystem csys2 = CoordinateSystem.from("cart", new Axis(0, "km"), new Axis(1, "mm"), new Axis(2, "um"));
		
		Random rnd = new Random(125L);
		int numPositionsToTest = 15;
		
		for (int i = 0; i < numPositionsToTest; i++) {
			double x = rnd.nextGaussian();
			double y = rnd.nextGaussian();
			double z = rnd.nextGaussian();
			double[] pos = {x, y, z};
			assertArrayEquals(new double[] {1e-3*x, 1e3*y, 1e6*z}, CoordinateSystems.transform(pos, csys, csys2), 1e-8);
		}

		csys2 = CoordinateSystem.from("cyl", new Axis(2, "um"));
		
		for (int i = 0; i < numPositionsToTest; i++) {
			double x = rnd.nextGaussian();
			double y = rnd.nextGaussian();
			double z = rnd.nextGaussian();
			double[] pos = {x, y, z};
			assertEquals(1e6*z, CoordinateSystems.transform(pos, csys, csys2)[2], 1e-8);
		}
	}
}
