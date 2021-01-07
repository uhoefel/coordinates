package eu.hoefel.coordinates;

import eu.hoefel.utils.Maths;

/**
 * Helper class for coordinate transformation methods. See for example the
 * <a href="https://doi.org/10.1007/978-94-009-8352-6">appendix</a> in the book
 * "Low Reynolds number hydrodynamics" by Happel and Brenner.
 * 
 * @author Udo Hoefel
 */
final class Transformations {

	/** Class is only a utility class! */
	private Transformations() {
	    throw new IllegalStateException("This is a pure utility class!");
	}

	/**
	 * Transforms <i>n</i>-D spherical coordinates to Cartesian coordinates.
	 * 
	 * @param coords the spherical coordinates
	 *               (r*cos(theta_1),r*sin(theta_1)*cos(theta_2),r*sin(theta_1)*sin(theta_2)*cos(theta_3),…,r*sin(theta_1)*…*sin(theta_(n-1)))
	 * @return the Cartesian coordinates (x,y,z,…)
	 */
	public static final double[] sphericalToCartesian(double[] coords) {
		double[] ret = new double[coords.length];
		for (int i = 0; i < ret.length; i++) {
			ret[i] = coords[0];
			for (int j = 0; j < i; j++) {
				ret[i] = ret[i] * Math.sin(coords[j + 1]);
			}
		}

		for (int i = 0; i < ret.length; i++) {
			if (i != ret.length - 1) {
				ret[i] = ret[i] * Math.cos(coords[i + 1]);
			}
		}

		return ret;
	}

	/**
	 * Transforms <i>n</i>-D Cartesian coordinates to spherical coordinates.
	 * 
	 * @param coords the Cartesian coordinates (x,y,z,…). Should at least be 3D.
	 * @return the spherical coordinates
	 *         (r*cos(theta_1),r*sin(theta_1)*cos(theta_2),r*sin(theta_1)*sin(theta_2)*cos(theta_3),…,r*sin(theta_1)*…*sin(theta_(n-1)))
	 */
	public static final double[] cartesianToSpherical(double[] coords) {
		double[] ret = new double[coords.length];
		double r = 0;
		for (int j = 0; j < coords.length; j++) {
			r += Math.pow(coords[j], 2);
		}
		ret[0] = Math.sqrt(r);

		for (int i = 0; i < coords.length - 1; i++) {
			double radicand = 0;
			for (int j = i + 1; j < coords.length; j++) {
				radicand += Math.pow(coords[j], 2);
			}
			
			ret[i + 1] = Maths.acot(coords[i] / Math.sqrt(radicand));
		}

		ret[coords.length - 1] = Maths.acot(coords[coords.length - 2] / coords[coords.length - 1]);

		return ret;
	}

	/**
	 * Transforms 3D 6-sphere coordinates to Cartesian coordinates.
	 * 
	 * @param coords the 6-sphere coordinates (u,v,w)
	 * @return the Cartesian coordinates (x,y,z)
	 */
	public static final double[] sixSphereToCartesian(double[] coords) {
		double[] ret = new double[3];
		double denominator = Math.pow(coords[0], 2) + Math.pow(coords[1], 2) + Math.pow(coords[2], 2);
		ret[0] = coords[0] / denominator;
		ret[1] = coords[1] / denominator;
		ret[2] = coords[2] / denominator;
		return ret;
	}

	/**
	 * Transforms 3D Cartesian coordinates to 6-sphere coordinates.
	 * 
	 * @param coords the Cartesian coordinates (x,y,z)
	 * @return the 6-sphere coordinates (u,v,w)
	 */
	public static final double[] cartesianToSixSphere(double[] coords) {
		double[] ret = new double[3];
		double denominator = Math.pow(coords[0], 2) + Math.pow(coords[1], 2) + Math.pow(coords[2], 2);
		ret[0] = coords[0] / denominator;
		ret[1] = coords[1] / denominator;
		ret[2] = coords[2] / denominator;
		return ret;
	}

	/**
	 * Transforms 3D oblate spheroidal coordinates to Cartesian coordinates.
	 * 
	 * @param coords the oblate spheroidal coordinates (mu,nu,phi)
	 * @param a      the distance of the focal points to the origin of the
	 *               generating 2D elliptic coordinates. One can imagine it like
	 *               this: Rotate the elliptical coordinate system about the
	 *               perpendicular bisector of the focal points (separated by 2a),
	 *               i.e. rotate around the y-axis. We take the focal points to lie
	 *               on the x-axis.
	 * @return the Cartesian coordinates (x,y,z)
	 */
	public static final double[] oblateSpheroidalToCartesian(double[] coords, double a) {
		double[] ret = new double[3];
		ret[0] = a * Math.cosh(coords[0]) * Math.cos(coords[1]) * Math.cos(coords[2]);
		ret[1] = a * Math.cosh(coords[0]) * Math.cos(coords[1]) * Math.sin(coords[2]);
		ret[2] = a * Math.sinh(coords[0]) * Math.sin(coords[1]);
		return ret;
	}

	/**
	 * Transforms 3D prolate spheroidal coordinates to Cartesian coordinates.
	 * 
	 * @param coords the prolate spheroidal coordinates (mu,nu,phi)
	 * @param a      the distance of the focal points to the origin of the
	 *               generating 2D elliptic coordinates. One can imagine it like
	 *               this: Rotate the elliptical coordinate system about the focal
	 *               axis of the ellipse, i.e., the symmetry axis on which the foci
	 *               are located (separated by 2a).
	 * @return the Cartesian coordinates (x,y,z)
	 */
	public static final double[] prolateSpheroidalToCartesian(double[] coords, double a) {
		double[] ret = new double[3];
		ret[0] = a * Math.sinh(coords[0]) * Math.sin(coords[1]) * Math.cos(coords[2]);
		ret[1] = a * Math.sinh(coords[0]) * Math.sin(coords[1]) * Math.sin(coords[2]);
		ret[2] = a * Math.cosh(coords[0]) * Math.cos(coords[1]);
		return ret;
	}

	/**
	 * Transforms 3D conical coordinates to Cartesian coordinates.
	 * 
	 * @param coords the conical coordinates (r,mu,nu)
	 * @param b      a scale factor
	 * @param c      a scale factor
	 * @return the Cartesian coordinates (x,y,z)
	 */
	public static final double[] conicalToCartesian(double[] coords, double b, double c) {
		double[] ret = new double[3];
		// check if input is sensible
		if (coords[2] >= c * c) throw new IllegalArgumentException("nu^2 has to be smaller than c^2!");
		if (c * c >= coords[1]) throw new IllegalArgumentException("c^2 has to be smaller than mu^2!");
		if (coords[1] >= b * b) throw new IllegalArgumentException("mu^2 has to be smaller than b^2!");

		ret[0] = coords[0] * coords[1] * coords[2] / (b * c);
		ret[1] = (coords[0] / b) * Math.sqrt((Math.pow(coords[1], 2) - b * b) * (Math.pow(coords[2], 2) - b * b) / (b * b - c * c));
		ret[2] = (coords[0] / c) * Math.sqrt((Math.pow(coords[1], 2) - c * c) * (Math.pow(coords[2], 2) - c * c) / (c * c - b * b));
		return ret;
	}

	/**
	 * Transforms 3D parabolic cylindrical coordinates to Cartesian coordinates.
	 * There is no inverse as the inverse transformation is not unique.
	 * 
	 * @param coords the parabolic cylindrical coordinates (sigma,tau,z)
	 * @return the Cartesian coordinates (x,y,z)
	 */
	public static final double[] parabolicCylindricalToCartesian(double[] coords) {
		double[] ret = new double[3];
		ret[0] = coords[0] * coords[1];
		ret[1] = 0.5 * (Math.pow(coords[1], 2) - Math.pow(coords[0], 2));
		ret[2] = coords[2];
		return ret;
	}

	/**
	 * Transforms 3D bispherical coordinates to Cartesian coordinates.
	 * 
	 * @param coords the bispherical coordinates (sigma,tau,theta)
	 * @param a      the position of the focal points of the generating 2D bipolar
	 *               coordinates. One can imagine it like this: Rotate the bipolar
	 *               coordinate system about the axis joining its two foci
	 *               (separated by 2a).
	 * @return the Cartesian coordinates (x,y,z)
	 */
	public static final double[] bisphericalToCartesian(double[] coords, double a) {
		double[] ret = new double[3];
		double factor = a / (Math.cosh(coords[1]) - Math.cos(coords[0]));
		ret[0] = factor * Math.sin(coords[0]) * Math.cos(coords[2]);
		ret[1] = factor * Math.sin(coords[0]) * Math.sin(coords[2]);
		ret[2] = factor * Math.sinh(coords[1]);
		return ret;
	}
}
