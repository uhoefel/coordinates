package eu.hoefel.coordinates.tensors;

import java.util.List;

/**
 * Describes how coordinate representations (like vectors, linear forms and
 * tensors) behave with respect to a base change in the underlying vector space.
 * <p>
 * <img src="doc-files/basischangebehavior.svg" alt="visualization of the
 * covariant and contravariant basis">
 * <p>
 * In this visualization, one can see the general curvilinear coordinate curves
 * (denoted <i>q</i><sup>1</sup>, <i>q</i><sup>2</sup> and <i>q</i><sup>3</sup>)
 * in black on the example of a 3D system. In the top right, the light gray
 * surfaces are the coordinate surfaces. The tangent basis is shown in yellow
 * (denoted <i>e</i><sub>1</sub>, <i>e</i><sub>2</sub> and <i>e</i><sub>3</sub>)
 * and the covector basis is shown in blue (denoted <i>e</i><sup>1</sup>,
 * <i>e</i><sup>2</sup> and <i>e</i><sup>3</sup>). The lower row shows an
 * example vector <i>v</i> in red, whose components can behave either
 * {@link #CONTRAVARIANT contravariantly} (using the tangent basis) or
 * {@link #COVARIANT covariantly} (using the covector basis). The contravariant
 * case is shown in the lower left. One can see that the vector addition of
 * <i>v</i><sup>1</sup><i>e</i><sub>1</sub>,
 * <i>v</i><sup>2</sup><i>e</i><sub>2</sub> and
 * <i>v</i><sup>3</sup><i>e</i><sub>3</sub> yields <i>v</i>. This is the most
 * common way to represent a vector. However, the same vector can also be
 * represented with covariant components, in which case the lower right
 * visualization may be helpful. In there, one can see that <i>v</i> in the
 * covector basis (i.e., the basis vectors are perpendicular to the coordinate
 * surfaces) corresponds to the vector addition of
 * <i>v</i><sub>1</sub><i>e</i><sup>1</sup>,
 * <i>v</i><sub>2</sub><i>e</i><sup>2</sup> and
 * <i>v</i><sub>3</sub><i>e</i><sup>3</sup>. In general, the tangent and
 * covector basis coincide only if the the coordinate system is orthogonal.
 * 
 * @author Udo Hoefel
 * 
 * @see <a href=
 *      "https://en.wikipedia.org/wiki/Covariance_and_contravariance_of_vectors">Covariance
 *      and contravariance of vectors</a>
 * @see <a href=
 *      "https://csm.mech.utah.edu/content/wp-content/uploads/2011/03/GoBagCurvilinear.pdf">Curvilinear
 *      Analysis in a Euclidean Space by Rebecca Moss Brannon</a>
 */
public enum TensorIndexType implements TensorTransformation {

    /**
     * The components of a vector are said to be covariant if they increase when the
     * basis vector length increases and vice versa. Covariancy is denoted by using
     * lower indices in Einstein notation.
     * 
     * @see <a href="https://www.youtube.com/watch?v=CliW7kSxxWU">Tensors Explained
     *      Intuitively: Covariant, Contravariant, Rank</a>
     */
    COVARIANT,

    /**
     * The components of a vector are said to be contravariant if they decrease when
     * the basis vector length increases and vice versa. Contravariancy is denoted
     * by using upper indices in Einstein notation. Most vectors have contravariant
     * components. For an example see also the picture in
     * {@link TensorIndexType}.
     * 
     * @see <a href="https://www.youtube.com/watch?v=CliW7kSxxWU">Tensors Explained
     *      Intuitively: Covariant, Contravariant, Rank</a>
     */
    CONTRAVARIANT;

    /**
     * Gets the other behavior with respect to a basis change, i.e. if the current
     * behavior with respect to a basis change is {@link #COVARIANT covariant}
     * return {@link #CONTRAVARIANT contravariant} and vice versa.
     * 
     * @return the other behavior with respect to a basis change 
     */
    public TensorIndexType flip() {
        return this == COVARIANT ? CONTRAVARIANT : COVARIANT;
    }

    @Override
    public List<TensorIndexType> indexTypes() {
        return List.of(this);
    }
}
