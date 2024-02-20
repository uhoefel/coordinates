package eu.hoefel.coordinates;

import java.lang.annotation.Documented;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * This annotation allows to decorate implementations of
 * {@link CoordinateSystem} to use the specified strings safely as their
 * representation in, e.g.,
 * {@link CoordinateSystems#transform(double[], String, String, java.util.Set)}.
 * Note that also anonymous inner classes can be annotated.
 * 
 * @author Udo Hoefel
 * @coordinates.apiNote The reason for using an annotation is the following: We
 *                      need the list of symbols to determine which coordinate
 *                      systems may be relevant for, e.g., a string based
 *                      coordinate transformation. However, at that point we may
 *                      not yet have an instance of the coordinate system
 *                      implementation, so we cannot simply call an instance
 *                      method to get the list of symbols. While we can try to
 *                      instantiate all potential classes and get the symbols
 *                      then, this has two drawbacks: First, performance: one
 *                      would have to use try-catch to control the program flow,
 *                      as some coordinate systems, for which the arguments were
 *                      not meant may throw an exception. This should be quite
 *                      slow compared to looking through some lists of strings.
 *                      Second, and that is quite important as well: if the user
 *                      gives wrong arguments, but the correct name of an
 *                      implementation, the method will not be able to know that
 *                      the arguments were wrong, as it cannot get the list of
 *                      symbols from the implementation (Which cannot be
 *                      instantiated as the arguments are illegal). This means
 *                      the user will get a way less helpful error message.
 *                      Static methods cannot be enforced via interfaces (and
 *                      would be more clumsy here in my opinion), so they offer
 *                      no great solution here.
 */
@Documented
@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.TYPE, ElementType.TYPE_USE})
public @interface CoordinateSystemSymbols {

    /**
     * The symbols that represent the coordinate system.
     * 
     * @return the representations of the coordinate system, e.g. "cartesian" and
     *         "cart" for {@link CartesianCoordinates}
     */
    public String[] value();
}
