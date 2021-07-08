package simulation.montecarlo;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Signifies that a method changes the lattice such that the magnetic moments
 * are incompatible with the magnetic fields they induce.
 * Means that either the lattice is not to be used for magnetic moment calculations
 * or that {@link Lattice#updateAllLocalFields()} should be used.
 * As a convention, the fields should always correspond to the magnetic moments which are then
 * adjusted using {@link Lattice#updateAllMagneticMoments(int, double, double)} to correspond
 * to the fields.
 */
@Retention(RetentionPolicy.SOURCE)
@Target({ElementType.METHOD,ElementType.CONSTRUCTOR})
public @interface CreatesInconsistency {
    public String value() default "";
}
