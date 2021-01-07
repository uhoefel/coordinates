package eu.hoefel.coordinates;

import java.util.List;
import java.util.stream.Stream;

import org.junit.jupiter.api.extension.ExtensionContext;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.ArgumentsProvider;

/**
 * Provides instances of the default coordinate systems.
 * 
 * @author Udo Hoefel
 */
class DefaultCoordinateSystemProvider implements ArgumentsProvider {

	@Override
	public Stream<? extends Arguments> provideArguments(ExtensionContext context) {
		return CoordinateSystems.DEFAULT_COORDINATE_SYSTEMS
				.stream()
				.map(DefaultCoordinateSystemProvider::instance)
				.map(Arguments::of);
	}

	/**
	 * Gets an instance of the specified coordinate system class.
	 * 
	 * @param clazz the class of the coordinate system
	 * @return a coordinate system instance of the specified class
	 */
	private static final CoordinateSystem instance(Class<? extends CoordinateSystem> clazz) {
		List<String> symbols = CoordinateSystems.symbolsFromClass(clazz);
		try {
			return CoordinateSystem.from(symbols.get(0), CoordinateSystems.DEFAULT_COORDINATE_SYSTEMS);
		} catch (Exception e) {
			// So we need one or two additional parameters, probably.
			// We will just provide some default values
			return CoordinateSystem.from(symbols.get(0), CoordinateSystems.DEFAULT_COORDINATE_SYSTEMS, "1 m", "1 m", 1.0);
		}
	}
}
