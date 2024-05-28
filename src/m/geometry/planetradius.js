function planetradius(planet) {//{{{
	/**
	 * PLANETRADIUS - return planet radius according to planetary body name
	 *
	 * Usage:
	 *     radius = planetradius(planet);
	 *
	 * Examples:
	 * earthradius = planetradius('earth');
	 */

	let radius = 0;
	if (planet === 'earth') {
		radius = 6.371012e6;
	} else if (planet === 'europa') {
		radius = 1.5008e6;
	} else {
		error('planet type ' + planet + ' not supported yet!');
	}
} //}}}
