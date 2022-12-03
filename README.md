# Coordinates

[![](https://img.shields.io/github/issues/uhoefel/coordinates?style=flat-square)](https://github.com/uhoefel/coordinates/issues)
[![](https://img.shields.io/github/stars/uhoefel/coordinates?style=flat-square)](https://github.com/uhoefel/coordinates/stargazers)
[![DOI](https://zenodo.org/badge/327729481.svg)](https://zenodo.org/badge/latestdoi/327729481)
[![Maven Central](https://img.shields.io/maven-central/v/eu.hoefel/coordinates.svg?label=Maven%20Central)](https://search.maven.org/search?q=g:%22eu.hoefel%22%20AND%20a:%22coordinates%22)
[![](https://img.shields.io/github/license/uhoefel/coordinates?style=flat-square)](https://choosealicense.com/licenses/mit/)

Coordinates is a [Java](https://openjdk.java.net/) library designed to handle coordinate systems (including curvilinear, non-orthogonal ones).

Some of the supported features include:
- coordinate systems take [units](https://github.com/uhoefel/units) into account. For example
  ```java
  CoordinateSystem csys  = new CartesianCoordinates(3); // use constructor directly. 3 for 3D
  CoordinateSystem csys2 = CoordinateSystem.from("cart", 3, new Axis(0, "km")); // string-based construction
  double[] pos = {1000, 2, 3};
  CoordinateSystems.transform(pos, csys, csys2); // yields {1,2,3} as "m" is the default
  ```
- [numerical implementations](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/CoordinateSystem.java) for a lot of functions, e.g. the metric tensor, gradients, divergences etc.
- for the coordinate systems shipped symbolic implementations are provided where possible
- a [python script](https://github.com/uhoefel/coordinates/blob/master/src/main/python/coordinate_system_implementation_generator.py) to easily create full symbolic implementations, including examples (e.g. [here](https://github.com/uhoefel/coordinates/blob/master/src/main/python/six_sphere_coordinates.py)) from the shipped implementations generated via the script
- full implementations:
    - [Cartesian coordinates](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/CartesianCoordinates.java) (*n*D)
    - [Polar coordinates](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/PolarCoordinates.java) (2D)
    - [Bipolar coordinates](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/BipolarCoordinates.java) (2D)
    - [Bispherical coordinates](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/BisphericalCoordinates.java) (3D)
    - [Cylindrical coordinates](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/CylindricalCoordinates.java) (3D)
    - [Oblate spheroidal coordinates](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/OblateSpheroidalCoordinates.java) (3D)
    - [Prolate spheroidal coordinates](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/ProlateSpheroidalCoordinates.java) (3D)
    - [6-sphere coordinates](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/SixSphereCoordinates.java) (3D)
    - [Spherical coordinates](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/SphericalCoordinates.java) (*n*D)
    - [Toroidal coordinates](https://github.com/uhoefel/coordinates/blob/master/src/main/java/eu/hoefel/coordinates/ToroidalCoordinates.java) (3D)

Installation
============

The artifact can be found at maven central:
```xml
<dependency>
    <groupId>eu.hoefel</groupId>
    <artifactId>coordinates</artifactId>
    <version>1.1.3</version>
</dependency>
```

Requirements
============
Coordinates is designed to work with Java 17+.
