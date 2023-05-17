# vincenty-core

Calculate distances between two coordinates using vincenty formulae.

## Overview

This Rust module provides functions for calculating the geodesic distance
between two points specified by latitude/longitude on the surface of an
ellipsoid. The computation is based on the Vincenty's formulae, known for their
high accuracy.

Features:

- High Precision: The module is designed to return results with a precision of
  up to 6 decimal places.

- Custom Error Handling: The module uses Rust's Result type for error handling,
  ensuring that any issues during computation are handled gracefully and
  informatively.

- Convenient Input Types: The module provides two functions for distance
  calculation that accept different input types: one takes
  `geo_types::geometry::Coord` objects, and the other takes points as `f64`
  parameters directly.

## Usage

Add the module to your Rust project and use it as follows:

```rust
use geo_types::geometry::Coord; use your_module_name;

let coord1 = Coord { x: 40.712776, y: -74.005974 };
let coord2 = Coord { x: 34.052235, y: -118.243683 };

let distance = your_module_name::distance_from_coordinates(&coord1, &coord2);
println!("Distance between New York and Los Angeles: {} km", distance.unwrap());
```

## Testing

The module includes unit tests for its main functions. To run these tests, use
the following command:

```bash
$ cargo test
```

## Limitations

Please note that the algorithm performs a maximum of 200 iterations to find a
result within a specified threshold `1e-9`. If the algorithm fails to converge
within these iterations, an error is returned.

## Contributing

Feel free to fork this module and submit pull requests for any enhancements or
bug fixes. Make sure to include unit tests for any new functionality.

## License

This module is available under the MIT License.
