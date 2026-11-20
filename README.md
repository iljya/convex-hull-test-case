# Manifold3D near-coplanar convex-hull reproducer

Minimal reproduction of a Manifold3D QuickHull result that is closed and has
positive volume, but is not convex and does not contain every input point.
The same point cloud is also passed directly to CGAL's 3D convex hull as a
reference.

## Build and run

Requirements:

- CMake 3.22 or newer
- A C++17 compiler
- Boost 1.74 or newer headers

Clone with submodules, configure, build, and run the test:

```sh
git clone --recurse-submodules <repository-url>
cd manifold-hull-repro
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

If Boost is installed outside the compiler's normal include paths, add
`-DBoost_ROOT=/path/to/boost` to the configure command.

The regression test intentionally fails with the affected Manifold revision.
Its output should resemble:

```text
Manifold3D: FAIL (... worst outside distance about 1.1082)
CGAL:       PASS (... worst outside distance near zero)
```

The test succeeds once both backends return outward-oriented hulls containing
every input point within the scale-relative tolerance.

To generate matching SVG views of the input points and both hull results, run:

```sh
./build/manifold_hull_repro --plot-dir plots
```

The command writes `plots/manifold.svg` and `plots/cgal.svg`. Both images use
the same camera and scale, aligned to show Manifold's largest hull violation.

The submodules are pinned to Manifold3D v3.5.4 commit
`ce50d78021d64507f89e8c9fc2c2e51018117857` and CGAL v6.2 commit
`cac3e9d75e254928db0e38a3161564216cb01919`.

## License

The reproduction code is MIT licensed. The Manifold3D and CGAL submodules keep
their respective upstream licenses.
