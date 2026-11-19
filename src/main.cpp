#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/number_utils.h>
#include <manifold/manifold.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Vec3 {
  double x;
  double y;
  double z;
};

using Triangle = std::array<std::size_t, 3>;

struct Mesh {
  std::vector<Vec3> vertices;
  std::vector<Triangle> triangles;
};

Vec3 operator-(Vec3 a, Vec3 b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 cross(Vec3 a, Vec3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

double norm(Vec3 v) { return std::sqrt(dot(v, v)); }

const std::vector<Vec3>& inputPoints() {
  static const std::vector<Vec3> points = {
      {-0x1.4ee13237d80d4p-4, -0x1.54f3a9d7d37aep-2,
       -0x1.f000b913035b4p-2},
      {-0x1.4ecbe1c958b75p-4, 0x1.8a96df0e83f3fp-4,
       -0x1.f923c9ca07652p-2},
      {-0x1.4ebe87b465330p-4, -0x1.5b630182f3de3p-2,
       -0x1.ee4c1ca0ce091p-2},
      {-0x1.4eaed7eb21207p-4, -0x1.5bed7764f1714p-2,
       -0x1.f01508a9dae71p-2},
      {-0x1.4e64e42d56f47p-4, -0x1.817cea083cc71p-2,
       -0x1.e1122db65c166p-2},
      {-0x1.4dcc64b7e7724p-4, 0x1.e10bf8bad71ddp-2,
       -0x1.9651092358dacp-2},
      {-0x1.4ca036112f1d9p-4, 0x1.e3c2acc267246p-2,
       -0x1.986db319ad40fp-2},
      {-0x1.49d80d63bca57p-4, -0x1.e3c2acc267246p-2,
       -0x1.b601139fa3bd3p-2},
      {-0x1.2e83f7e23700fp-4, -0x1.5161394d4e846p-2,
       -0x1.f27bc780a9880p-2},
      {-0x1.2940057f18871p-4, -0x1.5b782f364f7c1p-2,
       -0x1.f033ce1022f61p-2},
      {-0x1.0037fc2f3c812p-4, 0x1.5027560c5211bp-4,
       -0x1.fcdcd821a6618p-2},
      {-0x1.001ecde4cbf5ep-4, -0x1.5557ed4bbbb1dp-2,
       -0x1.f28d776dde517p-2},
      {-0x1.fcccae5ccc9ecp-5, 0x1.57e4f928d31e0p-4,
       -0x1.fce92572ae4fap-2},
      {-0x1.f5c5a180c3054p-5, -0x1.54bb5a13585f5p-2,
       -0x1.f2c911cd28d6ap-2},
      {-0x1.f1640ffde7936p-5, -0x1.549579bcb9dc5p-2,
       -0x1.f2d71ae97d164p-2},
      {-0x1.a78f0c6470bdbp-5, 0x1.db88ef574eff1p-2, 0x1.0000000000000p-1},
      {-0x1.a0cfe78e2f7c4p-5, 0x1.58263fc8244f9p-4,
       -0x1.fb9b75f231108p-2},
      {-0x1.976ef1b99b218p-5, -0x1.e1a6d9775d6a8p-2,
       -0x1.b9a4f5f92f7cfp-2},
      {-0x1.93f14a45d6362p-5, -0x1.8241de2e7a09bp-2,
       -0x1.e0bf05a3a2da6p-2},
      {0x1.3e3c5713660ddp-8, 0x1.81250fe42a059p-2,
       -0x1.b2d25bab55aa3p-2},
      {0x1.50788e898380ap-8, -0x1.4298ae086b33ap-2,
       -0x1.f980106770b01p-2},
      {0x1.5ff715462e815p-8, 0x1.1b74184752ba9p-4, -0x1.0000000000000p-1},
      {0x1.0e5a88700a148p-4, -0x1.b9ba3dd813082p-2, 0x1.dbcc6d295e7edp-2},
      {0x1.0f38dabae5c57p-4, -0x1.b5bebe5527465p-2,
       -0x1.c1f8e2b7e03eap-2},
      {0x1.4d5692e9101bep-4, 0x1.58b2b709e1ff3p-2, 0x1.f0251d0de687ap-2},
      {0x1.4d65e7e63f84fp-4, 0x1.537ea32dcc5f7p-2,
       -0x1.ad33626b40a29p-2},
      {0x1.4ecc8cfb8ef4ep-4, -0x1.98fb1ccbdb907p-2,
       -0x1.ccba48afbb19fp-2},
      {0x1.4ee13237d80d4p-4, 0x1.4c33de6a0fa11p-2,
       -0x1.a04a2b86d5249p-2}};
  return points;
}

Mesh manifoldHull(const std::vector<Vec3>& points) {
  std::vector<manifold::vec3> input;
  input.reserve(points.size());
  for (Vec3 point : points) input.push_back({point.x, point.y, point.z});

  const manifold::Manifold hull =
      manifold::Manifold::Hull(input).Simplify(0.0);
  if (hull.Status() != manifold::Manifold::Error::NoError)
    throw std::runtime_error("Manifold3D returned error status " +
                             std::to_string(static_cast<int>(hull.Status())));
  if (hull.IsEmpty() || !(hull.Volume() > 0.0))
    throw std::runtime_error("Manifold3D returned an empty or zero-volume hull");

  const manifold::MeshGL64 source = hull.GetMeshGL64();
  Mesh result;
  result.vertices.reserve(source.NumVert());
  for (std::size_t vertex = 0; vertex < source.NumVert(); ++vertex) {
    const std::size_t offset = vertex * source.numProp;
    result.vertices.push_back({source.vertProperties[offset],
                               source.vertProperties[offset + 1],
                               source.vertProperties[offset + 2]});
  }
  result.triangles.reserve(source.NumTri());
  for (std::size_t face = 0; face < source.triVerts.size(); face += 3) {
    result.triangles.push_back(
        {static_cast<std::size_t>(source.triVerts[face]),
         static_cast<std::size_t>(source.triVerts[face + 1]),
         static_cast<std::size_t>(source.triVerts[face + 2])});
  }
  return result;
}

Mesh cgalHull(const std::vector<Vec3>& points) {
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Point = Kernel::Point_3;
  using SurfaceMesh = CGAL::Surface_mesh<Point>;

  std::vector<Point> input;
  input.reserve(points.size());
  for (Vec3 point : points) input.emplace_back(point.x, point.y, point.z);

  SurfaceMesh source;
  CGAL::convex_hull_3(input.begin(), input.end(), source);
  if (source.is_empty()) throw std::runtime_error("CGAL returned an empty hull");

  Mesh result;
  std::vector<std::size_t> vertexMap(source.number_of_vertices());
  result.vertices.reserve(source.number_of_vertices());
  for (SurfaceMesh::Vertex_index vertex : source.vertices()) {
    const Point& point = source.point(vertex);
    vertexMap[vertex.idx()] = result.vertices.size();
    result.vertices.push_back({CGAL::to_double(point.x()),
                               CGAL::to_double(point.y()),
                               CGAL::to_double(point.z())});
  }

  result.triangles.reserve(source.number_of_faces());
  for (SurfaceMesh::Face_index face : source.faces()) {
    Triangle triangle{};
    std::size_t corner = 0;
    for (SurfaceMesh::Vertex_index vertex :
         CGAL::vertices_around_face(source.halfedge(face), source)) {
      if (corner == triangle.size())
        throw std::runtime_error("CGAL returned a non-triangular face");
      triangle[corner++] = vertexMap[vertex.idx()];
    }
    if (corner != triangle.size())
      throw std::runtime_error("CGAL returned a degenerate face");
    result.triangles.push_back(triangle);
  }
  return result;
}

double signedVolume(const Mesh& mesh) {
  double volume6 = 0.0;
  for (const Triangle& triangle : mesh.triangles) {
    const Vec3 a = mesh.vertices[triangle[0]];
    const Vec3 b = mesh.vertices[triangle[1]];
    const Vec3 c = mesh.vertices[triangle[2]];
    volume6 += dot(a, cross(b, c));
  }
  return volume6 / 6.0;
}

double worstOutsideDistance(const std::vector<Vec3>& points,
                            const Mesh& mesh) {
  double worst = 0.0;
  for (const Triangle& triangle : mesh.triangles) {
    const Vec3 a = mesh.vertices[triangle[0]];
    const Vec3 b = mesh.vertices[triangle[1]];
    const Vec3 c = mesh.vertices[triangle[2]];
    Vec3 normal = cross(b - a, c - a);
    const double length = norm(normal);
    if (!(length > 0.0))
      throw std::runtime_error("backend returned a degenerate triangle");
    normal = {normal.x / length, normal.y / length, normal.z / length};
    for (Vec3 point : points) worst = std::max(worst, dot(point - a, normal));
  }
  return worst;
}

bool report(const char* name, const std::vector<Vec3>& points, const Mesh& mesh,
            double tolerance) {
  const double outside = worstOutsideDistance(points, mesh);
  const double volume = signedVolume(mesh);
  const bool pass = volume > 0.0 && outside <= tolerance;
  std::cout << std::left << std::setw(12) << name << ": "
            << (pass ? "PASS" : "FAIL") << " (vertices=" << mesh.vertices.size()
            << ", faces=" << mesh.triangles.size() << ", volume=" << volume
            << ", worst outside=" << outside << ", tolerance=" << tolerance
            << ")\n";
  return pass;
}

}  // namespace

int main() {
  try {
    std::cout << std::setprecision(12);
    const std::vector<Vec3>& points = inputPoints();
    constexpr double tolerance = 1e-6;

    const bool manifoldPass =
        report("Manifold3D", points, manifoldHull(points), tolerance);
    const bool cgalPass = report("CGAL", points, cgalHull(points), tolerance);

    if (!manifoldPass && cgalPass)
      std::cerr << "Reproduced: Manifold3D violates containment while CGAL "
                   "passes the same check.\n";
    return manifoldPass && cgalPass ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 2;
  }
}
