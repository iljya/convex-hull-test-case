#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/number_utils.h>
#include <manifold/manifold.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
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

Vec3 operator+(Vec3 a, Vec3 b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator*(Vec3 v, double scale) {
  return {v.x * scale, v.y * scale, v.z * scale};
}

Vec3 normalized(Vec3 v) {
  const double length = norm(v);
  if (!(length > 0.0)) throw std::runtime_error("cannot normalize zero vector");
  return v * (1.0 / length);
}

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

Vec3 meshCenter(const Mesh& mesh) {
  Vec3 center{0.0, 0.0, 0.0};
  for (Vec3 vertex : mesh.vertices) center = center + vertex;
  return center * (1.0 / static_cast<double>(mesh.vertices.size()));
}

Vec3 outwardNormal(const Mesh& mesh, const Triangle& triangle) {
  const Vec3 center = meshCenter(mesh);
  const Vec3 a = mesh.vertices[triangle[0]];
  const Vec3 b = mesh.vertices[triangle[1]];
  const Vec3 c = mesh.vertices[triangle[2]];
  Vec3 normal = normalized(cross(b - a, c - a));
  if (dot(center - a, normal) > 0.0) normal = normal * -1.0;
  return normal;
}

struct Violation {
  double distance;
  Vec3 normal;
  Triangle triangle;
};

Violation worstVisualViolation(const std::vector<Vec3>& points,
                               const Mesh& mesh) {
  Violation worst{0.0, {0.0, 0.0, 1.0}, {0, 0, 0}};
  for (const Triangle& triangle : mesh.triangles) {
    const Vec3 a = mesh.vertices[triangle[0]];
    const Vec3 normal = outwardNormal(mesh, triangle);
    for (Vec3 point : points) {
      const double distance = dot(point - a, normal);
      if (distance > worst.distance)
        worst = {distance, normal, triangle};
    }
  }
  return worst;
}

struct ProjectedPoint {
  double x;
  double y;
  double depth;
};

struct Projection {
  Vec3 right;
  Vec3 up;
  Vec3 forward;
  double minX;
  double maxX;
  double minY;
  double maxY;
};

ProjectedPoint project(Vec3 point, const Projection& projection) {
  return {dot(point, projection.right), -dot(point, projection.up),
          dot(point, projection.forward)};
}

Projection makeProjection(const std::vector<Vec3>& points,
                          const Mesh& manifold) {
  const Violation violation = worstVisualViolation(points, manifold);
  const Vec3 a = manifold.vertices[violation.triangle[0]];
  const Vec3 b = manifold.vertices[violation.triangle[1]];
  const Vec3 c = manifold.vertices[violation.triangle[2]];
  const std::array<Vec3, 3> edges = {b - a, c - b, a - c};
  const Vec3 longestEdge = *std::max_element(
      edges.begin(), edges.end(),
      [](Vec3 lhs, Vec3 rhs) { return norm(lhs) < norm(rhs); });
  const Vec3 right = normalized(longestEdge);
  const Vec3 up = violation.normal;
  const Vec3 forward = normalized(cross(right, up));
  Projection projection{right,
                        up,
                        forward,
                        std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()};
  for (Vec3 point : points) {
    const ProjectedPoint p = project(point, projection);
    projection.minX = std::min(projection.minX, p.x);
    projection.maxX = std::max(projection.maxX, p.x);
    projection.minY = std::min(projection.minY, p.y);
    projection.maxY = std::max(projection.maxY, p.y);
  }
  return projection;
}

void writeHullPlot(const std::filesystem::path& path, const char* title,
                   const std::vector<Vec3>& points, const Mesh& mesh,
                   const Projection& projection) {
  constexpr double width = 900.0;
  constexpr double height = 680.0;
  constexpr double left = 70.0;
  constexpr double right = 50.0;
  constexpr double top = 78.0;
  constexpr double bottom = 72.0;
  const double plotWidth = width - left - right;
  const double plotHeight = height - top - bottom;
  const double spanX = projection.maxX - projection.minX;
  const double spanY = projection.maxY - projection.minY;
  const double scale = 0.9 * std::min(plotWidth / spanX, plotHeight / spanY);
  const double centerX = (projection.minX + projection.maxX) * 0.5;
  const double centerY = (projection.minY + projection.maxY) * 0.5;

  const auto screen = [&](Vec3 point) {
    const ProjectedPoint p = project(point, projection);
    return ProjectedPoint{width * 0.5 + (p.x - centerX) * scale,
                          top + plotHeight * 0.5 + (p.y - centerY) * scale,
                          p.depth};
  };

  std::vector<std::size_t> faceOrder(mesh.triangles.size());
  std::iota(faceOrder.begin(), faceOrder.end(), 0);
  std::sort(faceOrder.begin(), faceOrder.end(), [&](std::size_t a,
                                                    std::size_t b) {
    const Triangle& ta = mesh.triangles[a];
    const Triangle& tb = mesh.triangles[b];
    const double da = (screen(mesh.vertices[ta[0]]).depth +
                       screen(mesh.vertices[ta[1]]).depth +
                       screen(mesh.vertices[ta[2]]).depth) /
                      3.0;
    const double db = (screen(mesh.vertices[tb[0]]).depth +
                       screen(mesh.vertices[tb[1]]).depth +
                       screen(mesh.vertices[tb[2]]).depth) /
                      3.0;
    return da < db;
  });

  std::ofstream output(path);
  if (!output) throw std::runtime_error("cannot write " + path.string());
  output << std::fixed << std::setprecision(3);
  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" "
            "height=\"680\" viewBox=\"0 0 900 680\">\n"
         << "<rect width=\"900\" height=\"680\" fill=\"white\"/>\n"
         << "<text x=\"70\" y=\"34\" font-family=\"system-ui, sans-serif\" "
            "font-size=\"21\" font-weight=\"600\" fill=\"#111827\">"
         << title << "</text>\n"
         << "<text x=\"70\" y=\"56\" font-family=\"system-ui, sans-serif\" "
            "font-size=\"13\" fill=\"#6b7280\">"
         << points.size() << " input points · " << mesh.vertices.size()
         << " hull vertices</text>\n";

  for (std::size_t face : faceOrder) {
    const Triangle& triangle = mesh.triangles[face];
    const ProjectedPoint a = screen(mesh.vertices[triangle[0]]);
    const ProjectedPoint b = screen(mesh.vertices[triangle[1]]);
    const ProjectedPoint c = screen(mesh.vertices[triangle[2]]);
    output << "<path d=\"M " << a.x << ' ' << a.y << " L " << b.x << ' '
           << b.y << " L " << c.x << ' ' << c.y
           << " Z\" fill=\"#3b82f6\" fill-opacity=\"0.13\" "
              "stroke=\"#2563eb\" stroke-opacity=\"0.62\" "
              "stroke-width=\"1.2\" stroke-linejoin=\"round\"/>\n";
  }

  for (Vec3 point : points) {
    const ProjectedPoint p = screen(point);
    output << "<circle cx=\"" << p.x << "\" cy=\"" << p.y
           << "\" r=\"3.8\" fill=\"#111827\" stroke=\"white\" "
              "stroke-width=\"1.2\"/>\n";
  }

  const double axisX = 84.0;
  const double axisY = 606.0;
  const auto axis = [&](Vec3 direction, const char* label) {
    const double dx = dot(direction, projection.right);
    const double dy = -dot(direction, projection.up);
    const double length = std::sqrt(dx * dx + dy * dy);
    const double endX = axisX + dx / length * 47.0;
    const double endY = axisY + dy / length * 47.0;
    output << "<line x1=\"" << axisX << "\" y1=\"" << axisY
           << "\" x2=\"" << endX << "\" y2=\"" << endY
           << "\" stroke=\"#4b5563\" stroke-width=\"1.5\"/>\n"
           << "<circle cx=\"" << endX << "\" cy=\"" << endY
           << "\" r=\"2\" fill=\"#4b5563\"/>\n"
           << "<text x=\"" << endX + (dx >= 0 ? 7 : -12) << "\" y=\""
           << endY + (dy >= 0 ? 14 : -5)
           << "\" font-family=\"system-ui, sans-serif\" font-size=\"13\" "
              "fill=\"#374151\">"
           << label << "</text>\n";
  };
  axis({1.0, 0.0, 0.0}, "x");
  axis({0.0, 1.0, 0.0}, "y");
  axis({0.0, 0.0, 1.0}, "z");

  output << "<g font-family=\"system-ui, sans-serif\" font-size=\"13\" "
            "fill=\"#4b5563\">\n"
         << "<line x1=\"642\" y1=\"626\" x2=\"672\" y2=\"626\" "
            "stroke=\"#2563eb\" stroke-width=\"2\"/>\n"
         << "<text x=\"680\" y=\"631\">Hull surface</text>\n"
         << "<circle cx=\"642\" cy=\"650\" r=\"3.6\" fill=\"#111827\"/>\n"
         << "<text x=\"652\" y=\"655\">Input point</text>\n"
         << "</g>\n</svg>\n";
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

int main(int argc, char* argv[]) {
  try {
    std::filesystem::path plotDirectory;
    if (argc == 3 && std::string(argv[1]) == "--plot-dir") {
      plotDirectory = argv[2];
    } else if (argc != 1) {
      std::cerr << "Usage: " << argv[0] << " [--plot-dir DIRECTORY]\n";
      return 2;
    }

    std::cout << std::setprecision(12);
    const std::vector<Vec3>& points = inputPoints();
    constexpr double tolerance = 1e-6;

    const Mesh manifold = manifoldHull(points);
    const Mesh cgal = cgalHull(points);

    const bool manifoldPass =
        report("Manifold3D", points, manifold, tolerance);
    const bool cgalPass = report("CGAL", points, cgal, tolerance);

    if (!plotDirectory.empty()) {
      std::filesystem::create_directories(plotDirectory);
      const Projection projection = makeProjection(points, manifold);
      writeHullPlot(plotDirectory / "manifold.svg", "Manifold3D hull", points,
                    manifold, projection);
      writeHullPlot(plotDirectory / "cgal.svg", "CGAL hull", points, cgal,
                    projection);
      std::cout << "Wrote plots to " << plotDirectory << '\n';
    }

    if (!manifoldPass && cgalPass)
      std::cerr << "Reproduced: Manifold3D violates containment while CGAL "
                   "passes the same check.\n";
    return manifoldPass && cgalPass ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 2;
  }
}
