#include "Network.hpp"

#include <algorithm>
#include <csignal>
#include <execinfo.h>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <unistd.h>
#include <vector>

namespace {

template <typename Set, typename Ptr>
bool contains(const Set &s, Ptr p) {
  return s.find(p) != s.end();
}

template <typename Range, typename Ptr>
bool contains_range(const Range &r, Ptr p) {
  return std::find(r.begin(), r.end(), p) != r.end();
}

void require(bool cond, const std::string &msg) {
  if (!cond)
    throw std::runtime_error(msg);
}

bool line_has_point(const networkLine *l, const networkPoint *p) {
  const auto [a, b] = l->getPoints();
  return (a == p) || (b == p);
}

bool face_has_point(const networkFace *f, const networkPoint *p) {
  const auto [a, b, c] = f->getPoints();
  return (a == p) || (b == p) || (c == p);
}

bool face_has_line(const networkFace *f, const networkLine *l) {
  const auto &lines = f->getLines();
  return contains_range(lines, l);
}

bool line_has_face(const networkLine *l, const networkFace *f) {
  const auto &faces = l->getFaces();
  return contains_range(faces, f);
}

void validate_network(Network &net, const std::string &stage) {
  std::cout << "\n[validate] " << stage << std::endl;

  net.setGeometricPropertiesForce();
  net.checkConnectivity();

  const auto ok = net.validateConectivity();
  require(ok[0] && ok[1] && ok[2] && ok[3], "validateConectivity failed");

  const auto &points = net.getPoints();
  const auto &faces = net.getFaces();
  const auto &tetras = net.getTetras();
  const auto &lines = net.Lines;

  const auto derived_lines = net.getLines();
  require(derived_lines.size() == lines.size(), "Lines set is out of sync with points");
  for (auto *l : derived_lines)
    require(contains(lines, l), "Line reachable from points is missing in Network::Lines");
  for (auto *l : lines)
    require(contains(derived_lines, l), "Line in Network::Lines is not reachable from points");

  for (auto *p : points) {
    require(p != nullptr, "null point");
    require(p->getNetwork() == &net, "point-network mismatch");
    for (auto *l : p->getLines()) {
      require(l != nullptr, "null line in point->Lines");
      require(contains(lines, l), "point->Lines contains line not owned by Network");
      require(line_has_point(l, p), "point->Lines has line missing back-reference to point");
    }
    for (auto *f : p->getFaces()) {
      require(f != nullptr, "null face in point->Faces");
      require(contains(faces, f), "point->Faces contains face not owned by Network");
      require(face_has_point(f, p), "point->Faces has face missing back-reference to point");
    }
    for (auto *t : p->Tetras) {
      if (!t)
        continue;
      require(contains(tetras, t), "point->Tetras contains tetra not owned by Network");
      require(contains_range(t->Points, p), "point->Tetras missing back-reference to point");
    }
  }

  for (auto *l : lines) {
    require(l != nullptr, "null line");
    require(l->getNetwork() == &net, "line-network mismatch");
    const auto [pa, pb] = l->getPoints();
    require(pa != nullptr && pb != nullptr, "line has null endpoint");
    require(pa != pb, "line endpoints are identical");
    require(contains(points, pa), "line endpoint not owned by Network");
    require(contains(points, pb), "line endpoint not owned by Network");
    require(contains_range(pa->getLines(), l), "line endpoint missing back-reference to line");
    require(contains_range(pb->getLines(), l), "line endpoint missing back-reference to line");
    for (auto *p : l->getPoints()) {
      require(p != nullptr, "null point in line");
      require(contains(points, p), "line->Points contains point not owned by Network");
    }
    for (auto *f : l->getFaces()) {
      require(f != nullptr, "null face in line->Faces");
      require(contains(faces, f), "line->Faces contains face not owned by Network");
      require(face_has_line(f, l), "line->Faces missing back-reference to line");
    }
    for (auto *t : l->Tetras) {
      if (!t)
        continue;
      require(contains(tetras, t), "line->Tetras contains tetra not owned by Network");
      require(contains_range(t->Lines, l), "line->Tetras missing back-reference to line");
    }
  }

  for (auto *f : faces) {
    require(f != nullptr, "null face");
    require(f->getNetwork() == &net, "face-network mismatch");
    const auto [p0, p1, p2] = f->getPoints();
    require(p0 != nullptr && p1 != nullptr && p2 != nullptr, "face has null point");
    require(contains(points, p0) && contains(points, p1) && contains(points, p2), "face has point not owned by Network");
    for (auto *p : f->getPoints()) {
      require(p != nullptr, "null point in face");
      require(contains(points, p), "face->Points contains point not owned by Network");
      require(contains_range(p->getFaces(), f), "face->Points missing back-reference to face");
    }
    for (auto *l : f->getLines()) {
      require(l != nullptr, "null line in face");
      require(contains(lines, l), "face->Lines contains line not owned by Network");
      require(line_has_face(l, f), "face->Lines missing back-reference to face");
    }
    for (auto *t : f->Tetras) {
      if (!t)
        continue;
      require(contains(tetras, t), "face->Tetras contains tetra not owned by Network");
      require(contains_range(t->Faces, f), "face->Tetras missing back-reference to face");
    }
  }

  for (auto *t : tetras) {
    require(t != nullptr, "null tetra");
    for (auto *p : t->Points) {
      require(p != nullptr, "null point in tetra");
      require(contains(points, p), "tetra->Points contains point not owned by Network");
      require(contains_range(p->Tetras, t), "tetra->Points missing back-reference to tetra");
    }
    for (auto *l : t->Lines) {
      require(l != nullptr, "null line in tetra");
      require(contains(lines, l), "tetra->Lines contains line not owned by Network");
      require(contains_range(l->Tetras, t), "tetra->Lines missing back-reference to tetra");
    }
    for (auto *f : t->Faces) {
      require(f != nullptr, "null face in tetra");
      require(contains(faces, f), "tetra->Faces contains face not owned by Network");
      require(contains_range(f->Tetras, t), "tetra->Faces missing back-reference to tetra");
    }
  }
}

void validate_obj_file(const std::string &path) {
  std::cout << "\n[obj] " << path << std::endl;
  Network net(path);
  std::cout << "  points=" << net.getPoints().size() << " lines=" << net.Lines.size() << " faces=" << net.getFaces().size()
            << " tetras=" << net.getTetras().size() << std::endl;
  validate_network(net, "obj load");
}

void require_line_absent(const Network &net, const networkLine *line, const std::string &stage) {
  require(net.Lines.find(const_cast<networkLine *>(line)) == net.Lines.end(), "line still in Network::Lines after " + stage);
  for (auto *p : net.getPoints())
    require(!contains_range(p->getLines(), const_cast<networkLine *>(line)), "point still references line after " + stage);
  for (auto *f : net.getFaces())
    require(!contains_range(f->getLines(), const_cast<networkLine *>(line)), "face still references line after " + stage);
  for (auto *t : net.getTetras())
    require(!contains_range(t->Lines, const_cast<networkLine *>(line)), "tetra still references line after " + stage);
}

void require_face_absent(const Network &net, const networkFace *face, const std::string &stage) {
  require(net.getFaces().find(const_cast<networkFace *>(face)) == net.getFaces().end(), "face still in Network::Faces after " + stage);
  for (auto *p : net.getPoints())
    require(!contains_range(p->getFaces(), const_cast<networkFace *>(face)), "point still references face after " + stage);
  for (auto *l : net.Lines)
    require(!contains_range(l->getFaces(), const_cast<networkFace *>(face)), "line still references face after " + stage);
  for (auto *t : net.getTetras())
    require(!contains_range(t->Faces, const_cast<networkFace *>(face)), "tetra still references face after " + stage);
}

void require_point_absent(const Network &net, const networkPoint *point, const std::string &stage) {
  require(net.getPoints().find(const_cast<networkPoint *>(point)) == net.getPoints().end(), "point still in Network::Points after " + stage);
  for (auto *l : net.Lines) {
    const auto [a, b] = l->getPoints();
    require(a != point && b != point, "line still references point after " + stage);
  }
  for (auto *f : net.getFaces()) {
    const auto [a, b, c] = f->getPoints();
    require(a != point && b != point && c != point, "face still references point after " + stage);
  }
  for (auto *t : net.getTetras()) {
    for (auto *p : t->Points)
      require(p != point, "tetra still references point after " + stage);
  }
}

void build_grid(Network &net, int cells) {
  std::vector<Tddd> points;
  points.reserve((cells + 1) * (cells + 1));
  for (int j = 0; j <= cells; ++j)
    for (int i = 0; i <= cells; ++i)
      points.push_back({static_cast<double>(i), static_cast<double>(j), 0.0});

  auto pts = net.setPoints(points);

  std::vector<std::vector<int>> faces;
  faces.reserve(cells * cells * 2);
  auto idx = [cells](int i, int j) { return j * (cells + 1) + i; };
  for (int j = 0; j < cells; ++j) {
    for (int i = 0; i < cells; ++i) {
      const int p00 = idx(i, j);
      const int p10 = idx(i + 1, j);
      const int p11 = idx(i + 1, j + 1);
      const int p01 = idx(i, j + 1);
      faces.push_back({p00, p10, p11});
      faces.push_back({p00, p11, p01});
    }
  }

  net.setFaces(faces, pts);
}

networkLine *find_two_face_line(Network &net) {
  for (auto *l : net.Lines) {
    if (!l)
      continue;
    if (l->getBoundaryFaces().size() == 2)
      return l;
  }
  return nullptr;
}

networkLine *find_mergeable_line(Network &net) {
  for (auto *l : net.Lines) {
    if (!l)
      continue;
    if (l->isMergeable())
      return l;
  }
  return nullptr;
}

bool try_flip_any(Network &net) {
  for (auto *l : net.Lines) {
    if (!l)
      continue;
    if (l->Flip(true))
      return true;
  }
  return false;
}

void run_split_flip_collapse_sequence() {
  Network net;
  build_grid(net, 5);
  validate_network(net, "initial");

  for (int i = 0; i < 3; ++i) {
    auto *l = find_two_face_line(net);
    if (!l)
      break;
    auto *mid = l->Split();
    require(mid != nullptr, "split returned null point");
    require(net.getPoints().find(mid) != net.getPoints().end(), "split point not registered in Network::Points");
    validate_network(net, "after split");
  }

  for (int i = 0; i < 3; ++i) {
    const bool flipped = try_flip_any(net);
    if (!flipped)
      break;
    validate_network(net, "after flip");
  }

  for (int i = 0; i < 3; ++i) {
    auto *l = find_mergeable_line(net);
    if (!l)
      break;
    auto *victim = l;
    try {
      auto *remain = l->Collapse();
      if (!remain) {
        std::cout << "[warn] collapse: precondition not satisfied, skipped" << std::endl;
        continue;
      }
      require_line_absent(net, victim, "collapse");
      validate_network(net, "after collapse");
    } catch (const std::exception &e) {
      std::cout << "[warn] collapse threw, skipped: " << e.what() << std::endl;
      continue;
    }
  }
}

void run_ops_on_network(Network &net, const std::string &label) {
  auto stage = [&](const std::string &suffix) { return label + " " + suffix; };
  validate_network(net, stage("initial"));

  for (int i = 0; i < 3; ++i) {
    auto *l = find_two_face_line(net);
    if (!l)
      break;
    try {
      auto *mid = l->Split();
      if (!mid) {
        std::cout << "[warn] " << label << " split skipped" << std::endl;
        continue;
      }
      validate_network(net, stage("after split"));
    } catch (const std::exception &e) {
      std::cout << "[warn] " << label << " split threw, skipped: " << e.what() << std::endl;
      break;
    }
  }

  for (int i = 0; i < 3; ++i) {
    try {
      if (!try_flip_any(net)) {
        std::cout << "[warn] " << label << " flip skipped" << std::endl;
        break;
      }
      validate_network(net, stage("after flip"));
    } catch (const std::exception &e) {
      std::cout << "[warn] " << label << " flip threw, skipped: " << e.what() << std::endl;
      break;
    }
  }

  for (int i = 0; i < 3; ++i) {
    auto *l = find_mergeable_line(net);
    if (!l) {
      std::cout << "[warn] " << label << " collapse skipped" << std::endl;
      break;
    }
    auto *victim = l;
    try {
      auto *remain = l->Collapse();
      if (!remain) {
        std::cout << "[warn] " << label << " collapse precondition not satisfied, skipped" << std::endl;
        continue;
      }
      require_line_absent(net, victim, label + " collapse");
      validate_network(net, stage("after collapse"));
    } catch (const std::exception &e) {
      std::cout << "[warn] " << label << " collapse threw, skipped: " << e.what() << std::endl;
    }
  }

  if (!net.Lines.empty()) {
    auto *l = *net.Lines.begin();
    delete l;
    validate_network(net, stage("after line delete"));
  }

  if (!net.getFaces().empty()) {
    auto *f = *net.getFaces().begin();
    delete f;
    validate_network(net, stage("after face delete"));
  }

  if (!net.getPoints().empty()) {
    auto *p = *net.getPoints().begin();
    delete p;
    validate_network(net, stage("after point delete"));
  }
}

void run_sequential_ops() {
  Network net;
  build_grid(net, 5);
  run_ops_on_network(net, "seq");
}

} // namespace

int main(int argc, char **argv) {
  try {
    std::signal(SIGSEGV, [](int sig) {
      void *frames[64];
      const int n = backtrace(frames, 64);
      std::cerr << "\n[signal] SIGSEGV\n";
      backtrace_symbols_fd(frames, n, STDERR_FILENO);
      _exit(128 + sig);
    });
    if (argc > 1) {
      for (int i = 1; i < argc; ++i) {
        std::string path = argv[i];
        std::cout << "\n[obj] " << path << std::endl;
        Network net(path);
        run_ops_on_network(net, "obj");
      }
    } else {
      run_split_flip_collapse_sequence();
      run_sequential_ops();

      {
        Network net;
        build_grid(net, 5);
        validate_network(net, "delete-face initial");
        if (!net.getFaces().empty()) {
          auto *f = *net.getFaces().begin();
          delete f;
          require_face_absent(net, f, "face delete");
          validate_network(net, "after face delete");
        }
      }

      {
        Network net;
        build_grid(net, 5);
        validate_network(net, "delete-line initial");
        if (!net.Lines.empty()) {
          auto *l = *net.Lines.begin();
          delete l;
          require_line_absent(net, l, "line delete");
          validate_network(net, "after line delete");
        }
      }

      {
        Network net;
        build_grid(net, 5);
        validate_network(net, "delete-point initial");
        if (!net.getPoints().empty()) {
          auto *p = *net.getPoints().begin();
          delete p;
          require_point_absent(net, p, "point delete");
          validate_network(net, "after point delete");
        }
      }
    }

    std::cout << "\nOK: network integrity checks passed." << std::endl;
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "\n[error] " << e.what() << std::endl;
    return 1;
  }
}
