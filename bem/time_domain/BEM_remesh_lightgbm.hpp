#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace BEMMeshPipeline::RemeshAI {

struct RemeshScenarioFeatures {
  int op = -1;
  int reason = -1;
  int edge_category = -1;
  double len = -1.0;
  double len_ratio = -1.0;
  double theta = -1.0;
  double theta_ratio = -1.0;
  std::string scenario;
};

class LightGBMTextModel {
 public:
  bool load(const std::string& path) {
    trees_.clear();
    feature_names_.clear();
    last_error_.clear();

    std::ifstream in(path);
    if (!in) {
      last_error_ = "failed to open model: " + path;
      return false;
    }

    Tree current;
    bool in_tree = false;
    std::string line;
    while (std::getline(in, line)) {
      if (line.rfind("feature_names=", 0) == 0) {
        feature_names_ = parseStrings(afterEquals(line));
      } else if (line.rfind("Tree=", 0) == 0) {
        if (in_tree)
          trees_.push_back(std::move(current));
        current = Tree{};
        in_tree = true;
      } else if (in_tree && line.rfind("split_feature=", 0) == 0) {
        current.split_feature = parseInts(afterEquals(line));
      } else if (in_tree && line.rfind("threshold=", 0) == 0) {
        current.threshold = parseDoubles(afterEquals(line));
      } else if (in_tree && line.rfind("left_child=", 0) == 0) {
        current.left_child = parseInts(afterEquals(line));
      } else if (in_tree && line.rfind("right_child=", 0) == 0) {
        current.right_child = parseInts(afterEquals(line));
      } else if (in_tree && line.rfind("leaf_value=", 0) == 0) {
        current.leaf_value = parseDoubles(afterEquals(line));
      }
    }
    if (in_tree)
      trees_.push_back(std::move(current));

    if (feature_names_.empty()) {
      last_error_ = "model has no feature_names";
      return false;
    }
    if (trees_.empty()) {
      last_error_ = "model has no trees";
      return false;
    }
    for (std::size_t i = 0; i < trees_.size(); ++i) {
      const auto& t = trees_[i];
      const std::size_t n = t.split_feature.size();
      if (t.threshold.size() != n || t.left_child.size() != n ||
          t.right_child.size() != n || t.leaf_value.empty()) {
        last_error_ = "malformed tree at index " + std::to_string(i);
        return false;
      }
    }
    return true;
  }

  bool loaded() const { return !trees_.empty() && !feature_names_.empty(); }
  const std::string& lastError() const { return last_error_; }

  double predictRaw(const RemeshScenarioFeatures& f) const {
    const auto x = makeFeatureVector(f);
    double raw = 0.0;
    for (const auto& tree : trees_)
      raw += evalTree(tree, x);
    return raw;
  }

  double predictProbability(const RemeshScenarioFeatures& f) const {
    const double raw = predictRaw(f);
    if (raw >= 40.0)
      return 1.0;
    if (raw <= -40.0)
      return 0.0;
    return 1.0 / (1.0 + std::exp(-raw));
  }

 private:
  struct Tree {
    std::vector<int> split_feature;
    std::vector<double> threshold;
    std::vector<int> left_child;
    std::vector<int> right_child;
    std::vector<double> leaf_value;
  };

  std::vector<std::string> feature_names_;
  std::vector<Tree> trees_;
  std::string last_error_;

  static std::string afterEquals(const std::string& line) {
    const auto pos = line.find('=');
    return pos == std::string::npos ? std::string{} : line.substr(pos + 1);
  }

  static std::vector<std::string> parseStrings(const std::string& text) {
    std::vector<std::string> out;
    std::istringstream iss(text);
    std::string item;
    while (iss >> item)
      out.push_back(item);
    return out;
  }

  static std::vector<int> parseInts(const std::string& text) {
    std::vector<int> out;
    std::istringstream iss(text);
    int v = 0;
    while (iss >> v)
      out.push_back(v);
    return out;
  }

  static std::vector<double> parseDoubles(const std::string& text) {
    std::vector<double> out;
    std::istringstream iss(text);
    double v = 0.0;
    while (iss >> v)
      out.push_back(v);
    return out;
  }

  static bool contains(const std::string& s, const char* token) {
    return s.find(token) != std::string::npos;
  }

  static double countChar(const std::string& s, char c) {
    return static_cast<double>(std::count(s.begin(), s.end(), c));
  }

  static double finiteOrMinusOne(double v) {
    return std::isfinite(v) ? v : -1.0;
  }

  static double featureValue(const RemeshScenarioFeatures& f,
                             const std::string& name) {
    const std::string& s = f.scenario;
    if (name == "op") return static_cast<double>(f.op);
    if (name == "reason") return static_cast<double>(f.reason);
    if (name == "edge_category") return static_cast<double>(f.edge_category);
    if (name == "len") return finiteOrMinusOne(f.len);
    if (name == "len_ratio") return finiteOrMinusOne(f.len_ratio);
    if (name == "theta") return finiteOrMinusOne(f.theta);
    if (name == "theta_ratio") return finiteOrMinusOne(f.theta_ratio);
    if (name == "scenario_len") return static_cast<double>(s.size());
    if (name == "scenario_has_P") return contains(s, "P") ? 1.0 : 0.0;
    if (name == "scenario_has_C") return contains(s, "C") ? 1.0 : 0.0;
    if (name == "scenario_has_F") return contains(s, "F") ? 1.0 : 0.0;
    if (name == "scenario_has_Fv") return contains(s, "Fv") ? 1.0 : 0.0;
    if (name == "scenario_has_Fv2") return contains(s, "Fv2") ? 1.0 : 0.0;
    if (name == "scenario_has_Fd") return contains(s, "Fd") ? 1.0 : 0.0;
    if (name == "scenario_has_Fvs") return contains(s, "Fvs") ? 1.0 : 0.0;
    if (name == "scenario_has_Fv2s") return contains(s, "Fv2s") ? 1.0 : 0.0;
    if (name == "scenario_has_Fds") return contains(s, "Fds") ? 1.0 : 0.0;
    if (name == "scenario_has_S") return contains(s, "S") ? 1.0 : 0.0;
    if (name == "scenario_has_Ss") return contains(s, "Ss") ? 1.0 : 0.0;
    if (name == "scenario_smooth_count") return countChar(s, 'S');
    return -1.0;
  }

  std::vector<double> makeFeatureVector(const RemeshScenarioFeatures& f) const {
    std::vector<double> x;
    x.reserve(feature_names_.size());
    for (const auto& name : feature_names_)
      x.push_back(featureValue(f, name));
    return x;
  }

  static double evalTree(const Tree& tree, const std::vector<double>& x) {
    int node = 0;
    while (node >= 0) {
      if (static_cast<std::size_t>(node) >= tree.split_feature.size())
        return 0.0;
      const int feature_index = tree.split_feature[static_cast<std::size_t>(node)];
      const double value =
          (feature_index >= 0 && static_cast<std::size_t>(feature_index) < x.size())
              ? finiteOrMinusOne(x[static_cast<std::size_t>(feature_index)])
              : -1.0;
      const int child = (value <= tree.threshold[static_cast<std::size_t>(node)])
                            ? tree.left_child[static_cast<std::size_t>(node)]
                            : tree.right_child[static_cast<std::size_t>(node)];
      if (child < 0) {
        const int leaf = -child - 1;
        if (leaf >= 0 && static_cast<std::size_t>(leaf) < tree.leaf_value.size())
          return tree.leaf_value[static_cast<std::size_t>(leaf)];
        return 0.0;
      }
      node = child;
    }
    return 0.0;
  }
};

} // namespace BEMMeshPipeline::RemeshAI
