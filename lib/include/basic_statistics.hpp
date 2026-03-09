#ifndef basic_statistics_H
#define basic_statistics_H

#include "basic_vectors.hpp"
#include <cmath>
#include <vector>
/* ------------------------------------------------------ */
inline double Variance(const V_d &V) {
  if (V.empty())
    return 0.;
  double sum_of_squares = 0., total = 0., sum_of_variance = 0.;
  for (const auto &w : V)
    total += w;
  double mean = total / V.size();
  for (const auto &w : V)
    sum_of_variance += std::pow(w - mean, 2);
  return sum_of_variance / V.size();
};

template <std::size_t N> double Variance(const std::array<double, N> &V) {
  double sum_of_squares = 0., total = 0., sum_of_variance = 0.;
  for (const auto &w : V)
    total += w;
  double mean = total / N;
  for (const auto &w : V)
    sum_of_variance += std::pow(w - mean, 2);
  return sum_of_variance / N;
}

template <typename T> double RootMeanSquare(const T &TUPLE) {
  double ret = 0, div = 0;
  for_each(TUPLE, [&](const auto &t) {
    ret += t * t;
    div += 1.;
  });
  return std::sqrt(ret / div);
};

/* ------------------------------------------------------ */

/*DOC_EXTRACT

与えられたデータの数によって`sturges = log2(n) + 1`となるようにbinの数を決める．

| Vector             | Length Description                    | Meaning                                                    |
|--------------------|--------------------------------------|------------------------------------------------------------|
| `data`             | Length = n (total data points)       | Sorted raw data points.                                    |
| `bins`             | Length = sturges                     | Groups of data points within specific ranges.              |
| `count`            | Length = sturges                     | Number of data points in each bin.                         |
| `cumulative_count` | Length = sturges                     | Running total of data points up to each bin.               |
| `interval`         | Length = sturges + 1                 | Boundary values of each bin.                               |
| `mid_interval`     | Length = sturges                     | Midpoint value of each bin's range.                        |
| `diff`             | Length = sturges - 1                 | Difference in counts between consecutive bins.             |
| `diffdiff`         | Length = sturges - 2                 | Difference between consecutive `diff` values.              |
| `cumulative`       | Length = sturges                     | Cumulative distribution percentage up to each bin.         |

- **Data (`data`)**: The starting point, representing the original, sorted dataset.
- **Bins (`bins`)**: Containers that divide the entire range of data into equal intervals.
- **Count (`count`)**: Reflects the frequency of data points in each bin.
- **Cumulative Count (`cumulative_count`)**: Shows the accumulation of data points up to each bin, useful for understanding data distribution.
- **Interval (`interval`)**: Defines the range of each bin. It's crucial for determining which bin a data point belongs to.
- **Mid Interval (`mid_interval`)**: Useful for representing each bin by a single value.
- **Differences (`diff`, `diffdiff`)**: Provide insights into the rate of change in data frequency between bins.
- **Cumulative (`cumulative`)**: Expresses the cumulative proportion of data points up to each bin, helpful in statistical analyses.

*/

struct Histogram {
  V_d data;    // ソートした生データ
  int size;    // 全dataの長さ
  int sturges; // dataの数に応じてbinの幅を決める
  VV_d bins;
  double bin_width;
  V_i count;            // 角瓶の中にあるデータの数
  V_i cumulative_count; // 角瓶の中にあるデータの数を初めから累積してった結果
  V_d cumulative;
  V_i diff;
  V_i diffdiff;
  V_d interval;
  V_d mid_interval;

  static V_d sort(const V_d &dataIN) {
    V_d dataOUT = dataIN;
    std::sort(dataOUT.begin(), dataOUT.end(), [](const double a, const double b) { return a < b; });
    return dataOUT;
  }

  /*
  bins.size() = s = 5
         min <--------------------------->  max
  min+w*i ->| i=0 |  1  |  2  |  3  |  4  |<- max=min+w*s
            |<-w->|
            |     |     |<-- intervals
  w = bin_width = (max - min)/s

  i-th bin contains values between the range [min+w*i, min+w*(i+1)]
  or
  bin[(int)((value-min)/w)]

  v = w*i + min + [0,w]
  -> v - min - [0,w] = w*i
  -> (v-min)/w - [0,1] = i
  -> (int)((v-min)/w - [0,1]) = i
  -> (int)((v-min)/w) = i

  NOTE: v = max
  -> (int)((max-min)/w) = i
  -> (int)(s) = i
  -> s = i. out of bounds!
  -> if i == s -> i = s-1
  */

  void set(const V_d &dataIN) {
    this->data = (this->sort(dataIN));
    this->size = dataIN.size();
    this->sturges = (int)(std::log2(this->size) + 3);
    this->bins = VV_d(this->sturges, V_d{});
    this->bin_width = (*this->data.rbegin() - *this->data.begin()) / sturges;
    auto min_data = *this->data.begin();
    int n;
    for (const auto &d : this->data) {
      n = (int)((d - min_data) / this->bin_width);
      if (n == bins.size())
        n -= 1;
      this->bins[n].emplace_back(d);
    }
    this->diff.resize(this->bins.size() - 1, 0);
    for (auto i = 0; i < this->bins.size() - 1; i++)
      this->diff[i] = this->bins[i + 1].size() - this->bins[i].size();
    //
    this->diffdiff.resize(this->diff.size() - 1, 0);
    for (auto i = 0; i < this->diff.size() - 1; i++)
      this->diffdiff[i] = this->diff[i + 1] - this->diff[i];
    //
    this->interval.resize(this->bins.size() + 1, 0);
    for (auto i = 0; i < this->interval.size(); i++)
      this->interval[i] = i * bin_width + min_data;
    //
    this->mid_interval.resize(this->bins.size(), 0);
    this->count.resize(this->bins.size(), 0);
    for (auto i = 0; i < this->mid_interval.size(); i++) {
      this->mid_interval[i] = (this->interval[i] + this->interval[i + 1]) / 2.;
      this->count[i] = this->bins[i].size();
    }
    //
    this->cumulative.resize(this->bins.size());
    this->cumulative_count.resize(this->bins.size());
    this->cumulative_count[0] = this->count[0];
    this->cumulative[0] = (double)this->cumulative_count[0] / (double)this->size;
    for (auto i = 1; i < this->cumulative_count.size(); i++) {
      this->cumulative_count[i] = this->cumulative_count[i - 1] + this->count[i];
      this->cumulative[i] = (double)this->cumulative_count[i] / (double)this->size;
    }
  };

  Histogram() {};
  Histogram(const V_d &dataIN) : data(this->sort(dataIN)), size(dataIN.size()), sturges((int)(std::log2(this->size) + 3)), bins(VV_d(this->sturges, V_d{})), bin_width((*this->data.rbegin() - *this->data.begin()) / sturges) {
    // std::cout << "data " << this->data << std::endl;
    // std::cout << "size " << this->size << std::endl;
    // std::cout << "sturges " << this->sturges << std::endl;
    // std::cout << "bins " << this->bins << std::endl;
    // std::cout << "bin_width " << this->bin_width << std::endl;

    auto min_data = *this->data.begin();
    int n;
    for (const auto &d : this->data) {
      // n = floor((d - min_data) / this->bin_width);
      // this->bins[n >= bins.size() ? bins.size() - 1 : n].emplace_back(d);
      //
      n = (int)((d - min_data) / this->bin_width);
      if (n == bins.size())
        n -= 1;
      this->bins[n].emplace_back(d);
    }

    this->diff.resize(this->bins.size() - 1, 0);
    for (auto i = 0; i < this->bins.size() - 1; i++)
      this->diff[i] = this->bins[i + 1].size() - this->bins[i].size();
    //
    this->diffdiff.resize(this->diff.size() - 1, 0);
    for (auto i = 0; i < this->diff.size() - 1; i++)
      this->diffdiff[i] = this->diff[i + 1] - this->diff[i];
    //
    this->interval.resize(this->bins.size() + 1, 0);
    for (auto i = 0; i < this->interval.size(); i++)
      this->interval[i] = i * bin_width + min_data;
    //
    this->mid_interval.resize(this->bins.size(), 0);
    this->count.resize(this->bins.size(), 0);
    for (auto i = 0; i < this->mid_interval.size(); i++) {
      this->mid_interval[i] = (this->interval[i] + this->interval[i + 1]) / 2.;
      this->count[i] = this->bins[i].size();
    }
    //
    this->cumulative.resize(this->bins.size());
    this->cumulative_count.resize(this->bins.size());
    this->cumulative_count[0] = this->count[0];
    this->cumulative[0] = (double)this->cumulative_count[0] / (double)this->size;
    for (auto i = 1; i < this->cumulative_count.size(); i++) {
      this->cumulative_count[i] = this->cumulative_count[i - 1] + this->count[i];
      this->cumulative[i] = (double)this->cumulative_count[i] / (double)this->size;
    }
  };
};

#endif