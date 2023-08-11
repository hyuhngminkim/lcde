#ifndef LCDE_MODEL_H_
#define LCDE_MODEL_H_

#include "builder.h"

namespace lcde {

template <typename KeyType>
class Model {
  public:
  Model(const Model&) = delete;
  Model() = default;

  std::vector<int> rootMap;
  std::vector<typename Builder<KeyType>::cdfPoint> globalLayout;
  int LRANGE, RRANGE, fanout;
  double slope, intercept;

  void finalize(const Builder<KeyType>& b) {
    rootMap = std::vector<int>(b.rootMap);
    globalLayout = std::vector<typename Builder<KeyType>::cdfPoint>(b.globalLayout);
    LRANGE = b.LRANGE;
    RRANGE = b.RRANGE;
    fanout = b.fanout;
    slope = static_cast<double>(b.slope);
    intercept = static_cast<double>(b.intercept);

    // for (auto i : globalLayout) {
    //   std::cout << i.knot << " " << static_cast<double>(i.cdf) / 200000000 << std::endl;
    // }
  }

  std::pair<size_t, size_t> find(const KeyType& key) const {
    long rank = slope * key + intercept;
    rank = std::max(0L, std::min(static_cast<long>(fanout - 1), rank));
    long nrank = std::min(rank + 1, static_cast<long>(fanout - 1));
  
    auto it = std::upper_bound(globalLayout.begin() + rootMap[rank], 
                               globalLayout.begin() + rootMap[nrank], key,
    [](const KeyType& k, const typename Builder<KeyType>::cdfPoint& c) {
      return k < c.knot;
    });

    // Perform binary search on the entire `globalLayout`
    // Search speed is slower than using `rootMap`
    // auto it = std::upper_bound(globalLayout.begin(), globalLayout.end(), key,
    //   [](const KeyType& k, const typename Builder<KeyType>::cdfPoint& c) {
    //     return k < c.knot;
    // });
    
    auto LITER = std::max(std::prev(it, LRANGE), globalLayout.begin());;
    auto RITER = std::min(std::next(it, RRANGE), globalLayout.end() - 1);

    return std::make_pair(LITER->cdf, RITER->cdf);
  }

  size_t getSize() const {
    size_t int_size = sizeof(int);
    size_t double_size = sizeof(double);
    size_t key_size = sizeof(KeyType);

    // typename Builder<KeyType>::cdfPoint* ptr;

    return int_size * (rootMap.size() + 3)
          + double_size * (globalLayout.size() + 2)
          + key_size * globalLayout.size();
  }
  
};

}

#endif // LCDE_MODEL_H_