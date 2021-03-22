#ifndef FUZZY_DISTANCE_MATRIX_HPP
#define FUZZY_DISTANCE_MATRIX_HPP

#include <cmath>
#include <limits>
#include <vector>

template <typename T>
struct TMatrix {
  size_t _nRows, _nCols;
  std::vector<T> _values;

  typename std::vector<T>::reference operator()(size_t row, size_t col) {
    return _values[row * _nCols + col];
  }
  typename std::vector<T>::const_reference operator()(size_t row,
                                                      size_t col) const {
    return _values[row * _nCols + col];
  }

  TMatrix(size_t nRows, size_t nCols, T val = T())
    : _nRows(nRows), _nCols(nCols), _values(nRows * nCols, val) {
  }

  TMatrix() = default;
};

using Matrix = TMatrix<double>;
using BoolMatrix = TMatrix<bool>;

// Only square matrices for now !
struct FuzzyDistanceMatrix {
  Matrix minorM;
  Matrix majorM;
  BoolMatrix isExact;
  size_t nbUpdates, nbCandidateUpdates;

  FuzzyDistanceMatrix(size_t nRows, size_t nCols)
    : minorM(nRows, nCols, 0.),
      majorM(nRows, nCols, std::numeric_limits<double>::infinity()),
      isExact(nRows, nCols, false), nbUpdates(0), nbCandidateUpdates(0) {
    // TODO only for square/extended square matrices ?
    for(size_t elem = 0; elem < nRows; ++elem)
      majorM(elem, elem) = 0.;
  }

  /* // Not used for now
  struct UpdateElement {
    double priority;
    size_t row, col, med;
    bool operator<(const UpdateElement& other) {
      return this->priority < other.priority;
    }
  }; // */

  void updateMajors(size_t row, size_t col, double val) {
    // TODO fast update with Dijkstra ?
    std::vector<size_t> enhancedRows, enhancedCols;
    enhancedRows.reserve(majorM._nRows);
    enhancedCols.reserve(majorM._nCols);

    for(size_t otherRow = 0; otherRow < majorM._nRows; ++otherRow)
      if(majorM(otherRow, col) > majorM(otherRow, row) + val)
        enhancedRows.push_back(otherRow);

    for(size_t otherCol = 0; otherCol < majorM._nCols; ++otherCol)
      if(majorM(row, otherCol) > val + majorM(col, otherCol))
        enhancedCols.push_back(otherCol);

    for(size_t ior = 0; ior < enhancedRows.size(); ++ior) {
      size_t otherRow = enhancedRows[ior];
      for(size_t otherCol : enhancedCols) {
        nbCandidateUpdates++;
        if(majorM(otherRow, otherCol)
           > majorM(otherRow, row) + val + majorM(col, otherCol))
          nbUpdates++;
        majorM(otherRow, otherCol)
          = std::min(majorM(otherRow, otherCol),
                     majorM(otherRow, row) + val + majorM(col, otherCol));
        majorM(otherCol, otherRow) = majorM(otherRow, otherCol);
      }
    }
  }

  void updateMinors(size_t row, size_t col, double val) {
    // TODO find similar fast updates as majors' Dijkstra ?
    std::vector<int> enhancedRows, enhancedCols;
    enhancedRows.reserve(minorM._nRows);
    enhancedCols.reserve(minorM._nCols);

    for(size_t otherCol = 0; otherCol < majorM._nCols; ++otherCol)
      if(minorM(row, otherCol) < val - majorM(col, otherCol))
        enhancedCols.push_back(otherCol);

    for(size_t otherRow = 0; otherRow < majorM._nRows; ++otherRow)
      if(minorM(otherRow, col) < val - majorM(otherRow, row))
        enhancedRows.push_back(otherRow);

    for(size_t ior = 0; ior < enhancedRows.size(); ++ior) {
      size_t otherRow = enhancedRows[ior];
      for(size_t otherCol : enhancedCols) {
        nbCandidateUpdates++;
        if(minorM(otherRow, otherCol)
           < val - majorM(otherRow, row) - majorM(col, otherCol))
          nbUpdates++;
        minorM(otherRow, otherCol)
          = std::max(minorM(otherRow, otherCol),
                     val - majorM(otherRow, row) - majorM(col, otherCol));
        minorM(otherCol, otherRow) = minorM(otherRow, otherCol);
      }
    }
  }

  void updateAll(size_t row, size_t col, double val) {
    updateMajors(row, col, val);
    updateMinors(row, col, val);
    isExact(row, col) = true;
  }

  void initOnCluster(Matrix graph, std::vector<size_t> &clusterIndexes) {
    for(auto center : clusterIndexes)
      for(size_t row = 0; row < majorM._nRows; ++row)
        for(size_t col = 0; col < majorM._nCols; ++col) {
          majorM(row, col) = std::min(
            majorM(row, col), graph(row, center) + graph(center, col));
          minorM(row, col)
            = std::max(minorM(row, col),
                       std::abs(graph(row, center) - graph(center, col)));
          isExact(row, center) = true;
          isExact(center, col) = true;
        }
  }
};

#endif // FUZZY_DISTANCE_MATRIX_HPP
