/// \ingroup base
/// \class ttk::MultiTrackSource
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.09.2019
///
/// This filter creates an MultiTrackSource with a specified radius, center, and
/// number of subdivisions.

#pragma once

// base code includes
#include <Debug.h>

// std includes

namespace ttk {

  class MultiTrackSource : virtual public Debug {

  public:
    MultiTrackSource() {
      this->setDebugMsgPrefix("MultiTrackSource");
    };
    ~MultiTrackSource() override{};

  private:
  };
} // namespace ttk
