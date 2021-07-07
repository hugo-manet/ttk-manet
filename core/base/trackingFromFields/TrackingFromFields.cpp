#include "AbstractTriangulation.h"
#include "DataTypes.h"
#include "VineyardTracking.hpp"
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <ostream>

#include <TrackingFromFields.h>
#include <type_traits>
#include <unistd.h>
#include <vector>

/** If there is some (heavy ?) bug, consider having the max()
 *  computed only on the explicit structure, and not the implicit (link-cut)
 * one. I guess it should work and be less error-prone, but needs a more
 * systematic design.
 */

namespace ttk {

  double actualTime = 0.;
  double getVal(const MergeTreeLinkCutNode *const node) {
    return node->scalarStart
           + actualTime * (node->scalarEnd - node->scalarStart);
  }

  void setMax(MergeTreeLinkCutNode *&dest,
              const MergeTreeLinkCutNode *const option) {
    if(option == NULL)
      return;
    const double dVal = getVal(dest), oVal = getVal(option->actualMax);
    if(dVal < oVal || (dVal == oVal && (dest < option->actualMax)))
      dest = option->actualMax;
  }

  void updateMT(MergeTreeLinkCutNode *const x) {
    x->PT_max = x;
    for(auto son : x->PT_sons)
      setMax(x->PT_max, son);
  }

  void updateST(MergeTreeLinkCutNode *const x) {
    x->actualMax = x->PT_max;
    setMax(x->actualMax, x->ST_left);
    setMax(x->actualMax, x->ST_right);
  }

  void update(MergeTreeLinkCutNode *const x) {
    updateMT(x);
    updateST(x);
  }

  // /!\ does not conserve PT_sons/PT_parent ! That is done in
  // splay()
  void rotr(MergeTreeLinkCutNode *const x) {
    MergeTreeLinkCutNode *const y = x->ST_parent;
    MergeTreeLinkCutNode *const z = y->ST_parent;
    if((y->ST_left = x->ST_right))
      y->ST_left->ST_parent = y;
    x->ST_right = y, y->ST_parent = x;
    if((x->ST_parent = z)) {
      if(y == z->ST_left)
        z->ST_left = x;
      else
        z->ST_right = x;
    }
    updateST(y);
  }

  void rotl(MergeTreeLinkCutNode *const x) {
    MergeTreeLinkCutNode *const y = x->ST_parent;
    MergeTreeLinkCutNode *const z = y->ST_parent;
    if((y->ST_right = x->ST_left))
      y->ST_right->ST_parent = y;
    x->ST_left = y, y->ST_parent = x;
    if((x->ST_parent = z)) {
      if(y == z->ST_left)
        z->ST_left = x;
      else
        z->ST_right = x;
    }
    updateST(y);
  }

  void splay(MergeTreeLinkCutNode *const x) {
    if(x->ST_parent == NULL)
      return;
    while(x->ST_parent) {
      MergeTreeLinkCutNode *const y = x->ST_parent;
      if(y->ST_parent == 0) {
        if(y->PT_parent != NULL) {
          x->PT_parent = y->PT_parent;
          y->PT_parent = 0;
          x->PT_parent->PT_sons.erase(y);
          x->PT_parent->PT_sons.insert(x);
        }

        if(x == y->ST_left)
          rotr(x);
        else
          rotl(x);
      } else {
        MergeTreeLinkCutNode *const z = y->ST_parent;
        if(z->ST_parent == NULL) {
          if(z->PT_parent != NULL) {
            x->PT_parent = z->PT_parent;
            z->PT_parent = 0;
            x->PT_parent->PT_sons.erase(z);
            x->PT_parent->PT_sons.insert(x);
          }
        }

        if(y == z->ST_left) {
          if(x == y->ST_left)
            rotr(y), rotr(x);
          else
            rotl(x), rotr(x);
        } else {
          if(x == y->ST_right)
            rotl(y), rotl(x);
          else
            rotr(x), rotl(x);
        }
      }
    }
    updateST(x);
  }

  MergeTreeLinkCutNode *access(MergeTreeLinkCutNode *const x) {
    splay(x);
    if(x->ST_right) {
      // cut tree
      x->ST_right->PT_parent = x;
      x->PT_sons.insert(x->ST_right);
      x->ST_right->ST_parent = 0;
      x->ST_right = 0;
      update(x); // TODO maybe just update the PT_max ?
      // Maybe also check the actualMax, because this might change with the new
      // time ?
    }

    MergeTreeLinkCutNode *last = x;
    while(x->PT_parent) {
      MergeTreeLinkCutNode *const y = x->PT_parent;
      last = y;
      splay(y);
      if(y->ST_right) {
        y->ST_right->PT_parent = y;
        y->PT_sons.insert(y->ST_right);
        y->ST_right->ST_parent = 0;
      }
      y->ST_right = x;
      x->ST_parent = y;
      y->PT_sons.erase(x);
      x->PT_parent = 0;
      update(y);
      splay(x);
    }

    return last;
  }

  MergeTreeLinkCutNode *root(MergeTreeLinkCutNode *x) {
    access(x);
    while(x->ST_left)
      x = x->ST_left;
    splay(x);
    return x;
  }

  void cut(MergeTreeLinkCutNode *const x) {
    access(x);
    // assert(x->ST_left != NULL) // because x had a parent.
    // assert(x->MT_parent != NULL) // because x had a parent.
    x->ST_left->ST_parent = 0;
    x->ST_left = 0;

    x->MT_parent->MT_sons.erase(x);
    x->MT_parent = NULL;

    update(x);
  }

  void link(MergeTreeLinkCutNode *const x, MergeTreeLinkCutNode *const y) {
    access(x);
    access(y);
    // assert(x->ST_left == NULL) // because x had no parent.
    x->MT_parent = y;
    y->MT_sons.emplace(x);

    x->ST_left = y;
    y->ST_parent = x;
    update(x);
  }

  MergeTreeLinkCutNode *lca(MergeTreeLinkCutNode *const x,
                            MergeTreeLinkCutNode *const y) {
    access(x);
    return access(y);
  }

  EventQueue swapQueue;

  double crossTime(MergeTreeLinkCutNode *son, MergeTreeLinkCutNode *parent) {
    double crossing = -1.;
    if(son->scalarEnd - son->scalarStart - parent->scalarEnd
         + parent->scalarStart
       != 0.)
      crossing = (son->scalarStart - parent->scalarStart)
                 / (parent->scalarEnd - parent->scalarStart - son->scalarEnd
                    + son->scalarStart);

    if(crossing < 0. || crossing > 2.)
      crossing = 2.; // never gonna happen

    return crossing;
  }
  void createTreeSwapEvent(MergeTreeLinkCutNode *const son,
                           MergeTreeLinkCutNode *const parent) {
    double crossing = crossTime(son, parent);
    if(son->scalarEnd > parent->scalarEnd)
      crossing = 2.;

    if(crossing > 0. || (crossing == 0. && son->scalarEnd > parent->scalarEnd))
      swapQueue.insert({crossing, SwapEvent{son, parent, false, NULL
#ifndef NODEBUG
                                            ,
                                            actualTime, NULL
#endif
                                  }});
  }
  void newLink(MergeTreeLinkCutNode *const son,
               MergeTreeLinkCutNode *const parent) {
    link(son, parent);

    createTreeSwapEvent(son, parent);
  }

  /** WARNING : Must be done with actualTime < (!=) the time of the event.
   *  If not, you'll get a crossTime of +inf
   *  @return true if time was modified
   */
  bool updatePairEventTime(NodePair *ofThePair) {
    access(ofThePair->saddle);

    double newTime = crossTime(ofThePair->saddle->PT_max, ofThePair->max);
    if(ofThePair->saddle->PT_max->scalarEnd >= ofThePair->max->scalarEnd)
      newTime = 2.;

    if(newTime != ofThePair->itToEvent->first) {
#ifndef NODEBUG
      if(newTime <= actualTime) {
        auto val = &(cerr << "Updated a time into a time earlier than now !");
        cerr << val->bad();
        val = &((*val) << endl);
      }
#endif

#ifdef CPLUSPLUS17 // It's not. It should.
      auto handle = swapQueue.extract(ofThePair->itToEvent);
      handle.key() = newTime;
      ofThePair->itToEvent = swapQueue.insert(move(handle));
#else
      // Grrrrr. It's a whole copy for nothing :'(
      auto newIt = swapQueue.emplace(newTime, ofThePair->itToEvent->second);
      swapQueue.erase(ofThePair->itToEvent);
      ofThePair->itToEvent = newIt;
#endif
#ifndef NODEBUG
      ofThePair->itToEvent->second.timeOfLastUpdate = actualTime;
      ofThePair->itToEvent->second.wasLosingTo = ofThePair->saddle->PT_max;
#endif
      return true;
    }
    return false;
  }
  void MergeTreeLinkCutNode::swapWithSon(MergeTreeLinkCutNode *const son) {
    if(MT_sons.count(son) == 0)
      return; // This event is outdated

    const bool isLocal = upperLink.erase(son) == 1;
    if(isLocal) {
      son->upperLink.insert(this);
    }

    const int nbSons = MT_sons.size();
    const int nbGrandsons = son->MT_sons.size();

    if(nbSons == 1 && nbGrandsons == 1 && !isLocal) {
      // Optimization : Regular nodes swap. Times goes from 3h to 40m with this.
      // TODO see in benchmark if it's still the slowest part,
      // and if optimizing it further (ie drop access/splay or similar)
      // might make the code run faster
      access(son);
      splay(this);

      son->ST_left = this->ST_left;
      this->ST_left = NULL;
      son->ST_left->ST_parent = son;

      this->ST_right = NULL;
      son->ST_parent = NULL;

      son->ST_right = this;
      this->ST_parent = son;

      auto grandson = *son->MT_sons.begin();
      auto father = this->MT_parent;
      auto grandsonPT = *son->PT_sons.begin();

      grandsonPT->PT_parent = this;
      this->PT_sons.insert(grandsonPT);
      son->PT_sons.clear();

      father->MT_sons.erase(this);
      father->MT_sons.insert(son);
      son->MT_sons.clear();
      son->MT_sons.insert(this);
      this->MT_sons.clear();
      this->MT_sons.insert(grandson);
      son->MT_parent = father;
      this->MT_parent = son;
      grandson->MT_parent = this;

      update(this);
      update(son);

      createTreeSwapEvent(son, father);
      createTreeSwapEvent(grandson, this);

      return;
    }

    MergeTreeLinkCutNode *const myParent = this->MT_parent;
    cut(this);
    cut(son);

    // We need a temporary copy because iterators will be invalidated when we
    // cut.
    std::vector<MergeTreeLinkCutNode *> listOfGrandsons(
      son->MT_sons.begin(), son->MT_sons.end());
    for(auto grandson : listOfGrandsons)
      cut(grandson);

    std::set<MergeTreeLinkCutNode *> seenRoots;
    for(auto neigh : upperLink)
      seenRoots.insert(root(neigh));

    for(auto rt : seenRoots) {
      auto idx = std::find(listOfGrandsons.begin(), listOfGrandsons.end(), rt);
      if(idx != listOfGrandsons.end()) {
        newLink(*idx, this);
        listOfGrandsons.erase(idx);
      }
    }
    for(auto unseen : listOfGrandsons)
      link(unseen, son);

    newLink(son, myParent);
    link(this, son); // not a newLink because they've just swap

    const int nbGrandsonKept = son->MT_sons.size() - 1;

    // treat simple cases for PD following
    if(nbGrandsons == 0) {
      // Maximum transfer or destruction
      if(nbSons == 1) {
        this->pairOfMax = son->pairOfMax;
        this->pairOfMax->max = this;
        updatePairEventTime(this->pairOfMax);
      } else {
        swapQueue.erase(son->pairOfMax->itToEvent);
        delete son->pairOfMax;
      }

      son->pairOfMax = NULL;
      return;
    }
    if(nbGrandsons == 1 && nbGrandsonKept == 0)
      return; // node stays regular
    if(nbSons == 1 && nbGrandsonKept == nbGrandsons) {
      // New local max
      this->pairOfMax = new NodePair{this, son, -1, swapQueue.end()};
      this->pairOfMax->itToEvent
        = swapQueue.insert({2.,
                            {NULL, NULL, true, this->pairOfMax
#ifndef NODEBUG
                             ,
                             actualTime, son->PT_max
#endif
                            }});

      updatePairEventTime(this->pairOfMax);
      return;
    }

    // If we haven't returned yet, we need to rebuild saddle-maxima links
    // Here, all of *this and *son 's MT children are in a different splay
    // trees, so their actualMax are right !

    std::vector<MergeTreeLinkCutNode *> losingPairs;
    MergeTreeLinkCutNode *winningPair = NULL;
    for(auto branch : this->PT_sons) {
      if(winningPair == NULL) // Exactly the first iteration
        winningPair = branch;
      else {
        losingPairs.push_back(branch);
        if(getVal(branch->actualMax) > getVal(winningPair->actualMax))
          std::swap(winningPair, losingPairs.back());
      }
    }
    for(auto branch : losingPairs) {
      branch->actualMax->pairOfMax->saddle = this;
      updatePairEventTime(branch->actualMax->pairOfMax);
    }

    losingPairs.clear();

    for(auto branch : son->PT_sons) {
      losingPairs.push_back(branch);
      if(getVal(branch->actualMax) > getVal(winningPair->actualMax))
        std::swap(winningPair, losingPairs.back());
    }
    for(auto branch : losingPairs) {
      branch->actualMax->pairOfMax->saddle = son;
      updatePairEventTime(branch->actualMax->pairOfMax);
    }

    return;
  }

  std::vector<MergeTreeLinkCutNode> buildTree(AbstractTriangulation *grid,
                                              double *scalarsStart,
                                              double *scalarsEnd) {
    actualTime = 0.;
    swapQueue.clear();
    grid->preconditionVertexNeighbors();
    SimplexId nbNodes = grid->getNumberOfVertices();

#ifndef NONOISE
    srand(42); // Try 42 for bug at time 0.00996
    for(int i = 0; i < nbNodes; ++i) {
      scalarsStart[i] += (rand() - RAND_MAX / 2.) / (RAND_MAX * 1000.);
      scalarsEnd[i] += (rand() - RAND_MAX / 2.) / (RAND_MAX * 1000.);
    }
#endif

    std::vector<SimplexId> insertOrder(nbNodes);
    std::vector<MergeTreeLinkCutNode> treeData(nbNodes + 2);
    for(SimplexId i = 0; i < nbNodes; ++i)
      insertOrder[i] = i;
    auto compareIdx = [&](SimplexId a, SimplexId b) -> bool {
      return (scalarsStart[a] != scalarsStart[b])
               ? scalarsStart[a] < scalarsStart[b]
               : a < b;
    };
    std::sort(insertOrder.begin(), insertOrder.end(), compareIdx);
    std::reverse(insertOrder.begin(), insertOrder.end());

    double globMin = -1., globMax = 1.;
    auto globRoot = &treeData[nbNodes];
    auto globTop = &treeData[nbNodes + 1];
    for(auto idNode : insertOrder) {

      globMin = std::min(globMin, scalarsStart[idNode]);
      globMin = std::min(globMin, scalarsEnd[idNode]);
      globMax = std::max(globMax, scalarsStart[idNode]);
      globMax = std::max(globMax, scalarsEnd[idNode]);
      auto that = &treeData[idNode];
#ifndef NODEBUG
      std::cout << "inserting " << idNode << endl;
      that->numForDebug = idNode;
#endif
      that->scalarStart = scalarsStart[idNode];
      that->scalarEnd = scalarsEnd[idNode];
      that->actualMax = that->PT_max = that;
      SimplexId neigh;
      for(int iNeigh = 0; iNeigh < grid->getVertexNeighborNumber(idNode);
          ++iNeigh) {
        grid->getVertexNeighbor(idNode, iNeigh, neigh);
        if(compareIdx(idNode, neigh)) {

#ifndef NODEBUG
          std::cout << "iNeigh " << iNeigh
                    << " is to be inserted (and maybe linked)" << endl;
#endif
          that->upperLink.insert(&treeData[neigh]);
          auto neighRoot = root(&treeData[neigh]);
          if(neighRoot != that) {
#ifndef NODEBUG
            std::cout << "it is linked" << endl;
#endif

            newLink(neighRoot, that);
          }
        }
      }
      if(that->MT_sons.size() == 0) {
        that->pairOfMax = new NodePair{that, globRoot, idNode, swapQueue.end()};
      } else if(that->MT_sons.size() > 1) {
        access(that);
        auto maxSon = that->PT_max;
        for(auto son : that->MT_sons) {
          access(son);
          if(son->PT_max != maxSon) {
            auto np = son->PT_max->pairOfMax;
            np->saddle = that;
            np->itToEvent
              = swapQueue.insert({crossTime(son->PT_max, that->PT_max),
                                  {NULL, NULL, true, np
#ifndef NODEBUG
                                   ,
                                   actualTime, that->PT_max
#endif
                                  }});
          }
        }
      }
    }
    // Build sentinels
    globTop->scalarStart = globTop->scalarEnd = globMax * 1.01;
    globTop->actualMax = globTop->PT_max = globTop;
    globRoot->scalarStart = globRoot->scalarEnd
      = globMin * 1.01; // globMin < -1
    globRoot->actualMax = globRoot->PT_max = globRoot;

    link(globTop, globRoot);

    auto oldRoot = &treeData[insertOrder[nbNodes - 1]];
    link(oldRoot, globRoot);
    access(oldRoot);
    auto np = oldRoot->PT_max->pairOfMax;
    np->saddle = globRoot;
    np->itToEvent = swapQueue.insert({2.,
                                      {NULL, NULL, true, np
#ifndef NODEBUG
                                       ,
                                       actualTime, NULL
#endif
                                      }});

    return treeData;
  }

  // maybe this will be over some value yet to be inserted...
  // Not for maxSwap anymore, but maybe treeSwap
  void setActuTime(double eventTime) {
    const double oldActualTime[[maybe_unused]] = actualTime;

    actualTime = eventTime + (swapQueue.begin()->first - eventTime) / 16777216.;
    if(actualTime == eventTime) {
      actualTime = eventTime + (swapQueue.begin()->first - eventTime) / 2048.;
      if(actualTime == eventTime) {
        actualTime = eventTime + (swapQueue.begin()->first - eventTime) / 2.;
        if(actualTime == eventTime) {
          std::cerr << "ALERT : possible numerical instability at time "
                    << actualTime << std::endl;
          // TODO +eps is a anticipated hotfix ; to be done better someday
          actualTime += 0.0001;
        }
      }
    }
    if(ceil(actualTime * 10000) != ceil(oldActualTime * 10000))
      std::cout << "Time : " << actualTime << endl;
  }

  void loopQueue() {
#ifndef NODEBUG
    std::cout << "queue size " << swapQueue.size() << " starting at "
              << swapQueue.begin()->first << endl;
#endif
    swapQueue.insert({1., {NULL, NULL, false, NULL}});

    cerr << std::setprecision(15);
    // TODO switch back to 1.
    while(swapQueue.begin()->first < 1. / 1000.) {
      auto eventIt = swapQueue.begin();
      auto event = eventIt->second;
      double eventTime = eventIt->first;

      while(eventIt->first == eventTime && !event.isMaxSwap
            && event.rooting == eventIt->second.rooting
            && event.leafing == eventIt->second.leafing) {
        swapQueue.erase(eventIt);
        eventIt = swapQueue.begin();
      } // erase duplicates. Don't erase maxSwaps, those are persistent.

#ifndef NODEBUG
      if(event.isMaxSwap)
        std::cout << "queue size " << swapQueue.size() << ", event "
                  << (event.isMaxSwap ? "maxSwap" : "treeSwap") << " with time "
                  << eventTime << endl;
#endif

      if(event.isMaxSwap) {
        MergeTreeLinkCutNode *theSaddle = event.pairSwap->saddle;
        access(theSaddle);

#ifndef NODEBUG
        if(theSaddle->MT_sons.size() > 2) {
          cerr << "Multiple saddle pair swap !" << endl;
          double lol = crossTime(theSaddle->PT_max, event.pairSwap->max);
          cerr << "Has " << theSaddle->MT_sons.size() << " sons, swap time "
               << lol << endl;
        }
#endif

        MergeTreeLinkCutNode *theOldMax = theSaddle->PT_max;

        // TODO BUG This should explode.
        if(updatePairEventTime(event.pairSwap))
          continue;

#ifndef NODEBUG
        std::cout << "accomplish it" << endl;
#endif
        MergeTreeLinkCutNode *newMax = event.pairSwap->max;

        MergeTreeLinkCutNode *newSad = theOldMax->pairOfMax->saddle;
        access(newSad);
#ifndef NODEBUG
        if(theOldMax->pairOfMax->saddle->MT_sons.size() > 2) {
          cerr << "pair swap into a multiple saddle !" << endl;
          double lol
            = crossTime(theOldMax->pairOfMax->saddle->PT_max, theOldMax);
          double lol2 = crossTime(theOldMax->pairOfMax->saddle->PT_max, newMax);
          cerr << "Has " << theSaddle->MT_sons.size() << " sons, swap time was "
               << lol << " and will be " << lol2 << endl;
        }

        if(theOldMax == newMax) {
          std::cerr << "Tried swapping a node with himself ! WTF" << std::endl;
          continue;
        }
#endif
        // Maybe swap values and not ptrs ?
        std::swap(theOldMax->pairOfMax->saddle, newMax->pairOfMax->saddle);
        std::swap(
          theOldMax->pairOfMax->idFirstMax, newMax->pairOfMax->idFirstMax);
        updatePairEventTime(newMax->pairOfMax);
        setActuTime(eventTime);
        access(theSaddle);
        access(newSad);
        swapQueue.erase(theOldMax->pairOfMax->itToEvent);
        theOldMax->pairOfMax->itToEvent
          = swapQueue.insert({2., {NULL, NULL, true, theOldMax->pairOfMax}});
      } else {
        setActuTime(eventTime);
        event.leafing->swapWithSon(event.rooting);
      }
#ifndef NODEBUG
      if(actualTime >= swapQueue.begin()->first) {
        std::cerr << "ALERT : possible time inconsistency, inserted event "
                     "before actual time."
                  << endl;
        auto nextEv = swapQueue.begin()->second;
        if(nextEv.isMaxSwap)
          std::cerr << "Next event (max) :" << nextEv.pairSwap->max << " ^, "
                    << nextEv.pairSwap->saddle << " sad ; at time "
                    << swapQueue.begin()->first << endl;
        else
          std::cerr << "Next event (tree) :" << nextEv.leafing->numForDebug
                    << " ^, " << nextEv.rooting->numForDebug << " v ; at time "
                    << swapQueue.begin()->first << endl;

        if(event.isMaxSwap)
          std::cerr << "Old event (max) :" << event.pairSwap->max << " ^, "
                    << event.pairSwap->saddle << " sad ; at time " << eventTime
                    << endl;
        else
          std::cerr << "Old event (tree) :" << event.leafing->numForDebug
                    << " ^, " << event.rooting->numForDebug << " v ; at time "
                    << eventTime << endl;
        if(eventTime == swapQueue.begin()->first)
          cerr << "Time is the same" << endl;
      }
#endif
    }
  }

} // namespace ttk
