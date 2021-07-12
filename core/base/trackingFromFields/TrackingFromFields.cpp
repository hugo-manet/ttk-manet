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

namespace ttk {

  double getVal(const MergeTreeLinkCutNode *const node,
                const double actualTime) {
    return node->scalarStart
           + actualTime * (node->scalarEnd - node->scalarStart);
  }

  void setMax(MergeTreeLinkCutNode *&dest,
              const MergeTreeLinkCutNode *const option,
              const double actualTime) {
    if(option == NULL)
      return;
    const double dVal = getVal(dest, actualTime),
                 oVal = getVal(option->actualMax, actualTime);
    if(dVal < oVal || (dVal == oVal && (dest < option->actualMax)))
      dest = option->actualMax;
  }

  void updateMT(MergeTreeLinkCutNode *const x, const double actualTime) {
    x->PT_max = x;
    for(auto son : x->PT_sons)
      setMax(x->PT_max, son, actualTime);
  }

  void updateST(MergeTreeLinkCutNode *const x, const double actualTime) {
    x->actualMax = x->PT_max;
    setMax(x->actualMax, x->ST_left, actualTime);
    setMax(x->actualMax, x->ST_right, actualTime);
  }

  void update(MergeTreeLinkCutNode *const x, const double actualTime) {
    updateMT(x, actualTime);
    updateST(x, actualTime);
  }

  // /!\ does not conserve PT_sons/PT_parent ! That is done in
  // splay()
  void rotr(MergeTreeLinkCutNode *const x, const double actualTime) {
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
    updateST(y, actualTime);
  }

  void rotl(MergeTreeLinkCutNode *const x, const double actualTime) {
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
    updateST(y, actualTime);
  }

  void splay(MergeTreeLinkCutNode *const x, const double actualTime) {
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
          rotr(x, actualTime);
        else
          rotl(x, actualTime);
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
            rotr(y, actualTime), rotr(x, actualTime);
          else
            rotl(x, actualTime), rotr(x, actualTime);
        } else {
          if(x == y->ST_right)
            rotl(y, actualTime), rotl(x, actualTime);
          else
            rotr(x, actualTime), rotl(x, actualTime);
        }
      }
    }
    updateST(x, actualTime);
  }

  MergeTreeLinkCutNode *access(MergeTreeLinkCutNode *const x,
                               const double actualTime) {
    splay(x, actualTime);
    if(x->ST_right) {
      // cut tree
      x->ST_right->PT_parent = x;
      x->PT_sons.insert(x->ST_right);
      x->ST_right->ST_parent = 0;
      x->ST_right = 0;
      update(x, actualTime); // TODO maybe just update the PT_max ?
      // Maybe also check the actualMax, because this might change with the new
      // time ?
    }

    MergeTreeLinkCutNode *last = x;
    while(x->PT_parent) {
      MergeTreeLinkCutNode *const y = x->PT_parent;
      last = y;
      splay(y, actualTime);
      if(y->ST_right) {
        y->ST_right->PT_parent = y;
        y->PT_sons.insert(y->ST_right);
        y->ST_right->ST_parent = 0;
      }
      y->ST_right = x;
      x->ST_parent = y;
      y->PT_sons.erase(x);
      x->PT_parent = 0;
      update(y, actualTime);
      splay(x, actualTime);
    }

    return last;
  }

  MergeTreeLinkCutNode *root(MergeTreeLinkCutNode *x, const double actualTime) {
    access(x, actualTime);
    while(x->ST_left)
      x = x->ST_left;
    splay(x, actualTime);
    return x;
  }

  void cut(MergeTreeLinkCutNode *const x, const double actualTime) {
    access(x, actualTime);
    // assert(x->ST_left != NULL) // because x had a parent.
    // assert(x->MT_parent != NULL) // because x had a parent.
    x->ST_left->ST_parent = 0;
    x->ST_left = 0;

    x->MT_parent->MT_sons.erase(x);
    x->MT_parent = NULL;

    update(x, actualTime);
  }

  void link(MergeTreeLinkCutNode *const x,
            MergeTreeLinkCutNode *const y,
            const double actualTime) {
    access(x, actualTime);
    access(y, actualTime);
    // assert(x->ST_left == NULL) // because x had no parent.
    x->MT_parent = y;
    y->MT_sons.emplace(x);

    x->ST_left = y;
    y->ST_parent = x;
    update(x, actualTime);
  }

  MergeTreeLinkCutNode *lca(MergeTreeLinkCutNode *const x,
                            MergeTreeLinkCutNode *const y,
                            const double actualTime) {
    access(x, actualTime);
    return access(y, actualTime);
  }

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
                           MergeTreeLinkCutNode *const parent,
                           EventQueue &swapQueue,
                           const double actualTime) {
    if(parent->hasRecentSonElseInsert(son))
      return; // we've been here recently

    double crossing = crossTime(son, parent);
    if(son->scalarEnd > parent->scalarEnd || crossing >= 1.)
      return;
    if(crossing > 0. || (crossing == 0. && son->scalarEnd > parent->scalarEnd))
      swapQueue.insert({crossing, son, parent, false, NULL
#ifndef NODEBUG
                        ,
                        actualTime, NULL
#endif
      });
  }
  void newLink(MergeTreeLinkCutNode *const son,
               MergeTreeLinkCutNode *const parent,
               EventQueue &swapQueue,
               const double actualTime) {
    link(son, parent, actualTime);

    createTreeSwapEvent(son, parent, swapQueue, actualTime);
  }

  /** WARNING : Must be done with actualTime < (!=) the time of the event.
   *  If not, you'll get a crossTime of +inf
   *  @return true if time was modified
   */
  bool updatePairEventTime(NodePair *ofThePair,
                           EventQueue &swapQueue,
                           const double actualTime) {
    access(ofThePair->saddle, actualTime);

    double newTime = crossTime(ofThePair->saddle->PT_max, ofThePair->max);
    if(ofThePair->saddle->PT_max->scalarEnd >= ofThePair->max->scalarEnd)
      newTime = 2.;

    if(newTime != ofThePair->itToEvent->timestamp) {
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
      SwapEvent copy = *ofThePair->itToEvent;
      copy.timestamp = newTime;
      auto newIt = swapQueue.insert(copy).first;
      swapQueue.erase(ofThePair->itToEvent);
      ofThePair->itToEvent = newIt;
#endif
#ifndef NODEBUG
      ofThePair->itToEvent->timeOfLastUpdate = actualTime;
      ofThePair->itToEvent->wasLosingTo = ofThePair->saddle->PT_max;
#endif
      return true;
    }
    return false;
  }

  void repairNodeLinks(MergeTreeLinkCutNode *const node,
                       MergeTreeLinkCutNode *const other) {
    // If we see &*node in one of the links, we change it to &*other
    // and let repair(other,node) handle the other side
    if(node->ST_parent == node)
      node->ST_parent = other;
    else if(node->ST_parent) {
      if(node->ST_parent->ST_left == other)
        node->ST_parent->ST_left = node;
      else
        node->ST_parent->ST_right = node;
    }

    if(node->ST_left == node)
      node->ST_left = other;
    else if(node->ST_left)
      node->ST_left->ST_parent = node;

    if(node->ST_right == node)
      node->ST_right = other;
    else if(node->ST_right)
      node->ST_right->ST_parent = node;

    if(node->PT_parent == node)
      node->PT_parent = other;
    else if(node->PT_parent) {
      node->PT_parent->PT_sons.erase(other);
      node->PT_parent->PT_sons.emplace(node);
    }

    if(!node->PT_sons.empty()) {
      if(*node->PT_sons.begin() == node) {
        node->PT_sons.clear();
        node->PT_sons.emplace(other);
      } else
        (*node->PT_sons.begin())->PT_parent = node;
    }

    if(node->MT_parent == node)
      node->MT_parent = other;
    else if(node->MT_parent) {
      node->MT_parent->MT_sons.erase(other);
      node->MT_parent->MT_sons.emplace(node);
    }

    if(!node->MT_sons.empty()) {
      if(*node->MT_sons.begin() == node) {
        node->MT_sons.clear();
        node->MT_sons.emplace(other);
      } else
        (*node->MT_sons.begin())->MT_parent = node;
    }

    // we don't repair the values here : this can't be repaired easily
    // so we just update() the nodes afterwards.
  }
  void fastSwap(MergeTreeLinkCutNode *const that,
                MergeTreeLinkCutNode *const son,
                EventQueue &swapQueue,
                const double actualTime) {
    // swap every exterior links
    std::swap(that->ST_parent, son->ST_parent);
    std::swap(that->ST_left, son->ST_left);
    std::swap(that->ST_right, son->ST_right);
    std::swap(that->PT_parent, son->PT_parent);
    std::swap(that->PT_sons, son->PT_sons);
    std::swap(that->MT_parent, son->MT_parent);
    std::swap(that->MT_sons, son->MT_sons);
    std::swap(that->actualMax, son->actualMax);
    std::swap(that->PT_max, son->PT_max);

    repairNodeLinks(that, son);
    repairNodeLinks(son, that);
    update(son, actualTime);
    update(that, actualTime);
    update(son, actualTime); // Because son can depend on that

    createTreeSwapEvent(son, son->MT_parent, swapQueue, actualTime);
    createTreeSwapEvent(*that->MT_sons.begin(), that, swapQueue, actualTime);
  }
  void MergeTreeLinkCutNode::swapWithSon(MergeTreeLinkCutNode *const son,
                                         EventQueue &swapQueue,
                                         const double actualTime) {
    const bool isLocal = upperLink.erase(son) == 1;
    if(isLocal) {
      son->upperLink.insert(this);
    }

    const int nbSons = MT_sons.size();
    const int nbGrandsons = son->MT_sons.size();

    if(nbSons == 1 && nbGrandsons == 1 && !isLocal) {
      fastSwap(this, son, swapQueue, actualTime);
      return;
    }

    MergeTreeLinkCutNode *const myParent = this->MT_parent;
    cut(this, actualTime);
    cut(son, actualTime);

    // We need a temporary copy because iterators will be invalidated when we
    // cut.
    std::vector<MergeTreeLinkCutNode *> listOfGrandsons(
      son->MT_sons.begin(), son->MT_sons.end());
    for(auto grandson : listOfGrandsons)
      cut(grandson, actualTime);

    std::set<MergeTreeLinkCutNode *> seenRoots;
    for(auto neigh : upperLink)
      seenRoots.insert(root(neigh, actualTime));

    for(auto rt : seenRoots) {
      auto idx = std::find(listOfGrandsons.begin(), listOfGrandsons.end(), rt);
      if(idx != listOfGrandsons.end()) {
        newLink(*idx, this, swapQueue, actualTime);
        listOfGrandsons.erase(idx);
      }
    }
    for(auto unseen : listOfGrandsons)
      link(unseen, son, actualTime);

    newLink(son, myParent, swapQueue, actualTime);
    link(this, son, actualTime); // not a newLink because they've just swap

    const int nbGrandsonKept = son->MT_sons.size() - 1;

    // treat simple cases for PD following
    if(nbGrandsons == 0) {
      // Maximum transfer or destruction
      if(nbSons == 1) {
        this->pairOfMax = son->pairOfMax;
        this->pairOfMax->max = this;
        updatePairEventTime(this->pairOfMax, swapQueue, actualTime);
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
        = swapQueue
            .insert({2., NULL, NULL, true, this->pairOfMax
#ifndef NODEBUG
                     ,
                     actualTime, son->PT_max
#endif
            })
            .first;

      updatePairEventTime(this->pairOfMax, swapQueue, actualTime);
      return;
    }

    // If we haven't returned yet, we need to rebuild saddle-maxima links
    // Here, all of *this and *son 's MT children are in a different splay
    // trees, so their actualMax are right !

    std::vector<MergeTreeLinkCutNode *> losingPairs;
    std::vector<NodePair *> toUpdateTime;
    MergeTreeLinkCutNode *winningPair = NULL;
    for(auto branch : this->PT_sons) {
      if(winningPair == NULL) // Exactly the first iteration
        winningPair = branch;
      else {
        losingPairs.push_back(branch);
        if(getVal(branch->actualMax, actualTime)
           > getVal(winningPair->actualMax, actualTime))
          std::swap(winningPair, losingPairs.back());
      }
    }
    for(auto branch : losingPairs) {
      branch->actualMax->pairOfMax->saddle = this;
      toUpdateTime.push_back(branch->actualMax->pairOfMax);
    }

    losingPairs.clear();

    for(auto branch : son->PT_sons) {
      losingPairs.push_back(branch);
      if(getVal(branch->actualMax, actualTime)
         > getVal(winningPair->actualMax, actualTime))
        std::swap(winningPair, losingPairs.back());
    }
    for(auto branch : losingPairs) {
      branch->actualMax->pairOfMax->saddle = son;
      toUpdateTime.push_back(branch->actualMax->pairOfMax);
    }

    if(winningPair->actualMax->pairOfMax->saddle == this
       || winningPair->actualMax->pairOfMax->saddle == son)
      std::cerr << "found a saddle where noone wins" << endl;

    for(auto nodepair : toUpdateTime)
      updatePairEventTime(nodepair, swapQueue, actualTime);

    return;
  }

  std::pair<std::vector<MergeTreeLinkCutNode>, EventQueue>
    buildTree(const AbstractTriangulation *grid,
              double *scalarsStart,
              double *scalarsEnd) {
    std::pair<std::vector<MergeTreeLinkCutNode>, EventQueue> returnVals;
    EventQueue &swapQueue = returnVals.second;
    const double actualTime = 0.;
    SimplexId nbNodes = grid->getNumberOfVertices();

    std::vector<SimplexId> insertOrder(nbNodes);
    std::vector<MergeTreeLinkCutNode> &treeData = returnVals.first;
    treeData.resize(nbNodes + 2);
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
          auto neighRoot = root(&treeData[neigh], actualTime);
          if(neighRoot != that) {
#ifndef NODEBUG
            std::cout << "it is linked" << endl;
#endif

            newLink(neighRoot, that, swapQueue, actualTime);
          }
        }
      }
      if(that->MT_sons.size() == 0) {
        that->pairOfMax = new NodePair{that, globRoot, idNode, swapQueue.end()};
      } else if(that->MT_sons.size() > 1) {
        access(that, actualTime);
        auto maxSon = that->PT_max;
        for(auto son : that->MT_sons) {
          access(son, actualTime);
          if(son->PT_max != maxSon) {
            auto np = son->PT_max->pairOfMax;
            np->saddle = that;
            np->itToEvent = swapQueue
                              .insert({crossTime(son->PT_max, that->PT_max),
                                       NULL, NULL, true, np
#ifndef NODEBUG
                                       ,
                                       actualTime, that->PT_max
#endif
                              })
                              .first;
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

    link(globTop, globRoot, actualTime);

    auto oldRoot = &treeData[insertOrder[nbNodes - 1]];
    link(oldRoot, globRoot, actualTime);
    access(oldRoot, actualTime);
    auto np = oldRoot->PT_max->pairOfMax;
    np->saddle = globRoot;
    np->itToEvent = swapQueue
                      .insert({2., NULL, NULL, true, np
#ifndef NODEBUG
                               ,
                               actualTime, NULL
#endif
                      })
                      .first;

    return returnVals;
  }

  // maybe this will be over some value yet to be inserted...
  // Not for maxSwap anymore, but maybe treeSwap
  void
    setActuTime(double eventTime, EventQueue &swapQueue, double &actualTime) {
    const double oldActualTime[[maybe_unused]] = actualTime;

    actualTime
      = eventTime + (swapQueue.begin()->timestamp - eventTime) / 16777216.;
    if(actualTime == eventTime) {
      actualTime
        = eventTime + (swapQueue.begin()->timestamp - eventTime) / 2048.;
      if(actualTime == eventTime) {
        actualTime
          = eventTime + (swapQueue.begin()->timestamp - eventTime) / 2.;
        if(actualTime == eventTime) {
          std::cerr << "ALERT : possible numerical instability at time "
                    << actualTime << std::endl;
          // TODO +eps is a anticipated hotfix ; to be done better someday
          actualTime += 0.0000000001;
        }
      }
    }
#ifndef NOTIMEDISPLAY
    if(ceil(actualTime * 1000) != ceil(oldActualTime * 1000))
      std::cout << "Time : " << actualTime << endl;
#endif
  }

  void loopQueue(EventQueue &swapQueue) {
    double actualTime = 0.;
#ifndef NODEBUG
    std::cout << "queue size " << swapQueue.size() << " starting at "
              << swapQueue.begin()->timestamp << endl;
#endif
    swapQueue.insert({1., NULL, NULL, false, NULL});

    cerr << std::setprecision(15);
    // TODO switch back to 1.
    while(swapQueue.begin()->timestamp < 1.) {
      auto eventIt = swapQueue.begin();
      auto event = *eventIt;
      double eventTime = event.timestamp;
#ifndef NODEBUG
      if(event.isMaxSwap)
        std::cout << "queue size " << swapQueue.size() << ", event "
                  << (event.isMaxSwap ? "maxSwap" : "treeSwap") << " with time "
                  << eventTime << endl;
#endif

      if(event.isMaxSwap) {
        MergeTreeLinkCutNode *theSaddle = event.pairSwap->saddle;
        access(theSaddle, actualTime);

#ifndef NODEBUG
        if(theSaddle->MT_sons.size() > 2) {
          cerr << "Multiple saddle pair swap !" << endl;
          double lol = crossTime(theSaddle->PT_max, event.pairSwap->max);
          cerr << "Has " << theSaddle->MT_sons.size() << " sons, swap time "
               << lol << endl;
        }
#endif

        MergeTreeLinkCutNode *theOldMax = theSaddle->PT_max;

        if(updatePairEventTime(event.pairSwap, swapQueue, actualTime))
          continue;

#ifndef NODEBUG
        std::cout << "accomplish it" << endl;
#endif
        MergeTreeLinkCutNode *newMax = event.pairSwap->max;

        MergeTreeLinkCutNode *newSad = theOldMax->pairOfMax->saddle;
        access(newSad, actualTime);
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
        updatePairEventTime(newMax->pairOfMax, swapQueue, actualTime);
        setActuTime(eventTime, swapQueue, actualTime);
        access(theSaddle, actualTime);
        access(newSad, actualTime);
        swapQueue.erase(theOldMax->pairOfMax->itToEvent);
        theOldMax->pairOfMax->itToEvent
          = swapQueue.insert({2., NULL, NULL, true, theOldMax->pairOfMax})
              .first;
      } else {
        swapQueue.erase(eventIt);
        if(event.leafing->MT_sons.count(event.rooting) == 0)
          continue; // This event is outdated
        setActuTime(eventTime, swapQueue, actualTime);
        event.leafing->swapWithSon(event.rooting, swapQueue, actualTime);
      }
      if(actualTime >= swapQueue.begin()->timestamp) {
        std::cerr << "ALERT : possible time inconsistency, inserted event "
                     "before actual time."
                  << endl;
#ifndef NODEBUG
        auto nextEv = *swapQueue.begin();
        if(nextEv.isMaxSwap)
          std::cerr << "Next event (max) :" << nextEv.pairSwap->max << " ^, "
                    << nextEv.pairSwap->saddle << " sad ; at time "
                    << nextEv.timestamp << endl;
        else
          std::cerr << "Next event (tree) :" << nextEv.leafing->numForDebug
                    << " ^, " << nextEv.rooting->numForDebug << " v ; at time "
                    << nextEv.timestamp << endl;

        if(event.isMaxSwap)
          std::cerr << "Old event (max) :" << event.pairSwap->max << " ^, "
                    << event.pairSwap->saddle << " sad ; at time " << eventTime
                    << endl;
        else
          std::cerr << "Old event (tree) :" << event.leafing->numForDebug
                    << " ^, " << event.rooting->numForDebug << " v ; at time "
                    << eventTime << endl;
        if(eventTime == nextEv.timestamp)
          cerr << "Time is the same" << endl;
#endif
      }
    }
  }

} // namespace ttk
