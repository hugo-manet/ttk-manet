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
#define myAssert(x) //funcAssert(x, #x , __LINE__) 
  void funcAssert(bool val, const char* str, int numLigne) {
    if (!val) {
      std::cerr << "ca plante ! " << str << " en ligne " << numLigne << endl;
      std::cerr << ((MergeTreeLinkCutNode *)NULL)->actualMax->MT_sons.size();
    }
  }

  void doEraseOne(std::vector<MergeTreeLinkCutNode *> &theVec,
                  MergeTreeLinkCutNode *val) {
    switch(theVec.size()) {
      case 1:
        theVec.clear();
        return;
      case 2:
        if(theVec[0] == val)
          theVec[0] = theVec[1];
        theVec.resize(1);
        return;
      default:
        theVec.erase(std::find(theVec.begin(), theVec.end(), val));
    }
  }
  void doReplace(std::vector<MergeTreeLinkCutNode *> &theVec,
                 MergeTreeLinkCutNode *val,
                 MergeTreeLinkCutNode *newVal) {
    switch(theVec.size()) {
      case 1:
        theVec[0] = newVal;
        return;
      case 2:
        if(theVec[0] == val)
          theVec[0] = theVec[1];
        theVec[1] = newVal;
        return;
      default:
        *std::find(theVec.begin(), theVec.end(), val) = newVal;
    }
  }
  bool doHas(std::vector<MergeTreeLinkCutNode *> &theVec,
             MergeTreeLinkCutNode *val) {
    switch(theVec.size()) {
      case 0:
        return false;
      case 1:
        return theVec[0] == val;
      case 2:
        return theVec[0] == val || theVec[1] == val; /*
      case 3:
        return theVec[0] == val || theVec[1] == val || theVec[2] == val;
      case 4:
        return theVec[0] == val || theVec[1] == val || theVec[2] == val ||
      theVec[3] == val; // */
      default:
        return std::find(theVec.begin(), theVec.end(), val) != theVec.end();
    }
  }

  constexpr double getVal(const MergeTreeLinkCutNode *const node,
                          const double actualTime) {
    return node->scalarStart
           + actualTime * (node->scalarEnd - node->scalarStart);
  }
  constexpr double crossTimeSwapped(const MergeTreeLinkCutNode *const son,
                                    const MergeTreeLinkCutNode *const parent) {
#define DENOM \
  (parent->scalarEnd - parent->scalarStart - son->scalarEnd + son->scalarStart)
#define REALCROSS ((son->scalarStart - parent->scalarStart) / DENOM)
    return (DENOM != 0.)
             ? ((0. <= REALCROSS && REALCROSS <= 2.) ? REALCROSS : 2.)
             : 2.;
  }
  constexpr double crossTime(const MergeTreeLinkCutNode *const son,
                             const MergeTreeLinkCutNode *const parent) {
#define DENOM \
  (parent->scalarEnd - parent->scalarStart - son->scalarEnd + son->scalarStart)
#define REALCROSS ((son->scalarStart - parent->scalarStart) / DENOM)
    return (son->scalarStart > parent->scalarStart)
             ? crossTimeSwapped(son, parent)
             : crossTimeSwapped(parent, son);
  }
  double crossTimeNC(const MergeTreeLinkCutNode *const son,
                     const MergeTreeLinkCutNode *const parent) {

    return crossTime(son, parent);
  }

  void setMaxFinal(MergeTreeLinkCutNode *&dest,
                   const MergeTreeLinkCutNode *const option,
                   const double actualTime) {
    const double dVal = getVal(dest, actualTime),
                 oVal = getVal(option->actualMax, actualTime);

    if(dVal < oVal)
      dest = option->actualMax;
  }
  void setMax(MergeTreeLinkCutNode *&dest,
              const MergeTreeLinkCutNode *const option,
              const double actualTime) {
    if(option == NULL)
      return;
    double dVal = getVal(dest, actualTime),
           oVal = getVal(option->actualMax, actualTime);
    if(dVal - oVal < 1.e-10 && oVal - dVal < 1.e-10) {
      double cT = crossTime(dest, option->actualMax);
      if(cT != 2) {
        if(cT < actualTime) // <= ?
          setMaxFinal(dest, option, 1.);
        else
          setMaxFinal(dest, option, 0.);
        return;
      }
    }
    if(dVal < oVal)
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

  // /!\ does not conserve ST_first and PT_sons/PT_parent ! That is done in
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
        // Move ST_first to the new root
        myAssert(y->ST_first != NULL);
        x->ST_first = y->ST_first;
        y->ST_first = NULL;
        if(y->PT_parent != NULL) {
          x->PT_parent = y->PT_parent;
          y->PT_parent = 0;
          doReplace(x->PT_parent->PT_sons, y, x);
        }

        if(x == y->ST_left)
          rotr(x, actualTime);
        else
          rotl(x, actualTime);
      } else {
        MergeTreeLinkCutNode *const z = y->ST_parent;
        if(z->ST_parent == NULL) {
          // Move ST_first to the new root
          myAssert(z->ST_first != NULL);
          x->ST_first = z->ST_first;
          z->ST_first = NULL;
          if(z->PT_parent != NULL) {
            x->PT_parent = z->PT_parent;
            z->PT_parent = 0;
            doReplace(x->PT_parent->PT_sons, z, x);
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
      // cut linked list
      myAssert(x->ST_next != NULL);
      x->ST_right->ST_first = x->ST_next;
      x->ST_next = NULL;

      // cut tree
      x->ST_right->PT_parent = x;
      x->PT_sons.emplace_back(x->ST_right);
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
        myAssert(y->ST_next != NULL);
        y->ST_right->ST_first = y->ST_next;
        y->ST_next = NULL;

        y->ST_right->PT_parent = y;
        y->PT_sons.emplace_back(y->ST_right);
        y->ST_right->ST_parent = 0;
      }
      y->ST_right = x;
      myAssert(x->ST_first != NULL);
      y->ST_next = x->ST_first;
      x->ST_first = NULL;
      x->ST_parent = y;
      doEraseOne(y->PT_sons, x);
      x->PT_parent = 0;
      update(y, actualTime);
      splay(x, actualTime);
    }

    return last;
  }

  MergeTreeLinkCutNode *root(MergeTreeLinkCutNode *x, const double actualTime) {
    access(x, actualTime);
    
    return x->ST_first;
  }

  void cut(MergeTreeLinkCutNode *const x, const double actualTime) {
    access(x, actualTime);
    myAssert(x->ST_left != NULL); // because x had a parent.
    myAssert(x->MT_parent != NULL); // because x had a parent.
    myAssert(x->ST_first != NULL);
    x->ST_left->ST_first = x->ST_first;
    x->ST_first = x;
    myAssert(x->MT_parent->ST_next == x);
    x->MT_parent->ST_next = NULL;
      
    x->ST_left->ST_parent = 0;
    x->ST_left = 0;

    doEraseOne(x->MT_parent->MT_sons, x);
    x->MT_parent = NULL;

    update(x, actualTime);
  }

  void link(MergeTreeLinkCutNode *const x,
            MergeTreeLinkCutNode *const y,
            const double actualTime) {
    access(x, actualTime);
    access(y, actualTime);
    myAssert(x->ST_left == NULL); // because x had no parent.
    x->MT_parent = y;
    y->MT_sons.emplace_back(x);

    x->ST_left = y;
    myAssert(y->ST_next == 0);
    y->ST_next = x;
    myAssert(y->ST_first != 0);
    x->ST_first = y->ST_first;
    y->ST_first = NULL;
    y->ST_parent = x;
    update(x, actualTime);
  }

  MergeTreeLinkCutNode *lca(MergeTreeLinkCutNode *const x,
                            MergeTreeLinkCutNode *const y,
                            const double actualTime) {
    access(x, actualTime);
    return access(y, actualTime);
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
      swapQueue.push({crossing, son, parent, false, NULL
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

#ifndef NODEBUG
    if(ofThePair->saddle->PT_max == ofThePair->max) {
      auto val = &(cerr << "A pair is gonna try to swap with itself !");
      cerr << val->bad();
      val = &((*val) << endl);
      throw "MDR";
    }
#endif

    double newTime = crossTime(ofThePair->saddle->PT_max, ofThePair->max);
    if(ofThePair->saddle->PT_max->scalarEnd >= ofThePair->max->scalarEnd)
      newTime = 2.;

    if(newTime != (*ofThePair->itToEvent).timestamp) {
#ifndef NODEBUG
      if(newTime <= actualTime) {
        auto val = &(cerr << "Updated a time into a time earlier than now !");
        cerr << val->bad();
        val = &((*val) << endl);
      }
#endif

      // Grrrrr. It's a whole copy for nothing :'(
      auto &copy = ofThePair->itToEvent;
      (*copy).timestamp = newTime;
#ifndef NODEBUG
      (*copy).timeOfLastUpdate = actualTime;
      (*copy).wasLosingTo = ofThePair->saddle->PT_max;
#endif
      swapQueue.update(copy);
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
      doReplace(node->PT_parent->PT_sons, other, node);
    }

    if(!node->PT_sons.empty()) {
      if(*node->PT_sons.begin() == node) {
        node->PT_sons.clear();
        node->PT_sons.emplace_back(other);
      } else
        (*node->PT_sons.begin())->PT_parent = node;
    }

    if(node->MT_parent == node)
      node->MT_parent = other;
    else if(node->MT_parent) {
      doReplace(node->MT_parent->MT_sons, other, node);
    }

    if(!node->MT_sons.empty()) {
      if(*node->MT_sons.begin() == node) {
        node->MT_sons.clear();
        node->MT_sons.emplace_back(other);
      } else
        (*node->MT_sons.begin())->MT_parent = node;
    }

    /*
    // Now for the hard ones
    if(node->MT_parent->ST_next == other)
      node->MT_parent->ST_next = node;
    else
      for (auto implicitSon : node->MT_parent->PT_sons)
        if (implicitSon->ST_first == other) {
          implicitSon->ST_first = node;
          return;
        } // */

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
    std::swap(that->ST_next, son->ST_next);
    std::swap(that->ST_first, son->ST_first);
    std::swap(that->PT_parent, son->PT_parent);
    std::swap(that->PT_sons, son->PT_sons);
    std::swap(that->MT_parent, son->MT_parent);
    std::swap(that->MT_sons, son->MT_sons);
    std::swap(that->actualMax, son->actualMax);
    std::swap(that->PT_max, son->PT_max);

    repairNodeLinks(that, son);
    repairNodeLinks(son, that);
    //*
    if(son->MT_parent->ST_next == that)
      son->MT_parent->ST_next = son;
    else
      for (auto implicitSon : son->MT_parent->PT_sons)
        if (implicitSon->ST_first == that)
          implicitSon->ST_first = son;
    if(son->ST_next == son)
      son->ST_next = that;
    else
      for (auto implicitSon : son->PT_sons)
        if (implicitSon->ST_first == son)
          implicitSon->ST_first = that; // */
    update(son, actualTime);
    update(that, actualTime);
    update(son, actualTime); // Because son can depend on that

    createTreeSwapEvent(son, son->MT_parent, swapQueue, actualTime);
    createTreeSwapEvent(*that->MT_sons.begin(), that, swapQueue, actualTime);
  }
  void MergeTreeLinkCutNode::swapWithSon(MergeTreeLinkCutNode *const son,
                                         EventQueue &swapQueue,
                                         const double actualTime) {
    const bool isLocal = doHas(upperLink, son);
    if(isLocal) {
      doEraseOne(upperLink, son);
      son->upperLink.emplace_back(this);
    }

    const int nbSons = MT_sons.size();
    const int nbGrandsons = son->MT_sons.size();

    if(nbSons == 1 && nbGrandsons == 1 && !isLocal) {
      fastSwap(this, son, swapQueue, actualTime);
      return;
    }

    MergeTreeLinkCutNode *const myParent = this->MT_parent;
    cut(this, actualTime);

#ifndef NOCHECKSADDLES
    if(nbSons >= 2) {
      int nbDyingMe = 0;
      for(auto mySonImplicit : this->PT_sons)
        if(mySonImplicit->actualMax->pairOfMax->saddle == this)
          ++nbDyingMe;

      if(nbDyingMe != this->PT_sons.size() - 1)
        std::cerr << "Invariant invalid for me : " <<
#ifndef NODEBUG
          this->numForDebug
#else
          this
#endif
                  << " because " << nbDyingMe
                  << " die instead of nbGrandson - 1 = "
                  << (this->PT_sons.size() - 1) << endl;
    }
#endif

    cut(son, actualTime);

    // We need a temporary copy because iterators will be invalidated when we
    // cut.
    std::vector<MergeTreeLinkCutNode *> listOfGrandsons(
      son->MT_sons.begin(), son->MT_sons.end());
    for(auto grandson : listOfGrandsons)
      cut(grandson, actualTime);

#ifndef NOCHECKSADDLES
    if(nbGrandsons >= 2) {
      int nbDying = 0;
      for(auto grandson : listOfGrandsons)
        if(grandson->actualMax->pairOfMax->saddle == son)
          ++nbDying;

      if(nbDying != listOfGrandsons.size() - 1)
        std::cerr << "Invariant invalid for son : " <<
#ifndef NODEBUG
          son->numForDebug
#else
          son
#endif
                  << " because " << nbDying
                  << " die instead of nbGrandson - 1 = "
                  << (listOfGrandsons.size() - 1) << endl;
    }
#endif

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
      this->pairOfMax = new NodePair{this, son, -1};
      this->pairOfMax->itToEvent
        = swapQueue.push({2., NULL, NULL, true, this->pairOfMax
#ifndef NODEBUG
                          ,
                          actualTime, son->PT_max
#endif
        });

      updatePairEventTime(this->pairOfMax, swapQueue, actualTime);
      return;
    }

    // If we haven't returned yet, we need to rebuild saddle-maxima links
    // Here, all of *this and *son 's MT children are in a different splay
    // trees, so their actualMax are right !

    std::vector<MergeTreeLinkCutNode *> losingPairs;
    std::vector<NodePair *> toUpdateTime;
    MergeTreeLinkCutNode *winningPair = NULL;
    MergeTreeLinkCutNode *theRealMax = NULL;
    for(auto branch : this->PT_sons) {
      if(winningPair == NULL) { // Exactly the first iteration
        winningPair = branch;
        theRealMax = branch->actualMax;
      } else {
        losingPairs.emplace_back(branch);
        setMax(theRealMax, branch, actualTime);
        if(theRealMax == branch->actualMax)
          std::swap(winningPair, losingPairs.back());
      }
    }
    MergeTreeLinkCutNode *oldFarSaddle = NULL;
    for(auto branch : losingPairs) {
      if(branch->actualMax->pairOfMax->saddle != this
         && branch->actualMax->pairOfMax->saddle != son) {
        std::cerr << "found a saddle where we erase a saddle" << endl;
        oldFarSaddle = branch->actualMax->pairOfMax->saddle;
      }
      branch->actualMax->pairOfMax->saddle = this;
      toUpdateTime.emplace_back(branch->actualMax->pairOfMax);
    }

    losingPairs.clear();

    for(auto branch : son->PT_sons) {
      losingPairs.emplace_back(branch);
      setMax(theRealMax, branch, actualTime);
      if(theRealMax == branch->actualMax)
        std::swap(winningPair, losingPairs.back());
    }
    for(auto branch : losingPairs) {
      if(branch->actualMax->pairOfMax->saddle != this
         && branch->actualMax->pairOfMax->saddle != son) {
        std::cerr << "found a saddle where we erase a saddle, definitely"
                  << endl;
        oldFarSaddle = branch->actualMax->pairOfMax->saddle;
      }
      branch->actualMax->pairOfMax->saddle = son;
      toUpdateTime.emplace_back(branch->actualMax->pairOfMax);
    }

    if(winningPair->actualMax->pairOfMax->saddle == this
       || winningPair->actualMax->pairOfMax->saddle == son) {
      std::cerr << "found a saddle where noone pretends to win" << endl;
      if(oldFarSaddle) {
        std::cerr << "Found a candidate to correct that" << endl;
        winningPair->actualMax->pairOfMax->saddle = oldFarSaddle;
        toUpdateTime.emplace_back(winningPair->actualMax->pairOfMax);
      }
    }

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
#ifdef DEBUGTREE
      std::cout << "inserting " << idNode << endl;
#endif
#ifndef NODEBUG
      that->numForDebug = idNode;
#endif
      that->scalarStart = scalarsStart[idNode];
      that->scalarEnd = scalarsEnd[idNode];
      that->ST_first = that;
      that->actualMax = that->PT_max = that;
      SimplexId neigh;
      for(int iNeigh = 0; iNeigh < grid->getVertexNeighborNumber(idNode);
          ++iNeigh) {
        grid->getVertexNeighbor(idNode, iNeigh, neigh);
        if(compareIdx(idNode, neigh)) {

#ifdef DEBUGTREE
          std::cout << "iNeigh " << iNeigh
                    << " is to be inserted (and maybe linked)" << endl;
#endif
          that->upperLink.emplace_back(&treeData[neigh]);
          auto neighRoot = root(&treeData[neigh], actualTime);
          if(neighRoot != that) {
#ifdef DEBUGTREE
            std::cout << "it is linked" << endl;
#endif

            newLink(neighRoot, that, swapQueue, actualTime);
          }
        }
      }
      if(that->MT_sons.size() == 0) {
        that->pairOfMax = new NodePair{that, globRoot, idNode};
      } else if(that->MT_sons.size() > 1) {
        access(that, actualTime);
        auto maxSon = that->PT_max;
        for(auto son : that->MT_sons) {
          access(son, actualTime);
          if(son->PT_max != maxSon) {
            auto np = son->PT_max->pairOfMax;
            np->saddle = that;
            np->itToEvent = swapQueue.push(
              {crossTime(that->PT_max, son->PT_max), NULL, NULL, true, np
#ifndef NODEBUG
               ,
               actualTime, that->PT_max
#endif
              });
          }
        }
      }
    }
    // Build sentinels
    globTop->scalarStart = globTop->scalarEnd = globMax * 1.01;
    globTop->ST_first = globTop;
    globTop->actualMax = globTop->PT_max = globTop;
    globRoot->scalarStart = globRoot->scalarEnd
      = globMin * 1.01; // globMin < -1
    globRoot->ST_first = globRoot;
    globRoot->actualMax = globRoot->PT_max = globRoot;

    link(globTop, globRoot, actualTime);

    auto oldRoot = &treeData[insertOrder[nbNodes - 1]];
    link(oldRoot, globRoot, actualTime);
    access(oldRoot, actualTime);
    auto np = oldRoot->PT_max->pairOfMax;
    np->saddle = globRoot;
    np->itToEvent = swapQueue.push({2., NULL, NULL, true, np
#ifndef NODEBUG
                                    ,
                                    actualTime, NULL
#endif
    });

    return returnVals;
  }

  // maybe this will be over some value yet to be inserted...
  // Not for maxSwap anymore, but maybe treeSwap
  void
    setActuTime(double eventTime, EventQueue &swapQueue, double &actualTime) {
    const double oldActualTime[[maybe_unused]] = actualTime;

    actualTime
      = eventTime + (swapQueue.top().timestamp - eventTime) / 16777216.;
    if(actualTime == eventTime) {
      if(eventTime == swapQueue.top().timestamp) {
        auto nextEvent = swapQueue.top();
        if(!nextEvent.isMaxSwap) {
          // Probable situation : we had a maxSwap, & an old treeSwap event has
          // the same nodes. We check whether the next event(s) can be deleted,
          // remove those, and restart.
          if(!doHas(nextEvent.leafing->MT_sons, nextEvent.rooting)) {
            while(!(swapQueue.top() > nextEvent))
              swapQueue.pop();
            nextEvent.leafing->eraseRecentSons(nextEvent.rooting);
            setActuTime(eventTime, swapQueue, actualTime);
            return;
          }
        }
      }
      actualTime = eventTime + (swapQueue.top().timestamp - eventTime) / 2048.;
      if(actualTime == eventTime) {
        actualTime = eventTime + (swapQueue.top().timestamp - eventTime) / 2.;
        if(actualTime == eventTime) {
          std::cerr << "ALERT : possible numerical instability at time "
                    << actualTime << std::endl;
#ifndef NODEBUG
          auto nextEv = swapQueue.top();
          if(nextEv.isMaxSwap)
            std::cerr << "Next event (max) :"
                      << nextEv.pairSwap->max->numForDebug << " ^, "
                      << nextEv.pairSwap->saddle->numForDebug
                      << " sad ; at time " << nextEv.timestamp << endl;
          else
            std::cerr << "Next event (tree) :" << nextEv.leafing->numForDebug
                      << " ^, " << nextEv.rooting->numForDebug
                      << " v ; at time " << nextEv.timestamp << endl;

          if(eventTime == nextEv.timestamp)
            cerr << "Time is the same" << endl;
#endif
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
    double lastRunEvent = 0.;
    std::cout << "queue size " << swapQueue.size() << " starting at "
              << swapQueue.top().timestamp << endl;
    swapQueue.push({1., NULL, NULL, false, NULL});

    cerr << std::setprecision(15);
    while(swapQueue.top().timestamp < 1.) {
      SwapEvent event = swapQueue.top();
      double eventTime = event.timestamp;
#ifndef NODEBUG
#ifdef PRINTALLMAXSWAPS
      if(event.isMaxSwap)
        std::cout << "queue size " << swapQueue.size() << ", event "
                  << (event.isMaxSwap ? "maxSwap" : "treeSwap") << " with time "
                  << eventTime << endl;
#endif
#endif

      if(event.isMaxSwap) {
        if(actualTime >= eventTime)
          setActuTime(
            lastRunEvent, swapQueue, actualTime); // Trying to hotfix...
        MergeTreeLinkCutNode *theSaddle = event.pairSwap->saddle;
        access(theSaddle, actualTime);

#ifndef NODEBUG
#ifdef PRINTALLMAXSWAPS
        if(theSaddle->MT_sons.size() > 2) {
          cerr << "Multiple saddle pair swap !" << endl;
          double lol = crossTime(theSaddle->PT_max, event.pairSwap->max);
          cerr << "Has " << theSaddle->MT_sons.size() << " sons, swap time "
               << lol << endl;
        }
#endif
#endif

        MergeTreeLinkCutNode *theOldMax = theSaddle->PT_max;

        if(updatePairEventTime(event.pairSwap, swapQueue, actualTime)) {
#ifndef NODEBUG
#ifdef PRINTALLMAXSWAPS
          double newTime
            = (*event.pairSwap->max->pairOfMax->itToEvent).timestamp;
          std::cerr << "Changed time of " << event.pairSwap->max->numForDebug
                    << "^ ; sad " << event.pairSwap->saddle->numForDebug
                    << " from " << event.timestamp << " to " << newTime
                    << " because we prefer" << endl;

          if(event.timestamp == newTime)
            std::cerr << "Time are the same???" << endl;
#endif
#endif
          continue;
        }
#ifndef NODEBUG
#ifdef PRINTALLMAXSWAPS
        std::cout << "accomplish it" << endl;
#endif
#endif
        MergeTreeLinkCutNode *newMax = event.pairSwap->max;

        MergeTreeLinkCutNode *newSad = theOldMax->pairOfMax->saddle;
        access(newSad, actualTime);
#ifndef NODEBUG
#ifdef PRINTALLMAXSWAPS
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
#endif
        // Maybe swap values and not ptrs ?
        std::swap(theOldMax->pairOfMax->saddle, newMax->pairOfMax->saddle);
        std::swap(
          theOldMax->pairOfMax->idFirstMax, newMax->pairOfMax->idFirstMax);
        (*theOldMax->pairOfMax->itToEvent) = {2.,
                                              NULL,
                                              NULL,
                                              true,
                                              theOldMax->pairOfMax
#ifndef NODEBUG
                                              ,
                                              actualTime,
                                              newMax
#endif
        };
        swapQueue.update(theOldMax->pairOfMax->itToEvent);
        (*newMax->pairOfMax->itToEvent) = {2.,
                                           NULL,
                                           NULL,
                                           true,
                                           newMax->pairOfMax
#ifndef NODEBUG
                                           ,
                                           actualTime,
                                           NULL
#endif
        };
        swapQueue.update(newMax->pairOfMax->itToEvent);
        setActuTime(eventTime, swapQueue, actualTime);
        update(theSaddle, actualTime);
        update(newSad, actualTime);
        updatePairEventTime(newMax->pairOfMax, swapQueue, actualTime);
        lastRunEvent = eventTime;
      } else {
        while(!(swapQueue.top() > event))
          swapQueue.pop();
        if(!doHas(event.leafing->MT_sons, event.rooting)) {
          event.leafing->eraseRecentSons(event.rooting);
          continue; // This event is outdated
        }
        setActuTime(eventTime, swapQueue, actualTime);
        event.leafing->swapWithSon(event.rooting, swapQueue, actualTime);
        lastRunEvent = eventTime;
      }
      if(actualTime >= swapQueue.top().timestamp) {
        std::cerr << "ALERT : possible time inconsistency, inserted event "
                     "before actual time."
                  << endl;
#ifndef NODEBUG
        auto nextEv = swapQueue.top();
        if(nextEv.isMaxSwap)
          std::cerr << "Next event (max) :" << nextEv.pairSwap->max->numForDebug
                    << " ^, " << nextEv.pairSwap->saddle->numForDebug
                    << " sad ; at time " << nextEv.timestamp << endl;
        else
          std::cerr << "Next event (tree) :" << nextEv.leafing->numForDebug
                    << " ^, " << nextEv.rooting->numForDebug << " v ; at time "
                    << nextEv.timestamp << endl;

        if(event.isMaxSwap)
          std::cerr << "Old event (max) :" << event.pairSwap->max->numForDebug
                    << " ^, " << event.pairSwap->saddle->numForDebug
                    << " sad ; at time " << eventTime << endl;
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
