/* Link-cut structure is a courtesy of
 * https://github.com/saadtaame/link-cut-tree
 */

#pragma once

#include "AbstractTriangulation.h"
#include "DataTypes.h"
#include "FTMStructures.h"
#include "FTMTree_MT.h"
#include <boost/heap/d_ary_heap.hpp>
#include <boost/heap/policies.hpp>
#include <map>
#include <paraview-5.8/vtkDoubleArray.h>
#include <paraview-5.8/vtkIntArray.h>
#include <paraview-5.8/vtkType.h>
#include <paraview-5.8/vtkUnstructuredGrid.h>
#include <set>

#include <Triangulation.h>

#include <PersistenceDiagram.h>

#define NODEBUG

namespace ttk {
  struct MergeTreeLinkCutNode;
  struct SwapEvent;

  typedef boost::heap::d_ary_heap<SwapEvent,
                                  boost::heap::arity<2>,
                                  boost::heap::mutable_<true>,
                                  boost::heap::compare<std::greater<SwapEvent>>>
    EventQueue;

  struct NodePair {
    MergeTreeLinkCutNode *max, *saddle;
    SimplexId idFirstMax;

    EventQueue::handle_type itToEvent;
  };

  struct MergeTreeLinkCutNode {
#ifndef NODEBUG
    SimplexId numForDebug;
#endif
    /* Splay-tree */
    MergeTreeLinkCutNode *ST_parent;
    MergeTreeLinkCutNode *ST_left; // upper in path
    MergeTreeLinkCutNode *ST_right; // lower in path

    /* Splay-tree additionnal info */
    // next in the path, it's the favorite son in MT
    MergeTreeLinkCutNode *ST_next;
    // first of the path, or heavy path's "root"/start
    MergeTreeLinkCutNode *ST_first;

    /* Link between splay trees */
    MergeTreeLinkCutNode
      *PT_parent; // if ST root then real parent of path start, else NULL
    std::vector<MergeTreeLinkCutNode *> PT_sons; // inverse of PT_parent

    /* Explicit (merge) tree */
    MergeTreeLinkCutNode *MT_parent; // real parent
    std::vector<MergeTreeLinkCutNode *> MT_sons;

    /* Field geometric structure */
    std::vector<MergeTreeLinkCutNode *> upperLink;
    double scalarStart, scalarEnd;

    MergeTreeLinkCutNode *actualMax;
    MergeTreeLinkCutNode *PT_max; // actualMax of PT_sons. Is subtree max
                                  // excluding heavy path branch.
    NodePair *pairOfMax; // Max-saddle pairs on local maxima, NULL elsewhere

#define PAST_SON_SIZE 5
    MergeTreeLinkCutNode
      *pastSons[PAST_SON_SIZE]; // if the swap event already exists
    int actuSon;
    bool hasRecentSonElseInsert(MergeTreeLinkCutNode *theSon) {
      for(int i = 0; i < PAST_SON_SIZE; ++i)
        if(pastSons[i] == theSon) {
          actuSon = i;
          return true;
        }
      pastSons[++actuSon % PAST_SON_SIZE] = theSon;
      return false;
    }
    void eraseRecentSons(MergeTreeLinkCutNode *theSon) {
      for(int i = 0; i < PAST_SON_SIZE; ++i)
        if(pastSons[i] == theSon) {
          pastSons[i] = NULL;
          return;
        }
    }

    MergeTreeLinkCutNode()
      : ST_parent(NULL), ST_left(NULL), ST_right(NULL), PT_parent(NULL),
        PT_sons(), MT_parent(NULL), MT_sons(), upperLink(), scalarStart(42.),
        scalarEnd(-42.), actualMax(NULL), PT_max(NULL), pairOfMax(NULL),
        actuSon() {
      for(int i = 0; i < PAST_SON_SIZE; ++i)
        pastSons[i] = NULL;
    }

    void swapWithSon(MergeTreeLinkCutNode *son,
                     EventQueue &swapQueue,
                     const double actualTime);
  };

  struct SwapEvent {
    double timestamp;

    MergeTreeLinkCutNode *rooting;
    MergeTreeLinkCutNode *leafing;
    bool isMaxSwap;
    NodePair *pairSwap;

    bool operator>(const SwapEvent &other) const {
      if(timestamp != other.timestamp)
        return timestamp > other.timestamp;

      if(isMaxSwap && !other.isMaxSwap)
        return false;

      if(!isMaxSwap && other.isMaxSwap)
        return true;

      if(isMaxSwap) // both
        return pairSwap > other.pairSwap;

      if(rooting != other.rooting)
        return rooting > other.rooting;
      return leafing > other.leafing;
    }
#ifndef NODEBUG
    double timeOfLastUpdate;
    MergeTreeLinkCutNode *wasLosingTo;
#endif
  };

  /* builds a tree, returns a (symbolic) root.
   *
   */
  // TODO have a -inf node
  std::pair<std::vector<MergeTreeLinkCutNode>, EventQueue>
    buildTree(const AbstractTriangulation *grid,
              double *scalarsStart,
              double *scalarsEnd);

  void loopQueue(EventQueue &swapQueue);

} // namespace ttk
