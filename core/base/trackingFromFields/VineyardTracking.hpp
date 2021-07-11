/* Link-cut structure is a courtesy of
 * https://github.com/saadtaame/link-cut-tree
 */

#pragma once

#include "AbstractTriangulation.h"
#include "DataTypes.h"
#include "FTMStructures.h"
#include "FTMTree_MT.h"
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

  typedef std::multimap<double, SwapEvent> EventQueue;
  /*
  void update(MergeTreeLinkCutNode * const  x);
  void rotr(MergeTreeLinkCutNode * const x);
  void rotl(MergeTreeLinkCutNode * const x);
  void splay(MergeTreeLinkCutNode * const x);
  MergeTreeLinkCutNode *access(MergeTreeLinkCutNode * const x);
  MergeTreeLinkCutNode *root(MergeTreeLinkCutNode * const x);
  void cut(MergeTreeLinkCutNode * const x);
  void link(MergeTreeLinkCutNode * const x, MergeTreeLinkCutNode * const y);
  MergeTreeLinkCutNode *lca(MergeTreeLinkCutNode * const x, MergeTreeLinkCutNode
  * const y);

  void newLink(MergeTreeLinkCutNode* const son, MergeTreeLinkCutNode* const
  parent); // */

  struct NodePair {
    MergeTreeLinkCutNode *max, *saddle;
    SimplexId idFirstMax;

    EventQueue::iterator itToEvent;
  };

  struct MergeTreeLinkCutNode {
#ifndef NODEBUG
    SimplexId numForDebug;
#endif
    /* Splay-tree */
    MergeTreeLinkCutNode *ST_parent;
    MergeTreeLinkCutNode *ST_left; // upper in path
    MergeTreeLinkCutNode *ST_right; // lower in path

    /* Link between splay trees */
    MergeTreeLinkCutNode
      *PT_parent; // if ST root then real parent of path start, else NULL
    std::set<MergeTreeLinkCutNode *> PT_sons; // inverse of PT_parent

    /* Explicit (merge) tree */
    MergeTreeLinkCutNode *MT_parent; // real parent
    std::set<MergeTreeLinkCutNode *> MT_sons;

    /* Field geometric structure */
    std::set<MergeTreeLinkCutNode *> upperLink;
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

    /*
    MergeTreeLinkCutNode(double _scalarStart, double _scalarEnd,
    std::set<MergeTreeLinkCutNode*>& _upperLink) : ST_parent(NULL),
    ST_left(NULL), ST_right(NULL), PT_parent(NULL), MT_parent(NULL), MT_sons(),
    upperLink(_upperLink), scalarStart(_scalarStart), scalarEnd(_scalarEnd) {
      }
    void addSon(MergeTreeLinkCutNode* son) {
      MT_sons.insert(son);
      ttk::link(son, this);
    } // */
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
    MergeTreeLinkCutNode *rooting;
    MergeTreeLinkCutNode *leafing;
    bool isMaxSwap;
    NodePair *pairSwap;

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
