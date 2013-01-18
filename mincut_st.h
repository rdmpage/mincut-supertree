// $Id: mincut_st.h,v 1.1.1.1 2001/10/31 10:43:03 rdmp1c Exp $

/**
 * @file mincut_st.h
 *
 * Mincut of a graph
 *
 */

#ifndef MINIMUMCUTH
#define MINIMUMCUTH

#include <GTL/graph.h>

/**
 * @var node_pair
 * A pair of nodes
 */

typedef pair<node, node> node_pair;

/**
 * @fn mincut_st (const graph &G0, edge_map <int> &w0, list<node_pair> &st_list)
 * @brief Find mincut of a graph
 * Computes a global mincut for a graph, and stores a list of pairs of (s,t)
 * nodes that have this cut value.
 * @param G0 the input graph (unmodified by this function)
 * @param w0 the edges weights
 * @param st_list a list of (s,t) node pairs for which a (s,t) cut is a best cut
 *
 * @return the value of the minimum cut, and the list of (s,t) pairs in st_list
 
 */
int mincut_st (const graph &G0, edge_map <int> &w0, list<node_pair> &st_list);

#endif
