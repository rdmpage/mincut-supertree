/*
 * Supertree
 * A program for computing supertrees.
 * Copyright (C) 2001 Roderic D. M. Page <r.page@bio.gla.ac.uk>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307, USA.
 */
 
// $Id: stgraph.h,v 1.1 2002/03/19 15:57:34 rdmp1c Exp $

#ifndef STGRAPH_H
#define STGRAPH_H

// STL
#include <map>
#include <set>

// GTL
#include <GTL/graph.h>


typedef set<node, less<node> > NodeSet;


#define colour_uncontradicted 		0
#define colour_contradicted 		1
#define colour_adjto_contradicted 	2


/**
 * @class STGraph
 * STGraph extends the GTL class graph to provide support for the graphs
 * @f$S_T@f$ and @f$S_T /E_T^{\max }@f$. 
 *
 * STGraph overrides the GTL save_node_info and save_edge_info handlers to 
 * output labels for the nodes and edges when writing the graph to a GML
 * file. This simplifies debugging as the GML file can be displayed in
 * software such as Graphlet.
 *
 */
class STGraph : public graph
{
public:
	STGraph () { mShowLabels = false; mShowColours = false; };
		

	/** 
	 * A map between edges and the unweighted frequency of that edge in STEmax 
	 * @sa AddEdge
	 *
	 */
	edge_map<int> f;
	
	bool mShowLabels;
	bool mShowColours;


	/** 
	 * A map between edges and an integer weight, being the weight of that edge
	 * in the graph STEmax. 
	 * @sa AddEdge
	 *
	 */
	edge_map<int> w0;
	
		
		
	edge_map <int> edge_colour;

	/**
	 * An map between nodes and a set of nodes. By default the node set of a node n is {n},
	 * and this is set by post_new_node_handler. If one or more edges are compacted then 
	 * node sets from the adjacent nodes may be merged. 
	 * @sa post_new_node_handler()
	 */
	node_map <NodeSet> ns;
	/**
	 * A map between nodes and a label string. This allows us to access the label of a node.
	 */
	node_map <std::string> node_labels;
	/**
	 * A map between nodes and a label string. This allows us to refer to
	 * a node by its label.
	 */	
	map < std::string, node, std::less<std::string> > labelled_nodes; 

	/**
	 * Add an edge between two labelled nodes in the graph. If there is already
	 * an edge connecting the two nodes, increment the weight of that edge.
	 * @param s1 label of first node
	 * @param s2 label of second node
	 * @param weight weight of edge (default is 1)
	 */
	virtual void AddEdge (std::string s1, std::string s2, int weight = 1);
	/**
	 * Add a node labelled s to the graph, unless a node with this label already
	 * exists. 
	 * @param s the label of the node to add
	 */
	virtual void AddNode (std::string s);
	
	virtual bool EdgeExists (node n1, node n2);

	virtual int GetEdgeFreqFromNodeLabels (std::string s1, std::string s2);
	
	virtual int GetCommonNeighbours (edge e, NodeSet &common_set);
	
	/**
	 * Extends graph::post_new_node_handler to add the node n to n's node set.
	 * @param n the node just added to the graph
	 * 
	 */
	virtual void post_new_node_handler (node n);
	
	virtual void save_graph_info_handler (ostream *os) const;
	
	/**
	 * Extends graph::save_edge_info_handler to write the weight of the edge
	 * as a label when saving the graph to a GML file.
	 * @param ostream the stream being written to
	 * @param e the edge
	 * 
	 */
	virtual void save_edge_info_handler (ostream *os, edge e) const;
	/**
	 * Extends graph::save_node_info_handler to output a label
	 * for the node n when saving the graph to a GML file.
	 * @param ostream the stream being written to
	 * @param n the node
	 * 
	 */
	virtual void save_node_info_handler (ostream *os, node n) const;
	/**
      * @param s a node
      * @param t a node to merge with s
	  *
	  * Merge nodes s and t. t is added to s's node set. Any loops or 
	  * parallel edges resulting from the merge are deleted. If t is
	  * adjacent to node u then ensure there is an edge e(s,u). 
	  */
	virtual void mergeNodes (node s, node t);
	
	/**
	 * @param f output stream
	 *
	 * Write the graph in dot format (as used by programs in the
	 * <a href="http://www.research.att.com/sw/tools/graphviz/">GraphViz</A>
	 * package.
	 */
	virtual void WriteDotty (ostream &f);
	
	/**
	 * @param fname output file name
	 *
	 * Write the graph in dot format (as used by programs in the
	 * <a href="http://www.research.att.com/sw/tools/graphviz/">GraphViz</A>
	 * package.
	 */
	virtual void WriteDotty (const char *fname);

	
};	

#endif