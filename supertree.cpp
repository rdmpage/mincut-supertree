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
 
// $Id: supertree.cpp,v 1.15 2002/03/22 12:52:09 rdmp1c Exp $

/**
 * @file supertree.cpp
 *
 * Min cut supertree algorithm
 *
 */

#include "ntree.h"
#include "stree.h"
#include "profile.h"
#include "nodeiterator.h"
#include "quartet.h"


#include <GTL/graph.h>
#include <GTL/components.h>
#include <GTL/maxflow_ff.h> 
#include <GTL/biconnectivity.h>

#include <list>
#include <set>
#include <fstream>
#include <iomanip>

/*
#ifdef __GNUC__
	#include <strstream>
	# if __GNUC__ == 3
		#include <iterator>
	# endif
#else
	#include <sstream>
#endif
*/
#include <ctime>

#if __MWERKS__
	#if macintosh
		// Metrowerks support for Macintosh command line interface
		#include <console.h>
	#endif
#endif

#ifndef min
    #define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif
#ifndef max
    #define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#define MAJOR_VERSION "0"
#define MINOR_VERSION "5"
#define MINI_VERSION "0"

#include "mincut_st.h"
#include "strong_components.h"


// Modified SQUID code to handle command line options
#include "getoptions.h"
#define FILENAME_SIZE 256		// Maximum file name length

// Program options
static struct opt_s OPTIONS[] = {
	{ (char*)&"-p", true, ARG_STRING },
	{ (char*)&"-n", true, ARG_STRING },
	{ (char*)&"-k", true, ARG_STRING },
	{ (char*)&"-m", true, ARG_STRING },
	{ (char*)&"-v", true, ARG_NONE },
	{ (char*)&"-l", true, ARG_NONE },
	{ (char*)& "-b", true, ARG_NONE },
	{ (char*)&"-w", true, ARG_NONE },	
	{ (char*)&"-a", true, ARG_INT },
	{ (char*)&"-c", true, ARG_INT },		
	{ (char*)&"-d", true, ARG_NONE },
	{ (char*)&"-g", true, ARG_NONE }

};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static char usage[] = "\
Usage: supertree [-options] <tree-file> \n\
\n\
  Available options: \n\
     -v             show version information\n\
     -l             output taxon labels when writing GML files\n\
     -g             write ST and ST/Emax to GML file(s)\n\
     -n filename    write NEXUS file \n\
     -k filename    write Newick file \n\
     -p filename    write tree to Postscript file\n\
     -b             verbose\n\
     -d             write ST and ST/EMax to dot files\n\
     -w             use tree weights\n\
     -m filename    write MRP matrix to file \n\
     -a n           algorithm \n\
     -c n           compute cluster graph for k=n \n\
";


typedef set<node, less<node> > NodeSet;

// Debugging/verbose flags

// Algorithms
#define ALGORITHM_SEMPLE	0
#define ALGORITHM_ROD1		1
int use_algorithm; 

bool bSaveST 			= true;	// save ST to a GML file
bool bShowST			= true; // show ST
bool bSaveSTEmax		= true;	// save ST/Emax to a GML file
bool bShowSTEmax		= true; // show ST/Emax
bool bShowAllMinCuts 	= true; // list c(ST/e) for all e
bool bShowTrees 		= true; // show trees T1,..,Tk
bool bShowMinCutWeight	= true; // show c(ST)
bool bShowVertexSets	= true; // show the components of ST (or ST/Emax)
bool bShowRecursion		= true; // show when we call MinCutSupertree
bool bShowConnected		= true; // show whether ST is connected
bool bShowTS			= true; // show the tree after pruning leaves not in vertex set
bool bShowTStest		= true; // show tests for constructing tree from vertex set
bool bShowClusters		= true; // show clusters for input trees
bool bShowSTEdge		= true; // show edge in ST as it is created
bool bShowLevel			= true; // show level of recursion
bool bShowConstruct		= true; // show when we constrcut ST or ST/Emax
bool bShowLeafLabels	= true;
bool bSaveSTLabels		= false; // toggle whether we write taxon labels or node numbers to GML files
bool bShowContradicted	= false; // show whether a node is contradicted
bool bShowFreq			= false; // show frequency of edge in ST/Emax
bool bWeighted          = false;
bool bWritePostscript   = false;
bool bWriteNEXUS		= false;
bool bWriteNewick		= false;
bool bVerbose			= false;
bool bWriteDot			= false; // show ST and ST/Emax as dot files
bool bWriteGML          = false; // Write graphs using GML format
bool bWriteMRP			= false; // Write MRP
bool bClusterGraph		= false; // Make cluster graph
bool bWriteTS			= false; // Output trees at each step in the recursion
bool bShowFan			= false;

int level 				= 0;
int graph_count			= 0;

int cluster_k			= 2; // number shared leaves for cluster graph

int supertreeNodeNumber = 0;


#include "stgraph.h"

bool showColours = false;


typedef struct {
	int level;
	int	nodes;
	int trees;
	bool connected;
	int cut;
	int components;
} Info;

Info info;


/**
 * @var  vector <NTree> NTreeVector
 * @brief A vector of NTree's, used to store the input trees
 *
 */
typedef vector <NTree> NTreeVector;


/**
 * @var  STree superTree
 * @brief The supertree
 *
 */
STree superTree;

STGraph CO;


/** 
 * @fn void MakeSTEmax (STGraph &ST, int wsum, Profile<NTree> &p)
 * @brief Construct the graph @f$S_T /E_T^{\max }@f$ from @f$S_T@f$
 *
 * @param ST the graph @f$S_T@f$
 * @param wsum the sum of weights for all source trees
 * @param p the multiset of trees
 *
 * Construct the graph @f$S_T /E_T^{\max }@f$ from @f$S_T@f$ by contracting all edges
 * in @f$S_T@f$ that have the maximum weight @f$w_{{\rm sum}}  = \sum\limits_{T \in T} {w(T)}@f$.
 *
 * We construct @f$S_T /E_T^{\max }@f$ by identifying edges in the graph that represent
 * uncontradicted nestings. In the Semple and Steel algorithm these edges are just those edges
 * with weight @f$w_{{\rm sum}}  = \sum\limits_{T \in T} {w(T)}@f$. Any edge not meeting this
 * criterion is "hidden" using graph::hide. The resulting graph is disconected, and each component
 * represents the set of edges to collapse. We find the components, then restore all the hidden
 * edges. We then iterate over the nodes in each component, merging all edges connecting 
 * all pairs of nodes in that component. Each component is then represented by a single node. The set of
 * merged nodes is stored in the node's node set.
 */
void MakeSTEmax (STGraph &ST, int wsum, NTreeVector &T, Profile<NTree> &p);

/**
 * @fn MinCutEdges (STGraph &ST, int cG)
 * @brief find and delete all edges in at least one minimum weight cut set
 *
 * @param ST the graph @f$S_T /E_T^{\max }@f$
 * @param cG the minimum-weight cut of @f$S_T /E_T^{\max }@f$
 *
 * Proposition 4.1 of Semple and Steel states that an edge e of a graph G
 * is in a minimum-weight cut set if and only if 
 * @f${\rm c}(G\backslash e) + w(e) = {\rm c}(G)@f$.
 *
 * MinCutEdges finds all edges that are in at least one minimum weight cut 
 * of the graph @f$S_T /E_T^{\max }@f$, and "deletes" them by calling
 * graph::hide_edge.
 */
void MinCutEdges (STGraph &ST, int cG);
void MinCutSupertree (NTreeVector &T, Profile<NTree> &p);

/**
 * @fn AllMinCuts (STGraph &ST, list<node_pair> &st_list)
 * @brief find and delete all edges in at least one minimum weight cut set
 *
 * @param ST the graph @f$S_T /E_T^{\max }@f$
 * @param st_lst a list of all (s,t) pairs for which the minimum s-t cut
 * is also a global minimum cut
 *
 * AllMinCuts finds all edges that are in a least one minimum cut by doing
 * the following for each (s,t) pair:
 *   -# compute the maximum s-t flow 
 *   -# construct the residual graph @f$R$@f
 *   -# find all strongly connected components of @f$R$@f 
 *   -# any edge in @f$R$@f with its ends in different strongly
 *      connected components belongs to some minimum cut
 * 
 * This algorithm is based on Picard and Queryanne.
 *
 * Having found the edges we "delete" them by calling
 * graph::hide_edge.
 */
void AllMinCuts (STGraph &ST, list<node_pair> &st_list);
/**
 * @fn MakeClusterGraph 
 * @brief Make cluster graph for input trees
 *
 * Following Sanderson et al. 1998, we create a graph where the nodes
 * are the input trees, and a pair of nodes are connected by an edge if the number
 * of taxa shared by the two nodes is greater than or equal to the
 * threshold k. The default value of k=2 is the minimum necessary to build
 * a supertree. If the cluster graph is not connected we find and output
 * the components.
 *
 */
void MakeClusterGraph (NTreeVector &T, int k = 2);
/**
 * @fn MakeCOGraph 
 * @brief Make graph of co-occurrences of all leaves
 *
 */
void MakeCOGraph (NTreeVector &T, Profile<NTree> &p);
/**
 * @fn WriteMRP
 * @brief Write input trees to a NEXUS file encoded for MRP
 *
 * The binary encoded trees are written to a transposed matrix
 * in the same style as Joe Thorley's RadCon program. This gretaly
 * simplifies writing the matrix.
 *
 */
void WriteMRP (ostream &f, NTreeVector &T, Profile<NTree> &p);

/*
The idea that Semple and Steel only collapse edges that are unanimously supported.
Many edges are uncontradicted by any trees, but the corresponding nestings will not
be reflected in the mincut tree. The proposal is this:

1. We first perform SEmple and Steel collapsing (all edges with w == wsum are
merged)
2. We then identify all edges which are uncontradicted, that is, their frequency of
occurrence in ST is the same as the frequency which which the two nodes cooccur in
the input set of trees.
3. The goal is to merge as many uncontradicted edges as possible. 


*/
void CollapseGraph (STGraph &ST, STGraph &fan);


#define TEST_1				1
#define WRITE_CLIQUES 		0
#define SHOW_REDIRECTED 	0
#define SHOW_CONTRADICTED 	0
#define SHOW_COMMONS 		0	
#define SHOW_COLOURS		0
#define SHOW_COMPONENTS		0
#define SHOW_BROKEN			0

#define USE_CLIQUE			0
#define USE_OTHER			1

//------------------------------------------------------------------------------
void CollapseGraph (STGraph &ST, STGraph &fan)
{
	graph::edge_iterator cit = ST.edges_begin();
	graph::edge_iterator cend = ST.edges_end();
	int num_contradicted_edges = 0;

	// Partition edges into those that are uncontradicted and those that
	// are contradicted
	while (cit != cend)
	{
		node n1 = cit->source ();
		node n2 = cit->target ();

		int freq_nested = ST.f[*cit];
		int freq_co = CO.GetEdgeFreqFromNodeLabels (ST.node_labels[n1], ST.node_labels[n2]);
		int freq_fan = fan.GetEdgeFreqFromNodeLabels (ST.node_labels[n1], ST.node_labels[n2]);
		
		int conflict = freq_co - freq_nested - freq_fan;

		if (conflict != 0)
		{
			ST.edge_colour[*cit] = colour_contradicted;
#if SHOW_CONTRADICTED
			cout << "Contradicted edge " << (*cit) << endl;
#endif
			num_contradicted_edges++;
		}
		else
			ST.edge_colour[*cit] = colour_uncontradicted;
			
		cit++;
	}
	
//	cout << (ST.number_of_edges() - num_contradicted_edges) << " edges of graph are uncontradicted" << endl;
	
#if TEST_1
	// Edges that are adjacent to a contradicted edge should not
	// be merged into a component as they are conflicts that need to be
	// cut later. We colour these edges "colour_adjto_contradicted"
	cit = ST.edges_begin();
	while (cit != cend)
	{
		node s = cit->source ();
		node t= cit->target ();
		
		if (ST.edge_colour[*cit] == colour_contradicted)
		{
		
//			cout << "contradicted edge " << (*cit) << endl;
			NodeSet c;
		
		
			// List all nodes adjacent to s that are not linked
			// by a contradicted edge
			NodeSet adjacent_to_s;
			node::adj_edges_iterator nit = s.adj_edges_begin();
			node::adj_edges_iterator nend = s.adj_edges_end();
			while (nit != nend)
			{
				if (ST.edge_colour[*nit] != colour_contradicted)
					adjacent_to_s.insert (s.opposite (*nit));
				nit++;
			}	
			// List all nodes adjacent to t that are not linked
			// by a contradicted edge
			nit = t.adj_edges_begin();
			nend = t.adj_edges_end();
			while (nit != nend)
			{
				if (ST.edge_colour[*nit] != colour_contradicted)
				{
					NodeSet::iterator n = adjacent_to_s.find (t.opposite (*nit));
					if (n != adjacent_to_s.end())
					{
						c.insert (*n);
						ST.edge_colour[*nit] = colour_adjto_contradicted;
//						cout << " adj: " << (*nit) << endl;
					}
				}
				nit++;
			}	
			nit = s.adj_edges_begin();
			nend = s.adj_edges_end();
			while (nit != nend)
			{
				NodeSet::iterator n = c.find (s.opposite (*nit));
				if (n != c.end())
				{
					ST.edge_colour[*nit] = colour_adjto_contradicted;
//					cout << " adj: " << (*nit) << endl;
				}
				nit++;
			}	
			
						
		}
		
		cit++;
	}
	
#if SHOW_COLOURS
		{
//			char buf[256];
//			sprintf (buf, "colours%d.gml", (graph_count-1));
//			ST.save (buf);
			char buf[256];
			sprintf (buf, "colours%d.dot", (graph_count-1));
			ST.mShowColours = true;
			ST.WriteDotty (buf);
			ST.mShowColours = false;
		}
#endif
	
#endif	

#if USE_OTHER
	// 0. Hide all explicitly contradicted edges
	list<edge> hidden_edges;
	cit = ST.edges_begin();
	cend = ST.edges_end();
	while (cit != cend)
	{
		if (ST.edge_colour[*cit] == colour_contradicted)
		{
			edge tmp = *cit;
			cit++;
			hidden_edges.push_back (tmp);
			ST.hide_edge (tmp);
		}
		else
			cit++;
	}
	
	if (!ST.is_connected())
	{
		// We can cut only contradicted edges
		superTree.GetCurNode()->AppendLabel("h");
	}
	else
	{
		cit = ST.edges_begin();
		cend = ST.edges_end();
		while (cit != cend)
		{
			if (ST.edge_colour[*cit] != colour_uncontradicted)
			{
				edge tmp = *cit;
				cit++;
				hidden_edges.push_back (tmp);
				ST.hide_edge (tmp);
			}
			else
				cit++;
		}
	}
	
	if (!ST.is_connected())
	{
		// 2. Get components
	    components cp;
	    if (cp.check(ST) != algorithm::GTL_OK) 
		{
			cerr << "component check failed at line " << __LINE__ << endl;
			exit(1);
	    } 
		else 
		{
			if (cp.run(ST) != algorithm::GTL_OK) 
			{
		    	cerr << "component algorithm failed at line " << __LINE__ << endl;
				exit(1);
			} 
			else 
			{
//				cout << "Graph has " << cp.number_of_components() << " components" << endl;

				// 3. Restore hidden edges
				list<edge>::iterator lit = hidden_edges.begin();
				list<edge>::iterator lend = hidden_edges.end();
				while (lit != lend)
				{
					ST.restore_edge (*lit);
					lit++;
				}
				
				// 4. List components and merge member nodes
				components::component_iterator it = cp.components_begin ();
				components::component_iterator end = cp.components_end ();
				list <node> acomponent;
				while (it != end)
				{
					acomponent = (*it).first;
					list<node>::iterator anode = acomponent.begin();
					list<node>::iterator first = anode;					
					list<node>::iterator last_node = acomponent.end();
					set<std::string> vertices;
					while (anode != last_node)
					{
					
					#if SHOW_COMPONENTS
						std::copy (ST.ns[*anode].begin(), ST.ns[*anode].end(),
							std::ostream_iterator<node>(cout, " "));
					#endif		
										
						if (anode != first)
						{
							ST.mergeNodes (*first, *anode); 
						}							
										
						anode++;
					}
					
					// Hide the merged nodes
					
					// There was a subtle bug here for novice STL programmers like me.
					// I want to hide all nodes except *first. To do this I iterate over
					// the list of nodes in *first's node set and hide them, except for
					// *first itself. Originally I did this by skipping the first element
					// in the node set like this
					//
					//    sit++;
					//    while (sit != send)
					//
					// This assumes that the first element in the set is always the first added
					// and this need not be the case. New code starts from the begining of the set
					// and explicitly tests whether the current element is the node *first.

					NodeSet::iterator sit = ST.ns[*first].begin();
					NodeSet::iterator send = ST.ns[*first].end();
					while (sit != send)
					{
						if ((*sit) != (*first))
							ST.hide_node(*sit);
						sit++;
					}
					
					
					
		
					
					
					#if SHOW_COMPONENTS
					cout << endl;
					#endif
					
					it++;
				}
			}
		}
		
#if SHOW_BROKEN
		{
//			showColours = true;
			char buf[256];
			sprintf (buf, "broken%d.dot", (graph_count-1));
			ST.WriteDotty (buf);
//			cout << "broken ST bwritten to " << buf << endl;
//			showColours = false;
		}
#endif
		
	}		
		


#endif


#if USE_CLIQUE


	// 2. Use Tseng's algorithm to find maximal clique partitions
	int count = 0;
	bool done = false;
	NodeSet common_set;
	int mostcommons = -1;
	edge eij;
	
	while (!done)
	{
		if (mostcommons == -1)
		{
			graph::edge_iterator eit = ST.edges_begin();
			graph::edge_iterator eend = ST.edges_end ();
			// Find edge with most common neighbours
			while (eit != eend)
			{
				if (ST.edge_colour[*eit] == colour_uncontradicted)
				{
					NodeSet c;

					node s = eit->source();
					node t = eit->target();
					// If the two nodes connected by *eit have fewer
					// edges than the current highest number of common
					// neighbours then they can't have more common neighours
					// than the current maximum
					if (min (s.degree(), t.degree()) >= mostcommons)
					{				
						int cij = ST.GetCommonNeighbours (*eit, c);
						if (cij > mostcommons)
						{
							mostcommons = cij;
							eij = *eit;
							common_set = c;
						}
					}
				}
				eit++;
			}
		}
		if (mostcommons == -1)
			done = true;
		else
		{

#if SHOW_COMMONS
			cout << "Edge with most commons is " << eij << " " << mostcommons << endl;
#endif
			node s = eij.source();
			node t = eij.target();


			// Add t to s's node set
			NodeSet::iterator nsit = ST.ns[t].begin();
			NodeSet::iterator nsend = ST.ns[t].end();
			while (nsit != nsend)
			{
				ST.ns[s].insert (*nsit);
				nsit++;
			}

			// Visit all edges adjacent to s and delete those that do not connect
			// common neighbours of s and t 
			int number_neighbours = 0;
			list<edge> to_be_deleted;
			node::adj_edges_iterator it = s.adj_edges_begin();
			node::adj_edges_iterator end = s.adj_edges_end();
			while (it != end)
			{
				// Node v is adjacent to s
				node v = s.opposite(*it);

				NodeSet::iterator n = common_set.find (v);
				if (n == common_set.end())
				{
					// n is not a common neighbour of s and t. If this edge
					// is contradicted we keep it, otherwise delete it
					if (ST.edge_colour[*it] == colour_uncontradicted)
						to_be_deleted.push_back (*it);
						
//					cout << " adj to s " << (*it) << endl;
				}
				else
					number_neighbours++;
				it++;
			}


			// Visit all edges adjacent to t and delete them (t is being merged 
			// with s)
			it = t.adj_edges_begin();
			end = t.adj_edges_end();
			list<edge> to_redirect;
			while (it != end)
			{
				if (t.opposite(*it) != s)
				{		
//					cout << " adj to t " << (*it) << endl;
					switch (ST.edge_colour[*it])
					{
						case colour_uncontradicted: // merged as part of clique
							to_be_deleted.push_back (*it); 
							break;
						case colour_adjto_contradicted: 
							if (common_set.find (t.opposite(*it)) == common_set.end())
								to_redirect.push_back (*it);
							else
								to_be_deleted.push_back (*it); 
							break;
						case colour_contradicted:
							if (ST.EdgeExists (s, t.opposite(*it)))
								to_be_deleted.push_back (*it);
							else	
								to_redirect.push_back (*it); 
							break;
					}
							
/*						
					if (ST.edge_colour[*it] == colour_uncontradicted)
					{
						to_be_deleted.push_back (*it); 
					}
					else if (ST.edge_colour[*it] == colour_contradicted)
					else if (ST.edge_colour[*it] == colour_contradicted)
						to_redirect.push_back (*it);*/
				}
				it++;

			}

			// Any remaining edges adjacent to t are contradicted edges that
			// need to be redirected so that they link to s
			list<edge>::iterator lit = to_redirect.begin();
			list<edge>::iterator lend = to_redirect.end();
			while (lit != lend)
			{
#if SHOW_REDIRECTED
				cout << "redirect " << (*lit);
#endif
				edge e = (*lit);
				node v = t.opposite(e);
				e.change_source (s);
				e.change_target (v);
#if SHOW_REDIRECTED
				cout << " to " << (*lit) << endl;
#endif
				lit++;
			}


			// Delete the extra edges
			lit = to_be_deleted.begin();
			lend = to_be_deleted.end();
			while (lit != lend)
			{
				// Think about this some more
//				if (ST.edge_colour[*lit] == colour_contradicted)
//					cout << "*** Should not delete this edge " << (*lit) <<	endl;
				ST.del_edge (*lit);
				lit++;
			}
			ST.hide_node (t);
			
			// Select next edge to merge. If node s is not isolated then
			// we pick and edge adjacent to s that has the largest number
			// of common neighbours. If s is isolated then mostcommons 
			// remains at -1, triggering a search through the whole graph
			mostcommons = -1;
			if (number_neighbours > 0)
			{
			
				node::adj_edges_iterator eit = s.adj_edges_begin();
				node::adj_edges_iterator eend = s.adj_edges_end ();
				// Find edge with most common neighbours
				while (eit != eend)
				{
					if (ST.edge_colour[*eit] == colour_uncontradicted)
					{
						NodeSet c;
						int cij = ST.GetCommonNeighbours (*eit, c);
						if (cij > mostcommons)
						{
							mostcommons = cij;
							eij = *eit;
							common_set = c;
						}
					}
					eit++;
				}
			}

			// Output current partition (debugging)	
			count++;
#if WRITE_CLIQUES
			char buf[256];
			sprintf (buf, "clique%d.gml", count);
			ST.save (buf);
			sprintf (buf, "clique%d.dot", count);
			ST.WriteDotty (buf);
#endif

		}
	
	}
#endif 

}


// Sketch of algorithm


//------------------------------------------------------------------------------
void MakeSTEmax (STGraph &ST, int wsum, NTreeVector &T, Profile<NTree> &p)
{
	// Step 1: Simple Semple and Steel
	//
	// Any nodes connected by an edge e for which w(e) = wsum are
	// merged.
	
	// 1. Get list of all edges w(e) == wsum
	list <edge> edges;
	
	graph::edge_iterator eit = ST.edges_begin();
	graph::edge_iterator eend = ST.edges_end ();
	while (eit != eend)
	{
		if (ST.w0[*eit] == wsum)
		{
			edges.push_back (*eit);
		}										
		eit++;
	}


	// 2. Collapse these edges	
	while (edges.size() != 0)
	{
		edge e = edges.front();
		
		// Merge nodes s and t where lit=e(s,t)
		node s = e.source();
		node t = e.target();
		ST.ns[s].insert (t);
			
		node::adj_edges_iterator it;
		node::adj_edges_iterator end;
   		it = t.adj_edges_begin();
    	end = t.adj_edges_end();
    	
    	list<edge> to_be_deleted;
    	// When w(e) == sum, any edge adjacent to t is also adjacent to s,
    	// so we delete that edge
		while (it != end)
		{
			node v = t.opposite (*it);
			
			if (v != s)
			{
				// If this edge has w(e) == wsum then delete from list of 
				// edges to be considered
				if (ST.w0[*it] == wsum)
					edges.remove (*it);
					
				to_be_deleted.push_back (*it);								
										
			}
			it++;
		}
		
		// Delete the extra edges
		list<edge>::iterator lit = to_be_deleted.begin();
		list<edge>::iterator lend = to_be_deleted.end();
		while (lit != lend)
		{
			ST.del_edge (*lit);
			lit++;
		}
		
		edges.pop_front();
		
		// Delete the edge e(s,t)
		ST.del_edge (e);
		// Node t is now merged with s (i.e., is an element of s's node set).
		// Hide this node so the graph object ignores t
		ST.hide_node (t);				
	}
	

	
	
	if (use_algorithm == ALGORITHM_ROD1)
	{
	
		// fans
		STGraph fan;
		fan.make_undirected();
		
		for (int i = 0; i < T.size(); i++)
		{
			NNodePtr root = (NNodePtr)T[i].GetRoot();
			if (root->GetDegree() > 2)
			{
				if (bShowFan)
					cout << "fan" << endl;
			
    			T[i].BuildLabelClusters ();
				T[i].Update();
				if (bShowFan)
					T[i].Draw (cout);
	
				NNodePtr n1 = (NNodePtr)root->GetChild();	
				while (n1->GetSibling())
				{
					NNodePtr n2 = (NNodePtr)n1->GetSibling();
					while (n2)
					{
						IntegerSet::iterator n1it = n1->Cluster.begin();
						IntegerSet::iterator n1end = n1->Cluster.end();
						while (n1it != n1end)
						{					
							IntegerSet::iterator n2it = n2->Cluster.begin();
							IntegerSet::iterator n2end = n2->Cluster.end();
							while (n2it != n2end)
							{
								if (bShowFan)
									cout << (*n1it) << "-" << (*n2it) << endl;
								fan.AddEdge (p.GetLabelFromIndex ((*n1it)-1), p.GetLabelFromIndex ((*n2it)-1));
								n2it++;
							}
							n1it++;
						}
						n2 = (NNodePtr)n2->GetSibling();
					}
					n1 = (NNodePtr)n1->GetSibling();
				}
			}
		}
		{	
			char buf[256];
			sprintf (buf, "colours%d.gml", (graph_count-1));
			fan.save (buf);
		}	

		CollapseGraph (ST, fan);
	}
	

	
	if (bShowSTEmax)
		cout << "ST-EMax" << endl << ST << endl;
	if (bSaveSTEmax)
	{
		char buf[64];

		if (bWriteGML)
		{
			sprintf (buf, "STEmax%d.gml", (graph_count-1));
			ST.save (buf);
		}

		if (bWriteDot)
		{
			sprintf (buf, "STEmax%d.dot", (graph_count-1));
			ofstream f (buf);
			ST.WriteDotty (f);
			f.close ();
		}

	}	
}

//------------------------------------------------------------------------------
void MinCutEdges (STGraph &ST, int cG)
{
	
	list<node_pair> st_list;

	//  true if edge is in at least one minimum cut set of ST
	edge_map <bool> in_a_min_cut_set (ST, false);

	// Semple and Steel brute force		


	// To avoid any problems with hide and restore edge operations affecting
	// the edge iterators of the graph while we go through the graph, store 
	// a separate list of the egdes of ST
	list <edge> edges;
	graph::edge_iterator eit = ST.edges_begin();
	graph::edge_iterator eend = ST.edges_end ();
	while (eit != eend)
	{
		edges.push_back (*eit);
		eit++;
	}
	

	// Visit each edge in ST and compute c(ST\e). Do this by hiding the
	// edge, rather than actually deleting it.
	list<edge>::iterator lit = edges.begin();
	list<edge>::iterator lend = edges.end();
	while (lit != lend)
	{
		if (ST.w0[*lit] <= cG)
		{
			// minimum cut weight without this edge
			ST.hide_edge (*lit);
			int cGe = mincut_st (ST, ST.w0, st_list);
		
			if (bShowAllMinCuts)
				cout << " c(ST\\e) for " << *lit << "= " << cGe;	
			
			ST.restore_edge (*lit);

			// is edge in a minimum cut set?
			if (bShowAllMinCuts)
				cout << " w(e) = " << ST.w0[*lit];
			
			in_a_min_cut_set [*lit] = (cGe + ST.w0[*lit] == cG);
		
			if (bShowAllMinCuts)
			{
				if (in_a_min_cut_set [*lit])
					cout << " [in a cut set]";
				cout << endl;
			}
		}
		else
		{
				if (bShowAllMinCuts)
					cout << ".";
		}
		lit++;
	}	

	// Delete edges in at least one cut set
	lit = edges.begin();
	lend = edges.end();
	while (lit != lend)
	{
		if (in_a_min_cut_set [*lit])
			ST.hide_edge (*lit);
		lit++;
	}
}

bool bShowstlist = false;
bool bShowFlow = false;
bool bShowResiduals = false;
bool bShowStrong = false;

//------------------------------------------------------------------------------
void AllMinCuts (STGraph &ST, list<node_pair> &st_list)
{
	list<node_pair>::iterator st;
	
	if (bShowstlist)
	{
		// Show list of (s,t) pairs
		cout << st_list.size() << " (s,t) mincut pair(s) found:" << endl;
		st = st_list.begin();
		while (st != st_list.end())
		{
			cout << "   " << (*st).first << "," << (*st).second << endl;
			st++;
		}
		cout << endl;
	}
	
	// To find all minimum cuts of G we need to:
	// 1. For each (s,t) pair compute the maximum flow
	// 2. Get the strong components of the residual graph for the flow
	// 3. Identify the edges in the cut(s)


	// Boolean flag to indicate if an edge is in a cut
	edge_map <bool> in_a_min_cut_set (ST, false);
	
	// Map between edges and original edges in graph
	edge_map <edge> orig;

	// To compute a flow we need a directed graph, so we need to convert our undirected
	// graph into a directed graph by "splitting" the undirected edges into a pair
	// of directed edges pointing in opposite directions.
	ST.make_directed();
	edge_map<double> capacity;

	// Store list of original edges in the graph
	list<edge> l;
	edge e;
	forall_edges (e, ST)
	{
		l.push_back (e);
		orig[e] = e;
	}
	
	// Go through each edge and create a reversed copy
	list<edge>::iterator it = l.begin();
	while (it != l.end())
	{
		node source = it->source();
		node target = it->target();
		edge e = ST.new_edge (target, source); // reverse direction to *it
		ST.w0[e] = ST.w0[*it];
		
		orig[e] = *it;

		// Set edge capacity to the weight of the edge
		capacity[*it] = (double)ST.w0[e];
		capacity[e] = capacity[*it];
		
		it++;
	}

	//ST.save ("double.gml");
		
			
	// Iterate over each (s,t) pair 
	int st_count = 0;		
	st = st_list.begin();
	while (st != st_list.end())
	{
//		cout << (*st).first << "," << (*st).second << endl;
		
		// Find maximum (s,t) flow
		maxflow_ff ff;
	
		ff.set_vars (capacity, (*st).first, (*st).second);
		
	    if (ff.check(ST) != algorithm::GTL_OK) 
		{
			cerr << "maxflow_ff check failed" << endl;
			exit(1);
	    } 
		else 
		{
			if (ff.run(ST) != algorithm::GTL_OK) 
			{
		    	cerr << "maxflow_ff algorithm failed" << endl;
				exit(1);
			} 
			else 
			{
				if (bShowFlow)
					cout << "max flow = " << ff.get_max_flow() << endl;
				
				// ff.run restores the graph (sigh), which means the hidden nodes reappear.
				// Since some hidden nodes are merged nodes, I don't want these to 
				// be visible as they are unconnected and hence bugger the next iteration
				// of maximum flow
				{
					graph::node_iterator n_it, n_end, n_tmp;

					n_it = ST.nodes_begin();
					n_end = ST.nodes_end();

					while (n_it != n_end) 
					{
						n_tmp = n_it;
						++n_tmp;
						if (n_it->degree() == 0)
							ST.hide_node(*n_it);
						n_it = n_tmp;
					}
				}
				
				
				if (bShowFlow)
				{
					char buf[64];
					sprintf (buf, "flow_%d.gml", st_count);
					ST.save (buf);	
					st_count++;
				}
				
	
				// Residual capacities.
				// Hide any edge with a non positive residual capacity
	
				if (bShowFlow)
					cout << "   find residual flows..." << endl;
					
				edge e;
				list<edge> to_hide;
				forall_edges (e, ST)
				{
					if (ff.get_rem_cap (e) <= 0.0)
						to_hide.push_back(e);
				}
				list<edge>::iterator it = to_hide.begin();
				while (it != to_hide.end())
				{
					ST.hide_edge (*it);
					it++;
				}
				
				if (bShowResiduals)
				{
					char buf[64];
					sprintf (buf, "residual_%d.gml", st_count);
					ST.save (buf);	
					st_count++;
				}
							
				// get strong components
				node_map<int> cp(ST,-1);
				int scc = STRONG_COMPONENTS (ST, cp);
				
				if (bShowStrong)	
					cout << "   graph has " << scc << " strong components" << endl;
				
				/*{
					node n;
					forall_nodes (n, G0)
						cout << n << cp[n] << endl;
				}*/
					
				// Edges in mincut
				// Corollary 6 of Picard and Queryanne states that an edge
				// is in a mincut iff its ends do not lie in the same
				// strongly connected component
				if (bShowStrong)
					cout << "Identify edges in a cut" << endl;
				forall_edges (e, ST)
				{
					if (cp[e.source()] != cp[e.target()])
					{
						in_a_min_cut_set[orig[e]] = true;
						if (bShowStrong)
							cout << e << "is in a mincut" << endl;
					}
				}
				
				
				// Restore hidden edges
				it = to_hide.begin();
				while (it != to_hide.end())
				{
					ST.restore_edge (*it);
					it++;
				}
				
			}
		}
		
		
		st++;
	}
	

	// Restore original undirected graph

	// 1. Hide all the edges
	{
		edge e;
		list<edge> to_hide;
		forall_edges (e, ST)
		{
			to_hide.push_back(e);
		}
		list<edge>::iterator it = to_hide.begin();
		while (it != to_hide.end())
		{
			ST.hide_edge (*it);
			it++;
		}
	}
	
	// Restore only those original edges not in a cut set
	it = l.begin();
	while (it != l.end())
	{
		
		if (!in_a_min_cut_set [*it])
		{
			ST.restore_edge (*it);
		}	
		it++;
	}
	ST.make_undirected();
}


//------------------------------------------------------------------------------
void MinCutSupertree (NTreeVector &T, Profile<NTree> &p)
{
	int wsum = 0;
	
	STGraph ST;
	ST.make_undirected();

	if (bShowLevel)
		cout << "-------------- LEVEL " << level << "----------------" << endl;

	// 3. Construct ST
	if (bShowConstruct)
		cout << "Construct ST" << endl;
		
	for (int i = 0; i < T.size(); i++)
	{
    	T[i].BuildLabelClusters ();
		T[i].Update();
		
		wsum += (int)T[i].GetWeight(); // to do: this assumes integer weights!!
		
		if (bShowTrees)
			T[i].Draw (cout);
			
		if (bShowClusters)
		{
			cout << "Clusters" << endl;
			T[i].ShowClusters ();
		}

		NNodePtr n = (NNodePtr)T[i].GetRoot();
		n = (NNodePtr)n->GetChild();
		while (n != NULL)
		{
			if (n->IsLeaf())
			{
				ST.AddNode (n->GetLabel());
			}
			else
			{
				if (bShowSTEdge)
				{
				/*	cout << "{ ";
					std::copy (n->Cluster.begin(), n->Cluster.end(),
						std::ostream_iterator<int>(cout, " "));
					cout << "}" << endl;*/
				}

				IntegerSet::iterator iit = n->Cluster.begin();
				IntegerSet::iterator iend = n->Cluster.end();
				while (iit != iend)
				{
					IntegerSet::iterator jit = iit;
					jit++;
					while (jit != iend)
					{
						if (bShowSTEdge)
						{
							cout << p.GetLabelFromIndex ((*iit)-1) << "-" 
								<< p.GetLabelFromIndex ((*jit)-1) << endl;
						}
						// 5 Nov 2001
						// Edges are weighted by tree weights
						ST.AddEdge (p.GetLabelFromIndex ((*iit)-1), p.GetLabelFromIndex ((*jit)-1), (int)T[i].GetWeight());
						jit++;
					}
					iit++;
				}
			}
			n = (NNodePtr)n->GetSibling();
		}
	}
	
	info.level = graph_count; 
	info.nodes = ST.number_of_nodes();
	info.trees = T.size();

	if (bShowST)
		cout << "ST" << endl << ST << endl;
	if (bSaveST)
	{
		char buf[64];
		
		ST.mShowLabels = bSaveSTLabels;
		
		if (bWriteGML)
		{
			sprintf (buf, "ST%d.gml", graph_count);
			ST.save (buf);
		}
		
		
		if (bWriteDot)
		{
			sprintf (buf, "ST%d.dot", graph_count);
			ofstream f (buf);
			ST.WriteDotty (f);
			f.close ();
		}
		
	}
	graph_count++;

	int minimumCut = 0;
	if (ST.is_connected())
	{
		// If ST is connected then we construct ST/Emax. This merges nodes that are part of
		// a clique of nodes with maximally weighted edges. We delete all edges that are in
		// a minimal cut of the graph
		
		info.connected = true;
		
		// 5.
		if (bShowConnected)
			cout << "ST is connected so constructing ST/Emax" << endl;
		MakeSTEmax (ST, wsum, T, p);
		
		list<node_pair> st_list;

		minimumCut = mincut_st (ST, ST.w0, st_list);
		
		info.cut = minimumCut;
		
		if (bShowMinCutWeight)
			cout << "Minimum-weight cut of ST/Emax = " << minimumCut << " yielding ";
			
#if 1
		// All mincuts algorithm
		AllMinCuts (ST, st_list);
	
#else				
		// Semple and Steel brute force
		MinCutEdges (ST, minimumCut);		
#endif	

		char numbuf[16];
		sprintf (numbuf, "c%d", minimumCut);
		superTree.GetCurNode()->AppendLabel (numbuf);	

	}
	else
	{
		info.connected 	= false;
		info.cut 		= 0;
		
		// 4.
		if (bShowConnected)
			cout << "ST is not connected" << endl;
			
		superTree.GetCurNode()->AppendLabel ("c0");
	}

	// The vertex sets are the components of ST
    components cp;
    if (cp.check(ST) != algorithm::GTL_OK) 
	{
		cerr << "component check failed at " << __LINE__ << " in " << __FILE__  << endl;
		exit(1);
    } 
	else 
	{
		if (cp.run(ST) != algorithm::GTL_OK) 
		{
	    	cerr << "component algorithm failed at " << __LINE__ << " in " << __FILE__ << endl;
			exit(1);
		} 
		else 
		{
			// Show info -------------------------------------------------------
			info.components = cp.number_of_components();
			cout  << setiosflags (ios::right)
				<< setw (8) << info.level 
				<< setw (8) << info.nodes
				<< setw (8) << info.trees;
			if (info.connected)
				cout << "     yes" << setw (8) << info.cut;
			else
				cout << "      no        ";
			cout << setw (16) << info.components << endl;
			
			// Vist each component ---------------------------------------------
			if (bShowMinCutWeight)
				cout  << cp.number_of_components() << " components" << endl;
			if (bShowVertexSets)
				cout << endl << "ST has " << cp.number_of_components() << " vertex sets:" << endl;

			components::component_iterator it = cp.components_begin ();
			components::component_iterator end = cp.components_end ();
			list <node> acomponent;
			while (it != end)
			{
				acomponent = (*it).first;
				list<node>::iterator anode = acomponent.begin();
				list<node>::iterator last_node = acomponent.end();
				set<std::string> vertices;
				while (anode != last_node)
				{
//					cout << (*anode);
//					std::copy (ST.ns[*anode].begin(), ST.ns[*anode].end(),
//						std::ostream_iterator<node>(cout, " "));

					// Get set of nodes associated with this node in ST
					NodeSet::iterator nsit = ST.ns[*anode].begin();
					NodeSet::iterator nsend = ST.ns[*anode].end();
					while (nsit != nsend)
					{
						vertices.insert (ST.node_labels[*nsit]);
						nsit++;
					}

					anode++;
				}

				// Handle the vertex set
				set <std::string>::iterator nsit = vertices.begin();
				set <std::string>::iterator nsend = vertices.end();
				
				if (bShowVertexSets)
				{
					while (nsit != nsend)
					{
						cout << (*nsit) << " ";		
						nsit++;
					}
					cout << endl;
				}
				
				if (vertices.size() < 3)
				{
                	// The first component is a child of the current
                    // node in the growing supertree, the other
                    // components are siblings of the first component
                    if (it == cp.components_begin ())
                        superTree.MakeChild();
                    else
                        superTree.MakeSibling();

					nsit = vertices.begin ();
					
					if (vertices.size() == 1)
					{
						// a leaf
                        superTree.AddLeaf (*nsit);
					}
					else
					{
						// a cherry
                        std::string label1 = (*nsit);
                        nsit++;
                        std::string label2 = (*nsit);
                        superTree.AddCherry (label1, label2);
					}
				}
				else
				{
                	// Component has more than two leaves and hence needs
                    // further analysis

					// Construct the vector of trees T|S, i.e. the subtree of T that contains
					// only leaves in S
					NTreeVector TS;
					for (int i = 0; i < T.size(); i++)
					{
						NTree t = T[i];
						t.BuildLeafLabels ();
    					t.BuildLabelClusters ();

    					// Get set of leaves in T
						IntegerSet tset = ((NNodePtr)t.GetRoot())->Cluster;
						t.Update ();
						
    					if (bShowTS)
    					{
 							cout << " leaves = " << t.GetNumLeaves() << endl;
	    					cout << t << endl;
						/*	std::copy (tset.begin(),tset.end(),
								std::ostream_iterator<int>(cout, " "));*/
							cout << endl;
						}


						LabelMap::iterator it = p.Labels.begin();
						LabelMap::iterator end = p.Labels.end();
						while (it != end)
						{
							std::string s = (*it).first;	// label
							int index = (*it).second;		// 0-offset index of label
							index++;						// make index 1-offset
							
							if (bShowTStest)
								cout << s << "=" << index;
							
							if (tset.find (index) != tset.end())
							{
								// label is in original tree
								if (bShowTStest)
									cout << " in t";
								
								if (vertices.find (s) == vertices.end())
								{	
									if (bShowTStest)	
										cout << " not in v";										
									// but not in this vertex set	
									NodePtr q = t.leaf_labels[s];			
									t.RemoveNode (q);
									delete q;
									if (bShowTStest)
									{
										cout << " leaves = " << t.GetNumLeaves() << endl;
										cout << t << endl;
									}

								}
							}
							if (bShowTStest)
								cout << endl;
							it++;
						}
						
						if (bShowTStest)
							cout << " leaves = " << t.GetNumLeaves() << endl;
						
						// Only add tree to the list if it has some leaves
						if (t.GetNumLeaves() > 0)
						{				
							t.Update();
							if (bShowTS)
								t.Draw (cout);
						
							// Build clusters
    						t.BuildLabelClusters ();
							TS.push_back (t);
						}
					}
					
					// process T|S ---------------------------------------------
					if (TS.size () > 0)
					{
						if (TS.size() == 1)
						{
							// Here we can use a shortcut. If only a single tree
                            // has leaves in the current vertex set then we don't need to
                            // recursively do mincuts, we simply graft the corresponding
                            // subtree onto the growing supertree.
							superTree.AddSubtree (TS[0], (it == cp.components_begin ()));
						}
						else
						{				
							// More than one tree has leaves in the current vertex set,
                            // so we find the mincut supertree for the set of subtrees.
                            // This is where the algorithm becomes recursive.
							level++;
							if (bShowRecursion)
								cout << "--> MinCutSupertree" << endl;
								
                            // The first component is a child of the current
                            // node in the growing supertree, the other
                            // components are siblings of the first component
                            if (it == cp.components_begin ())
                                superTree.MakeChild();
                            else
                                superTree.MakeSibling();

							superTree.PushNode ();
							
							// Output trees for debugging
							if (bWriteTS)
							{
								char buf[256];
								sprintf (buf, "ts-%d.tre", (graph_count - 1));
								ofstream f (buf);
								for (int i = 0; i < TS.size(); i++)
								{
									NTree t = TS[i];
									f << t << endl;
								}
								f.close();
							}

							MinCutSupertree (TS, p);
							if (bShowRecursion)	
								cout << "<-- MinCutSupertree" << endl;
								
                            superTree.PopNode();

							level--;	
						}
					}		
				} // if (vertices.size() < 3)

				it++; // next component
			}
		}
	}
}


//------------------------------------------------------------------------------
void MakeCOGraph (NTreeVector &T, Profile<NTree> &p)
{
	CO.make_undirected();
		
	// For each tree insert an edge in CO between pairs of taxa that
	// cooccur in a tree
	for (int i = 0; i < T.size(); i++)
	{
   		T[i].BuildLabelClusters ();
		T[i].Update();		
		NNodePtr n = (NNodePtr)T[i].GetRoot();

		IntegerSet::iterator iit = n->Cluster.begin();
		IntegerSet::iterator iend = n->Cluster.end();
		while (iit != iend)
		{
			IntegerSet::iterator jit = iit;
			jit++;
			while (jit != iend)
			{
				CO.AddEdge (p.GetLabelFromIndex ((*iit)-1), p.GetLabelFromIndex ((*jit)-1));
				jit++;
			}
			iit++;
		}
	}

	{
		if (bWriteGML)
		{
			CO.save ("CO.gml");
			cout << "CO written to CO.gml"  << endl;
		}
		if (bWriteDot)
		{
			ofstream f ("CO.dot");
			CO.WriteDotty (f);
			f.close ();
			cout << "CO written to CO.dot"  << endl;
		}
	}

}

//------------------------------------------------------------------------------
void MakeClusterGraph (NTreeVector &T, int k)
{
	cout << "Computing cluster graph with k=" << k << endl;

	STGraph ClusterGraph;
	
	ClusterGraph.make_undirected();
	
	// A map between integers and nodes allows us to index nodes in 
	// ClusterGraph by the order of occurrence of the tree in the
	// set of input trees
	map <int, node, less <int> > t;
	
	// A map between nodes and integers allows us to recover the tree
	// associated with a node in ClusterGraph
	node_map<int> node_to_tree_number (ClusterGraph, 0);
	
	// Create a node for each tree
	for (int i = 0; i < T.size(); i++)
	{
		node n = ClusterGraph.new_node();
		t[i] = n;
		node_to_tree_number[n] = i;
		ClusterGraph.node_labels[n] = T[i].GetName();
	}
		
	
	// Store the leaves in each tree	
	vector <IntegerSet> leaves;
	for (int i = 0; i < T.size(); i++)
	{
   		T[i].BuildLabelClusters ();
		T[i].Update();	
		NNodePtr n = (NNodePtr)T[i].GetRoot();
		leaves.push_back (n->Cluster);
	}
	
	// Compare all trees and create a edge between corresponding nodes in
	// ClusterGraph if the number of shared leaves exceeds threshold
	for (int i = 1; i < T.size(); i++)
	{

		for (int j = 0; j < i; j++)
		{
			IntegerSet common;
			std::set_intersection (leaves[i].begin(), leaves[i].end(),
				leaves[j].begin(), leaves[j].end(),
				std::inserter (common, common.end()));
			int w = common.size();
			if (w >= k)
			{
				edge e = ClusterGraph.new_edge (t[i], t[j]);
				ClusterGraph.w0[e] = w;
			}
		
		}
	}
	
	
	{
		if (bWriteGML)
		{
			ClusterGraph.save ("cluster.gml");
			cout << "   Cluster graph written to \"cluster.gml\""  << endl;
		}
		if (bWriteDot)
		{
			ClusterGraph.WriteDotty ("cluster.dot");
			cout << "   Cluster graph  written to \"cluster.dot\""  << endl;
		}
	}


	// How many components are there?
	if (ClusterGraph.is_connected())
	{
		cout << "   Cluster graph is connected" << endl;	
	}
	else
	{
		cout << "Warning: Cluster graph for k=" << k << " is not connected" << endl;
		cout << "(if k<2 input trees cannot be combined into a supertree)" << endl;
		
    	components cp;
    	if (cp.check(ClusterGraph) != algorithm::GTL_OK) 
		{
			cerr << "component check failed for ClusterGraph" << endl;
			exit(1);
    	} 
		else 
		{
			if (cp.run(ClusterGraph) != algorithm::GTL_OK) 
			{
	    		cerr << "component algorithm failed for ClusterGraph" << endl;
				exit(1);
			} 
			else 
			{
				cout  << "Cluster graph has " << cp.number_of_components() << " components:" << endl;
				cout << "Comp.\tSize\tFile\tMembers" << endl;
				int count = 0;
				char fname[256];
				components::component_iterator it = cp.components_begin ();
				components::component_iterator end = cp.components_end ();
				list <node> acomponent;
				while (it != end)
				{
					cout << ++count << "\t";
					acomponent = (*it).first;
					cout << acomponent.size() << "\t";
					
					// Output trees
					sprintf (fname, "cluster.k%d.%d.%d.tre", k, count, (int)acomponent.size());
					cout << fname << "\t";

					ofstream f (fname);
					f << "#nexus" << endl;
					f << "[! Cluster number " << count << " containing " << acomponent.size() << " trees]" << endl;
					f << "begin trees;" << endl;
					
					list<node>::iterator anode = acomponent.begin();
					list<node>::iterator last_node = acomponent.end();
					while (anode != last_node)
					{
						cout << " " << ClusterGraph.node_labels[*anode];
						
						f << "\ttree " << ClusterGraph.node_labels[*anode] << " = [&R] ";
						f << "[&W " << T[node_to_tree_number[*anode]].GetWeight() << "] ";
						T[node_to_tree_number[*anode]].Write (f);
						
						f << endl;

						anode++;
					}
					
					cout << endl;
					
					f << "end;" << endl;
					f.close();
					it++; // next component
				}
			}
		}
	}

}


//------------------------------------------------------------------------------
// Use same format as Joe Thorley uses in RadCon
void WriteMRP (ostream &f, NTreeVector &T, Profile<NTree> &p)
{
	f << "#nexus" << endl << endl;
	f << "[MRP file written "; 
	time_t timer = time(NULL);
	struct tm* tblock = localtime(&timer);
	char time_buf[64];
	strncpy (time_buf, asctime(tblock), sizeof (time_buf));
	char *q = strrchr (time_buf, '\n');
	if (q)
		*q = '\0';
	f << time_buf << "]" << endl;
	f << endl;

	f << "begin taxa;" << endl;
	int ntax = p.GetNumLabels();
	f << "\tdimensions ntax=" << (ntax + 1) << ";" << endl;
	f << "\ttaxlabels" << endl;
	f << "\t\tmrp_outgroup" << endl;
	// Output leaf labels
	for (int i = 0; i < ntax; i++)
		f << "\t\t'" << p.GetLabelFromIndex (i) << "'" << endl;
	f << "\t\t;" << endl;
		
	f << "end;" << endl << endl;
	
		
	// Store info on character set corresponding to each tree
	std::vector <int> start;
	std::vector <int> end;

	
	// Count the number of characters
	int nchar = 0;
	for (int i = 0; i < T.size(); i++)
	{
		start.push_back (nchar+1);
		nchar += T[i].GetNumInternals();
		end.push_back (nchar);
	}


	f << "begin characters;" << endl;
	f << "\tdimensions nchar=" << nchar << ";" << endl;
	f << "\tformat symbols=\"01\" missing=? transpose nolabels;" << endl;
	f << "\tmatrix" << endl;
		
	for (int i = 0; i < T.size(); i++)
	{
   		T[i].BuildLabelClusters ();
		T[i].Update();	
		
		f << "[" << T[i].GetName() << "]" << endl;
		
		// Iterator over nodes in tree
		NNodePtr r = (NNodePtr)T[i].GetRoot();
		NodeIterator <NNode> ni (r);
		NNodePtr n = ni.begin();
		while (n != NULL)
		{
			if (!n->IsLeaf())
			{
				// Fill character string with ?
				std::string s (ntax, '?');
				
				// Only leaves in the root cluster will be coded for this tree
				IntegerSet::iterator iit = r->Cluster.begin();
				IntegerSet::iterator iend = r->Cluster.end();
				while (iit != iend)
				{
					s[*iit-1] = '0';
					iit++;
				}

				// Code any leaf in this cluster as "1", nodes not in this cluster
				// will be left as "0" 
				 iit = n->Cluster.begin();
				 iend = n->Cluster.end();
				 while (iit != iend)
				 {
				 	s[*iit-1] = '1';
					iit++;

				 }
				 
				 f << "0"; // mrp outgroup state is always 0
				 f << s << endl;
			}
			n = ni.next();
		}
	}
	f << "\t;" << endl;
	f << "end;" << endl;	
	f << endl;
	
	f << "begin sets;" << endl;
	f << "\t[charsets corresponding to binary codes for each tree]" << endl;
	for (int i = 0; i < T.size(); i++)
	{
		f << "\tcharset " << T[i].GetName() << " = " << start[i] << "-" << end[i] << ";" << endl;
	}
	f << "end;" << endl;
 	

}


//------------------------------------------------------------------------------
int main (int argc, char **argv)
{

#if __MWERKS__
	#if macintosh
		argc = ccommand(&argv);
	#endif
#endif

	// Parse options
    // Heavily borrowed from the squid library
    char *optname;
    char *optarg;
    int   optind;
    
    // default settings
    bSaveSTLabels	 = false;
	bSaveST 		 = false;
	bSaveSTEmax		 = false;
	bWritePostscript = false;
	bWriteNEXUS		 = false;
	bWriteNewick	 = false;	
	bVerbose		 = false;
	bWeighted		 = false;
	bWriteDot		 = false;
	bWriteMRP		 = false;
	use_algorithm    = ALGORITHM_ROD1;
	bClusterGraph	 = false;
	cluster_k		 = 2;
	
	char ps_name[FILENAME_SIZE];
	char nxs_name[FILENAME_SIZE];	
	char nwk_name[FILENAME_SIZE];
	char mrp_name[FILENAME_SIZE];


    while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
    &optind, &optname, &optarg))
    {
    	if (strcmp(optname, "-l") == 0) {  bSaveSTLabels = true; }
    	if (strcmp(optname, "-b") == 0) {  bVerbose = true; }
      	if (strcmp(optname, "-w") == 0) {  bWeighted = true; }
    	if (strcmp(optname, "-d") == 0) {  bWriteDot = true; bWriteGML = false; }    	
    	else if (strcmp(optname, "-g") == 0) {  bSaveST = bSaveSTEmax = true; }	
    	else if (strcmp(optname, "-p") == 0) 
    	{  
    		bWritePostscript = true; 
			strcpy( ps_name, optarg);   		
    	}
		else if (strcmp(optname, "-a") == 0)
		{
			use_algorithm = atoi(optarg);
		}
		else if (strcmp(optname, "-c") == 0)
		{
			bClusterGraph = true;
			cluster_k = atoi(optarg);
			if (cluster_k < 0)
			{
				cerr << "k must be greater than zero" << endl;
				exit (0);				
			}
		}
    	else if (strcmp(optname, "-n") == 0) 
    	{  
    		bWriteNEXUS = true; 
		strcpy( nxs_name, optarg);   		
    	}
		else if (strcmp(optname, "-m") == 0) 
    	{  
    		bWriteMRP = true; 
		strcpy( mrp_name, optarg);   		
    	}
     	else if (strcmp(optname, "-k") == 0) 
    	{  
    		bWriteNewick = true; 
			strcpy( nwk_name, optarg);   		
    	}
        else if (strcmp(optname, "-v") == 0) 
        { 
        	cout << "Supertree version 0.1" << endl;
			cout << "   Uses GTL " << GTL_MAJOR_VERSION << "." << GTL_MINOR_VERSION << "." << GTL_MINI_VERSION
				<< " (c) University of Passau" << endl;
            exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind != 1)
    {
        cerr << "Incorrect number of arguments:" << usage << endl;
        exit (0);
    }

	// Get options from command line
	char fname[FILENAME_SIZE];
	strcpy( fname, argv[optind++] );
	
	// Check file exists
	FILE* file = fopen( fname, "r" );
	if( file != NULL ) 
	{
		fclose(file);
		file = NULL;
	}
	else
	{
		cerr << "File \"" << fname << "\" does not exist." << endl;
        exit (0);
    }

	cout << "Mincut supertree version " << MAJOR_VERSION << "." MINOR_VERSION << "." << MINI_VERSION << endl;

	ifstream f (fname);

	bShowST				= false;
	bShowSTEmax			= false;
	bShowAllMinCuts 	= true;
	bShowTrees 			= false;
	bShowMinCutWeight	= true;
	bShowVertexSets		= false;
	bShowRecursion		= false;
	bShowTStest			= false;
	bShowConnected		= false;
	bShowTS				= false;
	bShowClusters		= false;
	bShowSTEdge			= false;
	bShowLevel			= false;
	bShowConstruct		= false;
	bShowLeafLabels		= false;

	if (!bVerbose)
	{
		bShowAllMinCuts = false;	
		bShowMinCutWeight = false;
	}
	
	Profile<NTree> p;

    if (!p.ReadTrees (f))
	{
		cerr << "Failed to read trees, bailing out" << endl;
		exit(0);
	}
	
    if (bVerbose)
    	p.ShowTrees (cout);
    p.MakeLabelFreqList ();
    
    if (bShowLeafLabels)
    	p.ShowLabelList (cout);

	// Create initial multiset of trees T
	NTreeVector T;
	int minLeaves = 1000; // ugh
	int maxLeaves = 0;

	for (int i = 0; i < p.GetNumTrees(); i++)
	{
		NTree t = p.GetIthTree (i);

		// We need to build clusters using the same labels across
		// all trees
		t.MakeNodeList();
		for (int j = 0; j < t.GetNumLeaves(); j++)
			t[j]->SetLabelNumber (p.GetIndexOfLabel (t[j]->GetLabel()) + 1);
		minLeaves = min (minLeaves, t.GetNumLeaves());
		maxLeaves = max (maxLeaves, t.GetNumLeaves());
		
		T.push_back (t);
	}
	
	// k-cluster graph
	if (bClusterGraph)
		MakeClusterGraph (T, cluster_k);	

	// Get co-occurrences of taxa
	MakeCOGraph (T, p);
	
	// MRP matrix
	if (bWriteMRP)
	{
		ofstream mrpfile (mrp_name);
		WriteMRP (mrpfile, T, p);
		mrpfile.close ();
		cout << "MRP file written to \"" << mrp_name << "\"" << endl;
		return 0;	
	}
	


	// Output the options
	if (bSaveST)
	{
		cout << "   ST graphs are saved in GML";
		if (bWriteDot)
			cout << " and .dot";
		cout << " format with nodes labelled";
		if (bSaveSTLabels)
			cout << " with leaf names";
		else
			cout << " with numbers";
		cout << " to files ST[n].gml (.dot)" << endl;
	}
	if (bSaveSTEmax)
	{
		cout << "   ST/Emax graphs are saved in GML";
		if (bWriteDot)
			cout << " and .dot";
		cout << " format with nodes labelled";
		if (bSaveSTLabels)
			cout << " with leaf names";
		else
			cout << " with numbers";
		cout << " to files STEmax[n].gml (.dot)" << endl;
	}
	
	
	
	
//	cout << "i\tnodes\tconnected\tcut\tcomponents" << endl;
	
//	        "-------+-------+-------+-------+-------+-------+------+
	cout << endl;
	cout << "       i   nodes  trees     con?     cut     components" << endl;   
	cout << "-------------------------------------------------------" << endl;


	clock_t t1 = clock();
	
	

    // Initialise supertree
    // Create the root and push it on the stack
    superTree.MakeRoot();
    superTree.PushNode ();
    superTree.SetInternalLabels (true);

	MinCutSupertree (T, p);

    superTree.PopNode ();
    
	cout << "-------------------------------------------------------" << endl;
	cout << "         i: step in algorithm" << endl;
	cout << "     nodes: number of nodes in ST graph" << endl;
	cout << "     trees: number of trees with leaves in T|S" << endl;
	cout << "      con?: is ST connected? (yes/no)" << endl;
	cout << "       cut: weight of minimum cut of ST/Emax" << endl;
	cout << "components: components in ST" << endl;

	clock_t t2 = clock();
	
	if (bVerbose)
	{
		// Show the tree
		cout << "Supertree:" << endl;
		cout << superTree << endl;
	}
	
	if (bWriteNewick)
	{
		ofstream nwkfile (nwk_name);
		nwkfile << superTree;
		nwkfile.close ();	
	}	

	
	if (bWriteNEXUS)
	{
		ofstream nxsfile (nxs_name);

		// Output the results in NEXUS format
		nxsfile << "#nexus" << endl;
		nxsfile << endl;

		nxsfile << "begin trees;";
		
		// Date the file
		nxsfile << " [Treefile written ";
		time_t timer = time(NULL);
		struct tm* tblock = localtime(&timer);
		char time_buf[64];
		strncpy (time_buf, asctime(tblock), sizeof (time_buf));
		char *q = strrchr (time_buf, '\n');
		if (q)
			*q = '\0';
		nxsfile << time_buf << "]" << endl;
		
		// Source tree info
		nxsfile << "[!" << endl;
		nxsfile << "   Mincut supertree computed from " << T.size() << " trees in file \"" << fname << "\"" << endl;
		nxsfile << "      Smallest source tree had " << minLeaves << " leaves" << endl;
		nxsfile << "      Largest source tree had " << maxLeaves << " leaves" << endl;
		if (bWeighted)
			nxsfile << "   Input tree weights used" << endl;

		switch (use_algorithm)
		{		
			case ALGORITHM_SEMPLE:
				nxsfile << "   Semple and Steel (2000)";
				break;
			case ALGORITHM_ROD1:
				nxsfile << "   Rod Page unpublished";
				break;
			default:
				nxsfile << "   Undocumented";
				break;
		}
		nxsfile << " algorithm used" << endl;
		
		if (bWeighted)
			nxsfile << "   Trees weighted using user-supplied weights" << endl;
		else
			nxsfile << "   Trees all equally weighted" << endl;

		nxsfile << endl;
		nxsfile << "   supertree version " << MAJOR_VERSION << "." MINOR_VERSION << "." << MINI_VERSION << endl;
		// GTL details
       	nxsfile << "   GTL version " << GTL_MAJOR_VERSION << "." << GTL_MINOR_VERSION << "." << GTL_MINI_VERSION
				<< " (c) University of Passau" << endl;
 		nxsfile << endl;
	
		// Time to compute tree
		nxsfile << "   CPU time used = " << (t2 - t1)/CLOCKS_PER_SEC << " seconds" << endl;
		nxsfile << "]" << endl;
		nxsfile << "\ttree min_cut_supertree = [&R] ";
		nxsfile << superTree << endl;
		nxsfile << "end;" << endl;
		
		nxsfile.close();
	}
	
	// Postscript tree
	if (bWritePostscript)
	{
		superTree.Update();
	    Port.StartPicture (ps_name);
	    GBaseFont fo;
	    GRect r;
	    Port.GetPrintingRect (r);
		r.SetRight (r.GetRight() - 100);
		superTree.Plot (r, fo, TS_LEFT | TS_CLADOGRAM | TS_SLANT, 1);
	    Port.EndPicture ();
	    cout << "Supertree written to \"" << ps_name << "\"" << endl;

	}
	
	cout << endl <<  "CPU time used = " << (t2 - t1)/CLOCKS_PER_SEC << " seconds" << endl;

	
	// Compute measure of similarity between supertree and each input tree
	// using triplets
	if (1)
	{
		double sum_fit = 0.0;
		double weighted_count = 0.0;
		bool lots = true;
		
		// Crucial step. We need to ensure that leaves in the supertree are numbered
		// sensibly (i.e., 1..n). If we have grafted a subtree on then this need not be the case
	    NodeIterator <Node> n (superTree.GetRoot());
	    Node *q = n.begin();
	    int count = 0;
	    while (q)
	    {
	    	if (q->IsLeaf ())
	        {
	       		count++;
	        	q->SetLeafNumber (count);
	        }
	        q = n.next();
	    }
	    superTree.MakeNodeList();
		
		if (lots)
		{
			cout << endl << "Fit statistics" << endl;
			cout << "    tree  weight  leaves       n       d       s      r1      r2     fit   name" << endl;
			cout << "-------------------------------------------------------------------------------" << endl;
		}
			
		for (int i = 0; i < p.GetNumTrees(); i++)
		{
			// We are comparing each input tree with the subtree induced in the supertree
			// by the set of leaves in the input tree. This is destructive so we make a 
			// copy of the supertree.
			NTree t1 (superTree);
			// Get ith input tree
			NTree t2 = p.GetIthTree (i);
			
			t1.MakeNodeList();
			t2.MakeNodeList();
			
			
			// We need to prune excess leaves from supertree, and ensure that leaves in the pruned
			// supertree have their LeafNumber field in the range 1...n so the triplet comparison
			// will work correctly. Note that leaves in t2 are already numbered 1...n
	        IntegerSet toPrune;
	        int count = 0;
	        int j = 0;
	        while (j < t1.GetNumLeaves())
	        {
	        	NodePtr matchingLeaf = t2.GetLeafWithLabel (t1[j]->GetLabel());
	        	if (matchingLeaf != NULL)
	            {
	            	// This leaf is in t1 AND t2. Set the LeafNumber of this leaf
	            	// in t1 to match that in t2
	            	count++;
					t1[j]->SetLeafNumber(matchingLeaf->GetLeafNumber());
	            }
	            else
	            {
	            	// This leaf is not in t2 so we will prune it from the supertree
	            	toPrune.insert (j);
	            }
	            j++;
	        }

			// Prune excess leaves from t1
            IntegerSet::iterator nit = toPrune.begin();
            IntegerSet::iterator nend = toPrune.end();
            while (nit != nend)
            {
            	t1.RemoveNode (t1[*nit]);
                delete t1[*nit];
                nit++;
            }

			// Build clusters in the two trees, and rebuild node list so we
			// can look up leaves
   			t1.BuildLeafClusters();
   			t2.BuildLeafClusters();
			t1.MakeNodeList();
			t2.MakeNodeList();
				
			QTValues Q;
			CompareTriplets (t1, t2, Q);
			SummaryStats (Q);

				
			double different_in_t2 = Q.d + Q.r2;
			double resolved_in_t2 = Q.d + Q.s + Q.r2;
			double t2_fit = 1.0 - different_in_t2/resolved_in_t2;
			
			
			if (lots)					
				cout  << setiosflags (ios::right)
					<< setw (8) << (i+1)
					<< setw (8) << t2.GetWeight()
					<< setw (8) << t2.GetNumLeaves()
					<< setw (8) << Q.n
					<< setw (8) << Q.d
					<< setw (8) << Q.s
					<< setw (8) << Q.r1
					<< setw (8) << Q.r2
					<< setprecision (3) << setw (8) << t2_fit
					<< "   " << t2.GetName()
					<< endl;
			
			sum_fit += (double)(t2.GetWeight()) * t2_fit;
			weighted_count += (double)(t2.GetWeight());
			
		}
		if (lots)
		{
			cout << "-------------------------------------------------------------------------------" << endl;
			cout << "weight: weight of tree i";
			if (!bWeighted) 
				cout << " (note that tree weights were not used)";
			cout << endl;
			cout << "leaves: number of leaves in tree i" << endl;
			cout << "     n: total number of triplets" << endl;
			cout << "     d: triplets resolved differently in supertree and tree i" << endl;
			cout << "     s: triplets resolved identically in supertree and tree i" << endl;
			cout << "    r1: triplets resolved in supertree but not in tree i" << endl;
			cout << "    r2: triplets resolved in tree i but not in supertree" << endl;
			cout << "   fit: 1 - (d + r2)/(d + s + r2)" << endl;
			cout << endl;
		}
		cout << "Average fit: " << sum_fit/weighted_count << endl;
		
	}
	
	
    return 0;
}
