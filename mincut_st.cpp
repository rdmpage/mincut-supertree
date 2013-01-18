// $Id: mincut_st.cpp,v 1.1.1.1 2001/10/31 10:43:03 rdmp1c Exp $

#include "mincut_st.h"

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <queue>
#include <list>
#include <map>
#include <set>

#ifdef __GNUC__
#include <algorithm>
#endif

#ifdef __BORLANDC__
	#include <values.h>
#endif
#if (defined __MWERKS__) || (defined __GNUC__)
	#include <climits>
	#define MAXINT INT_MAX
#endif

#include "fheap.h"


bool bShowOriginal 	= false;	
bool bShowCopy		= false;
bool bShowEdges		= false;	
bool bShowCut		= false;


int mincut_st (const graph &G0, edge_map <int> &w0, list<node_pair> &st_list)
{
	// If graph is not connected then minimum weight cut is zero, and we do not need to
	// compute a minimum cut. 
	if (!G0.is_connected())
		return 0;


	graph G;
    G.make_undirected();

	if (bShowOriginal)
	{
	    cout << "Original graph" << endl << G0 << endl;
	 	   
	}



    // Make a local copy of the graph as MinimumCut modifies the original graph

	// List of nodes in the original graph
	node_map <node> partner (G0);
	node_map <node> orig (G);

	node x;
	int counter = 0;
	forall_nodes (x, G0)
	{
		partner[x] = G.new_node(); 
		orig[partner[x]] = x; // so we can look up original node
	}

	// Create edges and associated weights
	edge_map<int> w(G, 0);
	edge e;
	forall_edges (e, G0)
	{
		if (e.source() != e.target())
		{
			edge ec = G.new_edge (partner[e.source()], partner[e.target()]);
			w[ec] = w0[e];
		}
	}

	// Display the graph
	if (bShowCopy)
	{
		cout << "Copy of graph" << endl << G << endl;
		G.save ("start.gml");
	}

	// Display edge weights
	if (bShowEdges)
	{
		cout << "Edge weights" << endl;
		graph::edge_iterator it = G.edges_begin();
		graph::edge_iterator end = G.edges_end();
	
		while (it != end)
		{
			cout << *it << " " << w[*it] << endl;
			it++;
		} 
	}

	// Start of algorithm. $a$ is an arbitrary single node in $G$. The set $A$
	// of nodes initially comprises $a$
	graph::node_iterator na = G.nodes_begin();
	node a = *na;
	int n = G.number_of_nodes();
	int cut_weight = MAXINT;
	int best_value = MAXINT;
	while (n >= 2 )
	{
		node t = a;
		node s, v;
		edge e;
   		node::adj_edges_iterator it;
		node::adj_edges_iterator end;
		
		fheap_t *pq = fh_alloc (n);
		node_map<int> vertex_number (G, 0);
		map <int, node, less<int> > nv;
		int vertex_count = 0;
			
		// Nodes in $A$ are not in the queue
		node_map<bool> in_PQ(G, false);
		forall_nodes (v, G)
		{
			vertex_number[v] = vertex_count;
			nv[vertex_count] = v;
			vertex_count++;
			if (v != a)
			{
				in_PQ[v] = true;
				fh_insert (pq, vertex_number[v], 0);	
			}
		}
		node_map<int> inf (G, 0); 
		// Get weight of edges adjacent to $a$
		it = a.adj_edges_begin();
		end = a.adj_edges_end();
		while (it != end)
		{
			v = a.opposite (*it);
			inf[v] += w[*it];	
			it++;
		}
		// Store weights in a queue
		it = a.adj_edges_begin();
		end = a.adj_edges_end();
		while (it != end)
		{
			v = a.opposite (*it);
			fh_decrease_key (pq, vertex_number[v], -inf[v]);  
			it++;
		}

		while (pq->n > 0)
		{
			s = t;

			// Get the node that is most tightly connected to $A$
			t = nv[fh_delete_min (pq)];
			cut_weight = inf[t];
			in_PQ[t] = false;

			// Increase the key of nodes adjacent to t and not in $A$ by adding the
			// weights of edges connecting t with nodes not in $A$ 
			it = t.adj_edges_begin();
			end = t.adj_edges_end();
			while (it != end)
			{
				v = t.opposite (*it);
				if (in_PQ[v])
				{
					inf[v] += w[*it];
					fh_decrease_key (pq, vertex_number[v], -inf[v]);  
				}
				it++;
			}	
		}
		fh_free (pq);
		
		
		if (bShowCut)
			cout << "   cut-of-the-phase = " << cut_weight << endl;
		
		if (cut_weight <= best_value)
		{
			if (cut_weight < best_value)
			{
				// Clear list of (s,t) pairs
				st_list.erase (st_list.begin(), st_list.end());
				best_value = cut_weight;
			}
			st_list.push_back (node_pair (orig[s], orig[t]));
		}

		// Nodes s and t are the last two nodes to be added to A

		if (bShowCut)
		{
			cout << "s=" << s << " t=" << t << endl;
					 
		}

		// Get list of edges adjacent to s
		edge dummy;
		node_map<edge> s_edge(G, dummy);
		it = s.adj_edges_begin();
		end = s.adj_edges_end();
		while (it != end)
		{
			s_edge[s.opposite(*it)] = *it;
			it++;
		}

		// Merge s and t
   		it = t.adj_edges_begin();
    	end = t.adj_edges_end();


		// Iterate over edges adjacent to t. If a node v adjacent to
		// t is also adjacent to s, then add w(it) to e(s,v)
		// otherwise make a new edge e(s,v)
		while (it != end)
		{
			v = t.opposite (*it);

			if (s_edge[v] != dummy)
			{
				w[s_edge[v]] += w[*it];
			}
			else if (s != v)
			{
				edge ne = G.new_edge (s, v);
				w[ne] = w[*it];
			}				
			it++;

			
		}

		// Delete node t from graph
		G.del_node(t);

		
		if (bShowCut)
		{
			char fname[256];
			sprintf (fname, "step%d.gml", n);
			G.save (fname);
			 
		}

		n--;
		
		

	}
		
    return best_value;

}
