#include "stgraph.h"

#include <fstream>

#ifdef __GNUC__
	# if __GNUC__ == 3
		#include <iterator>
	# endif
#endif


//------------------------------------------------------------------------------
void STGraph::AddNode (std::string s)
{
	if (labelled_nodes.find (s) == labelled_nodes.end ())
	{
		node n = new_node ();
		labelled_nodes[s] = n;
		node_labels[n] = s;
	}
}

//------------------------------------------------------------------------------
void STGraph::AddEdge (std::string s1, std::string s2, int weight)
{
	node n1, n2;
	bool both_nodes_exist = true;

	if (labelled_nodes.find (s1) == labelled_nodes.end ())
	{
		both_nodes_exist = false;
		n1 = new_node ();
		labelled_nodes[s1] = n1;
		node_labels[n1] = s1;
	}
	else
	{
		n1 = labelled_nodes[s1];
	}

	if (labelled_nodes.find (s2) == labelled_nodes.end ())
	{
		both_nodes_exist = false;
		n2 = new_node ();
		labelled_nodes[s2] = n2;
		node_labels[n2] = s2;
	}
	else
	{
		n2 = labelled_nodes[s2];
	}

	if (both_nodes_exist)
	{	
		// Is there already an edge connecting n1 and n2?
	   	node::adj_edges_iterator eit = n1.adj_edges_begin();
		node::adj_edges_iterator eend = n1.adj_edges_end();
		node found = n1;
		while ((found == n1) && (eit != eend))
		{
			if (n1.opposite (*eit) == n2)
				found = n2;
			else
				eit++;
		}
		if (found == n2)
		{
			// Yes, so increment weight of edge
			w0[*eit] += weight;	
			f[*eit] += 1;
		
		}
		else
		{
			// No, so create new edge in graph
			edge e = new_edge (n1, n2);
			w0[e] = weight;
			f[e] = 1;
		}
	}
	else
	{
		// One or both nodes are new, so create new edge in graph
		edge e = new_edge (n1, n2);
		w0[e] = weight;
		f[e] = 1;
	}
	
}

//------------------------------------------------------------------------------
bool STGraph::EdgeExists (node n1, node n2)
{
	bool result = false;
	
	// Is there an edge connecting n1 and n2?
	node::adj_edges_iterator eit = n1.adj_edges_begin();
	node::adj_edges_iterator eend = n1.adj_edges_end();
	node found = n1;
	while ((found == n1) && (eit != eend))
	{
		if (n1.opposite (*eit) == n2)
			found = n2;
		else
			eit++;
	}
	if (found == n2)
	{
		result = true;		
	}
	return result;
}

//------------------------------------------------------------------------------
int STGraph::GetEdgeFreqFromNodeLabels (std::string s1, std::string s2)
{
	int result = 0;
	bool both_nodes_exist = true;
	node n1, n2;

	if (labelled_nodes.find (s1) == labelled_nodes.end ())
	{
		both_nodes_exist = false;
	}
	else
	{
		n1 = labelled_nodes[s1];
	}

	if (labelled_nodes.find (s2) == labelled_nodes.end ())
	{
		both_nodes_exist = false;
	}
	else
	{
		n2 = labelled_nodes[s2];
	}

	if (both_nodes_exist)
	{	
		// Is there an edge connecting n1 and n2?
	   	node::adj_edges_iterator eit = n1.adj_edges_begin();
		node::adj_edges_iterator eend = n1.adj_edges_end();
		node found = n1;
		while ((found == n1) && (eit != eend))
		{
			if (n1.opposite (*eit) == n2)
				found = n2;
			else
				eit++;
		}
		if (found == n2)
		{
			result = f[*eit];
		}
	}
	return result;
}



//------------------------------------------------------------------------------
void STGraph::post_new_node_handler (node n)
{
	graph::post_new_node_handler (n);
	ns[n].insert (n);
}

//------------------------------------------------------------------------------
void STGraph::save_graph_info_handler (ostream *os) const
{
	graph::save_graph_info_handler (os);
}

//------------------------------------------------------------------------------
void STGraph::save_edge_info_handler (ostream *os, edge e) const
{
	graph::save_edge_info_handler (os, e);	
	*os << "label \"" << w0[e];
//	if (bShowFreq)
//		*os << "[" << f[e] << "]";
	*os << "\"" << endl;

	// Line width 1 pt
	*os << "graphics [" << endl;
	
	if (mShowColours)
	{
	
		switch (edge_colour[e])
		{
			case colour_uncontradicted:
				*os << "width 1.0" << endl;
				break;
			case colour_contradicted:
				*os << "fill \"#ff0000\"";
				*os << "width 2.0" << endl;
				break;
			case colour_adjto_contradicted:
				*os << "fill \"#0000ff\"";
				*os << "width 2.0" << endl;
				break;
		}
	}
	else
		*os << "width 1.0" << endl;
		

	*os << "]" << endl;
	
	// Use standard Postscript font
	*os << "LabelGraphics [" << endl;
	*os << "type \"text\"" << endl;
	*os << "font \"Helvetica\"" << endl;
	*os << "]" << endl;	
}

//------------------------------------------------------------------------------
void STGraph::save_node_info_handler (ostream *os, node n) const
{
	graph::save_node_info_handler (os, n);	
	*os << "label \"";

	if (mShowLabels)
	{
		// Full names
		NodeSet nset = ns[n];
		NodeSet::iterator sit =nset.begin();
		NodeSet::iterator send = nset.end();
		while (sit != send)
		{
			*os << node_labels[*sit] << " ";
			sit++;
		}
	}
	else
	{
		// Numbers
	/*	std::copy (ns[n].begin(), ns[n].end(),
			std::ostream_iterator<node>(*os, " ")); */
	}
	*os << "\"" << endl;

	// Use standard Postscript font
	*os << "LabelGraphics [" << endl;
	*os << "type \"text\"" << endl;
	*os << "font \"Helvetica\"" << endl;
	*os << "]" << endl;
	
}


//------------------------------------------------------------------------------
int STGraph::GetCommonNeighbours (edge e, NodeSet &common_set)
{
	int count = 0;
	
	node s = e.source();
	node t = e.target();
	
	// List all nodes adjacent to s that are not linked
	// by a contradicted edge
	NodeSet adjacent_to_s;
	node::adj_edges_iterator nit = s.adj_edges_begin();
	node::adj_edges_iterator nend = s.adj_edges_end();
	while (nit != nend)
	{
		if (edge_colour[*nit] != colour_contradicted)
			adjacent_to_s.insert (s.opposite (*nit));
		nit++;
	}	
	// List all nodes adjacent to t that are not linked
	// by a contradicted edge
	nit = t.adj_edges_begin();
	nend = t.adj_edges_end();
	while (nit != nend)
	{
		if (edge_colour[*nit] != colour_contradicted)
		{
			NodeSet::iterator n = adjacent_to_s.find (t.opposite (*nit));
			if (n != adjacent_to_s.end())
			{
				count++;
				common_set.insert (*n);
			}
		}
		nit++;
	}	
	return count;
}

//------------------------------------------------------------------------------
void STGraph::mergeNodes (node s, node t)
{

	NodeSet::iterator nsit = ns[t].begin();
	NodeSet::iterator nsend = ns[t].end();
	while (nsit != nsend)
	{
		ns[s].insert (*nsit);
		nsit++;
	}

	
	// List all nodes adjacent to s
	NodeSet adjacent_to_s;
	node::adj_nodes_iterator nit = s.adj_nodes_begin();
	node::adj_nodes_iterator nend = s.adj_nodes_end();
	while (nit != nend)
	{
		adjacent_to_s.insert (*nit);
		nit++;
	}
	
	
	// Visit all edges adjacent to t. 
	list<edge> to_be_deleted;
	node::adj_edges_iterator it = t.adj_edges_begin();
	node::adj_edges_iterator end = t.adj_edges_end();
	while (it != end)
	{
		// Node u is adjacent to t
		node u = t.opposite(*it);
		
		if (u == s)
		{
			// Edge (s,t) will become a loop when s and t are merged
					
			// Action: delete edge (s,t)
			to_be_deleted.push_back (*it); 
			
//			cout << "   deleted loop " << (*it) << endl;
		
		}
		else
		{
			NodeSet::iterator in_as = adjacent_to_s.find (u);
			if (in_as != adjacent_to_s.end())
			{
				// Node u is adjacent to both s and t, and so edge (t,u) will become parallel with 
				// edge (s, u) when s and t are merged.
			
				// Action: delete edge (t,u)
				to_be_deleted.push_back (*it); 
				
//				cout << "   deleted parallel edge " << (*it) << endl;

			}
			else 
			{
				// Node u is adjacent to t but not s. When s and t are merged, u will become adjacent
				// to s. Hence we need an edge (s,u)
				
				// Action: delete edge (t, u)
				int saved_weight = w0[*it];
				to_be_deleted.push_back (*it); 
				
//				cout << "   Redirected edge " << (*it);			
				
				// Action: create edge (s,u)
				edge e = new_edge (s, u);
				w0[e] = saved_weight;
				
//				cout << " to  " << e << endl;
			}
		}		
		it++;

	}
	// Delete the extra edges
	list<edge>::iterator lit = to_be_deleted.begin();
	list<edge>::iterator lend = to_be_deleted.end();
	while (lit != lend)
	{
		del_edge (*lit);
		lit++;
	}

}

//------------------------------------------------------------------------------
void STGraph::WriteDotty (const char *fname)
{
	std::ofstream f (fname);
	WriteDotty (f);
	f.close ();
}

//------------------------------------------------------------------------------
void STGraph::WriteDotty (ostream &f)
{
	node_map <int> index;
	graph::node_iterator nit = nodes_begin();
	graph::node_iterator nend = nodes_end();
	int count = 0;
	while (nit != nend)
	{
		index[*nit] = count++;
		nit++;
	}
	
	f << "graph G {" << endl;

	// Try and make the graph look nice
	f << "   node [width=.1,height=.1,fontsize=8,style=filled,color=lightblue2];" << endl;
   	f << "   edge [fontsize=8,len=2];" << endl;
	
	// Write node labels
	nit = nodes_begin();
	while (nit != nend)
	{
		f << "   " << index[*nit] << " [label=\"";
		 
		if (mShowLabels)
		{
			// Full names
			NodeSet nset = ns[*nit];
			NodeSet::iterator sit =nset.begin();
			NodeSet::iterator send = nset.end();
			while (sit != send)
			{
				f << node_labels[*sit];
				sit++;
				if (sit != send)
					f << " ";
			}
		}
		else
		{
			// Numbers
		/*	std::copy (ns[*nit].begin(), ns[*nit].end(),
				std::ostream_iterator<node>(f, " ")); */
		}
		f << "\"];" << endl;

		nit++;
	}
	
	
	// Write edges labelled with edge weight
	graph::edge_iterator it = edges_begin();
	graph::edge_iterator end = edges_end();
	while (it != end)
	{
		f << "   " << index[it->source()] << " -- " << index[it->target()] 
			<< "[label=\"" << w0[*it] << "\"";
		if (mShowColours)
		{
		
			switch (edge_colour[*it])
			{
				case colour_uncontradicted:
					f << ",color=black";
					break;
				case colour_contradicted:
					f << ",color=red";
					break;
				case colour_adjto_contradicted:
					f << ",color=blue";
					break;
			}
		}

			
		f << "];" << endl;
		it++;
	}
	
	f << "}" << endl;
}
