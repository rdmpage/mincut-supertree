#include "strong_components.h"

#include <stack>

typedef stack<node> node_list;

static void scc_dfs(const graph& G, node v, node_map<int>& compnum,
                                            node_map<int>& dfsnum,
                                            node_list& unfinished,
                                            node_list & roots,
                                            int& count1, int& count2, node_map<bool> &is_unfinished );

static void scc_dfs(const graph& G, node v, node_map<int>& compnum,
                                            node_map<int>& dfsnum,
                                            node_list& unfinished,
                                            node_list & roots,
                                            int& count1, int& count2, node_map<bool> &is_unfinished )
{
	node w;

	dfsnum[v] = ++count1;
	unfinished.push(v);
	is_unfinished[v] = true;
	roots.push(v);

	node::adj_nodes_iterator it = v.adj_nodes_begin();
	while (it != v.adj_nodes_end())
    { 
    	node w = *it;
		if (dfsnum[w]==-1) 
			scc_dfs(G,w,compnum,dfsnum,unfinished,roots,count1,count2, is_unfinished);
		else 
			if (is_unfinished[w])
				while (dfsnum[roots.top()] > dfsnum[w])  
					roots.pop();
		it++;
	}

	if (v == roots.top()) 
	{ 
		do { 
			w = unfinished.top();
			unfinished.pop();
			is_unfinished[w] = false;
          	/* w is an element of the scc with root v */
          	compnum[w] = count2;
		} while (v != w);
     	count2++;
     	roots.pop(); 
	}
}




  // int STRONG_COMPONENTS(graph& G, node_m<int>& compnum)
  // computes strong connected components (scc) of digraph G
  // returns m = number of scc 
  // returns in node_array<int> compnum for each node an integer with
  // compnum[v] = compnum[w] iff v and w belong to the same scc
  // 0 <= compnum[v] <= m-1 for all nodes v

int STRONG_COMPONENTS(const graph& G, node_map<int>& compnum)
{

	node_list      roots;
	node_list      unfinished;
	node_map<bool> is_unfinished (G, false);
	node_map<int> dfsnum(G,-1);

 	int count1 = 0; 
	int count2 = 0;

	node v;

	forall_nodes(v,G) 
		if (dfsnum[v] == -1) 
			scc_dfs(G,v,compnum,dfsnum,unfinished,roots,count1,count2, is_unfinished);

	return count2;
}
