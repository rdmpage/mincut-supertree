mincut-supertree
================

C++ program to compute modified mincut supertrees. For background see

Page, R. D. M. (2002). Lecture Notes in Computer Science. (R. Guigó & D. Gusfield, Eds.) (pp. 537–551). Springer-Verlag. [doi:10.1007/3-540-45784-4_41](http://dx.doi.org/10.1007/3-540-45784-4_41)

You can get a PDF of this article [here](http://darwin.zoology.gla.ac.uk/~rpage/supertree/wabi.pdf).

The original software was written in 2002, see [http://darwin.zoology.gla.ac.uk/~rpage/supertree/](http://darwin.zoology.gla.ac.uk/~rpage/supertree/). This version as been edited to compile under gcc 4. NEXUS parsing has been removed as that code doesn't compile (yet).

### Building

First install GTL (there is a version available [here](https://github.com/rdmpage/graph-template-library)). Then:

   make Makefile

### Command line options

-v	show version information	 
-g	write ST and ST/Emax to GML file(s)
-l	output taxon labels when writing graph files
-d	write ST and ST/EMax to dot files
-b	verbose	Print the input trees, and extra information about the progress of the algorithm to the screen
-n filename	 write NEXUS file
-w	use tree weights	By default each tree has equal weight. If weights are included in the NEXUS file, the -w option will use those weights
-m filename	write MRP matrix to file <filename>
-k filename	write Newick file
-a n choose algorithm, either 1 for Rod Page modification, or 0 for original Semple and Steel	 
-c n compute cluster graph for k=n	
The program will construct the cluster graph for the input value of k, then finds all components of that graph and outputs the trees belong to each component to a file called

cluster.<k>.<count>.<size>.tre

where count is the running count of the number of components of the graph, and size is the number of trees in that component. You can use this option to generate well connect sets of tres for further analysis. If the -g or -d options are on, then the cluster graph is written to the file "cluster.gml" or "cluster.dot", respectively.

-p filename	write tree to Postscript file <filename>

