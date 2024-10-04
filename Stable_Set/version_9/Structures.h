#ifndef Structures_h
#define Structures_h

// VERTEX
struct Vertex{
    int gNumber;    // global number
    int lNumber;    // local number --> when graph is reduced we renumber the vertices
    Vertex(int g, int l): gNumber(g), lNumber(l){}
};

// EDGE
struct Edge{
    int vertex1;
    int vertex2;
    Edge(int i,int j): vertex1(i), vertex2(j){}
};

// Parameters for the algorithm --> values read from file params
struct Params{
    
    // smallest graph for which we doo full enumerattion to find stable set
    int minSizeGraphSizeForFullEnum;
    
	int strategyBranchVar;  		// 1: random
		                            // 2: min
		                            // 3: max
		                            // 4: closest to 0.5
                            
    int initialBranchingStrategy; 	// 0: MIN --> select vertex with lowest degree
    								// 1: MAX --> select vertex with highest degree           
    
    int maxExecutionTime; 			// in seconds
};


#endif /* Structures_h */
