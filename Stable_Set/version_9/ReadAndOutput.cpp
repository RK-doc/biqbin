#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "Structures.h"
#include "ReadAndOutput.h"


void read_params(Params &params){

    std::string s;
    std::ifstream paramFile("params_ss");
    
    paramFile >> s >> params.minSizeGraphSizeForFullEnum;
    paramFile >> s >> params.strategyBranchVar;
    paramFile >> s >> params.initialBranchingStrategy;
    paramFile >> s >> params.maxExecutionTime;
    
    paramFile.close();
}

int read_input(char *argv[], int &n, int &m, std::vector<Vertex> &vertices, std::vector<Edge> &edges){
    
    FILE *file = fopen(argv[1], "r");

    if (!file){
        std::cerr << "Problem opening the file.\n";
        return -1;
    }
    
    // read number of vertices
    if ( fscanf(file, "%d", &n) != 1 ){
        std::cerr << "Problem reading number of vertices.\n";
        return -1;
    }

    if( n <= 0 ){
        std::cerr << "Error! Number of vertices has to be positive.\n";
       return -1;
    }
    
    // read number of edges
    if ( fscanf(file, "%d", &m) != 1 ){
        std::cerr << "Problem reading number of edges.\n";
        return -1;
    }
    
    if( m < 0 ){
        std::cerr << "Error! Number of edges has to be nonnegative.\n";
        return -1;
    }

    // create list of vertices
    for(int i = 1; i <= n; i++){
        vertices.emplace_back(i,i);
    }
    
    // read edges and store in list, ignore weights
    int v1, v2, weight;
    for(int i = 0; i < m; i++){
        
        if ( fscanf(file, "%d %d %d", &v1, &v2, &weight ) != 3 ){
            std::cerr << "Problem reading edges of the graph.\n";
            return -1;
        }

        // check edge
        if(v1 < 1 || v1 > n){
            std::cerr << "Problems with edge. Vertex has to be between 1 and " << n << ".\n";
            return -1;
        }
        if(v2 < 1 || v2 > n){
            std::cerr << "Problems with edge. Vertex has to be between 1 and " << n << ".\n";
            return -1;
        }

        // add edge to list
        if(v1 < v2)
        {
            edges.emplace_back(v1,v2);
        }
        else
        {
            edges.emplace_back(v2,v1);
        }
    }

    fclose(file);
    return 1; 		// everything OK
}

void print_solution(int n, int m, int sol[], int val, int counterBaBNodes, double timeTotal, double density, bool optimal) {
	
	std::cout << "{" << std::endl;

	std::cout << "\t\"GraphMD\": {" << std::endl;
        std::cout << "\t\t\"ObjectTypeString\": \"Graph\"," << std::endl;
        std::cout << "\t\t\"Vertices\": " << n << "," << std::endl;
        std::cout << "\t\t\"Edges\": " << m << "," << std::endl;
        printf("\t\t\"Density\": %.3lf\n", density);

        std::cout << "\t}," << std::endl;

	std::cout << "\t\"ExecutionMD\": {" << std::endl;
	printf("\t\t\"ExecutionTime\": %.3lf,\n", timeTotal);
	
	if (optimal)
		std::cout << "\t\t\"SolutionType\": \"Optimal\"," << std::endl;
	else
		std::cout << "\t\t\"SolutionType\": \"Approximate\"," << std::endl;
	
	std::cout << "\t\t\"Solution\": " << val << "," << std::endl;
	
	int counter = 1;	
	std::cout << "\t\t\"SolutionSet\": \"(";
	for(int i = 0; i < n; i++){
        if ( (sol[i] == 1) ){
        	if (counter != val){
            	std::cout << i+1 << ", ";
            	++counter;
            }
            else
            	std::cout << i+1 << ")\"," << std::endl;	
        }    
    }
   
    std::cout << "\t\t\"BabNodes\": " << counterBaBNodes << std::endl;
    std::cout << "\t}" << std::endl;
	std::cout << "}" << std::endl;
	
}

void outputFile(char *argv[], int n, int m, int sol[], int val, int counterBaBNodes, double timeTotal, double density, bool optimal){

	std::string nameOfFile = argv[2];
	if (optimal)
		nameOfFile.append("_solution_optimal.txt");
	else
		nameOfFile.append("_solution_approximate.txt");

	std::ofstream fileOutput(nameOfFile);
	
	fileOutput << "{" << std::endl;


	fileOutput << "\t\"GraphMD\": {" << std::endl;
        fileOutput << "\t\t\"ObjectTypeString\": \"Graph\"," << std::endl;
        fileOutput << "\t\t\"Vertices\": " << n << "," << std::endl;
        fileOutput << "\t\t\"Edges\": " << m << "," << std::endl;
        fileOutput << "\t\t\"Density\": " << std::setprecision(3) << density << "\n";
        fileOutput << "\t}," << std::endl;


	fileOutput << "\t\"ExecutionMD\": {" << std::endl;
	fileOutput << "\t\t\"ExecutionTime\": " << std::setprecision(3) << timeTotal << ",\n";
	
	if (optimal)
		fileOutput << "\t\t\"SolutionType\": \"Optimal\"," << std::endl;
	else
		fileOutput << "\t\t\"SolutionType\": \"Approximate\"," << std::endl;
	
	fileOutput << std::setprecision(6) << "\t\t\"Solution\": " << val << "," << std::endl;
	
	int counter = 1;	
	fileOutput << "\t\t\"SolutionSet\": \"(";
	for(int i = 0; i < n; i++){
        if ( (sol[i] == 1) ){
        	if (counter != val){
            	fileOutput << i+1 << ", ";
            	++counter;
            }
            else
            	fileOutput << i+1 << ")\"," << std::endl;	
        }    
    }
   
    fileOutput << "\t\t\"BabNodes\": " << counterBaBNodes << std::endl;
    fileOutput << "\t}" << std::endl;
    
	fileOutput << "}" << std::endl;
	
	fileOutput.close();
}


void outputFile_after_start(char *argv[], int n, int m, int sol[], int val,  int counterBaBNodes, double timeTotal, double density, bool begin){
	
	std::string nameOfFile = argv[2];

	if (begin){
		nameOfFile.append("_start_info.txt");
		std::ofstream fileOutput(nameOfFile);


		 fileOutput << "{" << std::endl;


        	fileOutput << "\t\"GraphMD\": {" << std::endl;
        	fileOutput << "\t\t\"ObjectTypeString\": \"Graph\"," << std::endl;
        	fileOutput << "\t\t\"Vertices\": " << n << "," << std::endl;
        	fileOutput << "\t\t\"Edges\": " << m << "," << std::endl;
        	fileOutput << "\t\t\"Density\": " << std::setprecision(3) << density << "\n";
        	fileOutput << "\t}," << std::endl;


        	fileOutput << "\t\"ExecutionMD\": {" << std::endl;
        	fileOutput << "\t\t\"ExecutionTime\": " << std::setprecision(3) << timeTotal << ",\n";

        	fileOutput << "\t\t\"SolutionType\": \"Approximate\"," << std::endl;

        	fileOutput << std::setprecision(6) << "\t\t\"Solution\": " << val << "," << std::endl;

        	int counter = 1;
        	fileOutput << "\t\t\"SolutionSet\": \"(";
        	for(int i = 0; i < n; i++){
        		if ( (sol[i] == 1) ){
                		if (counter != val){
                			fileOutput << i+1 << ", ";
                			++counter;
            		}
            		else
                		fileOutput << i+1 << ")\"," << std::endl;
        		}
    		}	

    		fileOutput << "\t\t\"BabNodes\": " << counterBaBNodes << std::endl;
   		fileOutput << "\t}" << std::endl;

	        fileOutput << "}" << std::endl;
		
		fileOutput.close();

	}
	else {
		static int count = 1;
		nameOfFile.append("_temporary_solution_");
		nameOfFile.append(std::to_string(count));
		nameOfFile.append(".txt");

		std::ofstream fileOutput(nameOfFile);

		fileOutput << "{" << std::endl;

                fileOutput << "\t\"GraphMD\": {" << std::endl;
                fileOutput << "\t\t\"ObjectTypeString\": \"Graph\"," << std::endl;
                fileOutput << "\t\t\"Vertices\": " << n << "," << std::endl;
                fileOutput << "\t\t\"Edges\": " << m << "," << std::endl;
                fileOutput << "\t\t\"Density\": " << std::setprecision(3) << density << "\n";
                fileOutput << "\t}," << std::endl;


                fileOutput << "\t\"ExecutionMD\": {" << std::endl;
                fileOutput << "\t\t\"ExecutionTime\": " << std::setprecision(3) << timeTotal << ",\n";

                fileOutput << "\t\t\"SolutionType\": \"Approximate\"," << std::endl;

                fileOutput << std::setprecision(6) << "\t\t\"Solution\": " << val << "," << std::endl;

                int counter = 1;
                fileOutput << "\t\t\"SolutionSet\": \"(";
                for(int i = 0; i < n; i++){
                        if ( (sol[i] == 1) ){
                                if (counter != val){
                                        fileOutput << i+1 << ", ";
                                        ++counter;
                        }
                        else
                                fileOutput << i+1 << ")\"," << std::endl;
                        }
                }

                fileOutput << "\t\t\"BabNodes\": " << counterBaBNodes << std::endl;
                fileOutput << "\t}" << std::endl;

                fileOutput << "}" << std::endl;
		
		fileOutput.close();
		count++;
	}
}

