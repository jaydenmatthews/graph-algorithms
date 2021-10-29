// Centrality Measures API implementation


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "PQ.h"

/***** ReachableNode *****/
//Returns number of reachable nodes from v
static int ReachableNode(Graph g, Vertex v) {
	ShortestPaths sps = dijkstra(g, v);
	int i = 0, result = 0;
	
	
	while(i < sps.numNodes) {
		if (sps.dist[i] != INT_MAX) {
			result += 1;
		}
		i++;
	}
	
	freeShortestPaths(sps);
	
	return result;
}

/***** TotalSPDist *****/
//Returns total distance of shortest paths from v
static int TotalSPDist(Graph g, Vertex v) {
	ShortestPaths sps = dijkstra(g, v);
	int i = 0;
	double result = 0; 
	
	
	while(i < sps.numNodes) {
		if (sps.dist[i] != INT_MAX) {
			result += sps.dist[i];
		}
		i++;
	}
	
	freeShortestPaths(sps);
	
	return result;
}

NodeValues closenessCentrality(Graph g) {
	NodeValues nvs = {0};
	nvs.values = malloc(sizeof(double) * GraphNumVertices(g));
	nvs.numNodes = GraphNumVertices(g);
	
	double N = nvs.numNodes;
	
	for (int v = 0; v < nvs.numNodes; v++) {
	
		double n = ReachableNode(g, v);
		double sum = TotalSPDist(g, v);
		
		//Wasserman and Faust Formula
		double ans = ((n - 1)/(N - 1)) * ((n - 1)/sum);
	
		if (ans != ans) { //check for NaN
			ans = 0.0;
		}
		
		nvs.values[v] = ans;
	}
	
	
	return nvs;
}

/***** pathfinder *****/
//Recurses backwards from destination vertex to src, through all predecessors
//Finds the total number of shortest paths and shortest paths through v
static double pathfinder(ShortestPaths sps_src, Vertex dest, Vertex src, Vertex v, bool check, double *SPThroughV) {
	
	if (dest == src) { // we have reached the src vertex
		if (check) {
			// If we have reached intermediate v vertex 
			// then increment shortest paths through v
			*SPThroughV += 1;	
		}
		return 1;
	}
	
	if (dest == v) { // we have reached the intermediate "v" vertex
		check = true;
	}
	
	PredNode *curr = sps_src.pred[dest];
	double total = 0.0;
	
	// Recursively call pathinder for for the each predecssors node
	while (curr != NULL) { 
		total += pathfinder (sps_src, curr->v, src, v, check, SPThroughV);
		curr = curr->next;
	}
	
	return total;
}


NodeValues betweennessCentrality(Graph g) {
	NodeValues nvs = {0};
	nvs.numNodes = GraphNumVertices(g);
	nvs.values = malloc(sizeof(double) * GraphNumVertices(g));
	
	int n = 0, i = 0, j = 0;
	int nV = GraphNumVertices(g);
	bool passed = false;	
	double sptv = 0;
	double TotalNumSP;
	double btw = 0.0;
	
	while (n < nV) {//loops for each vertex
		i = 0;
		btw = 0.0;
		
		//These 2 loops find the shortest paths through v
		// for each pair of vertices
		while (i < nV) {
			j = 0;
			ShortestPaths sps_s = dijkstra(g, i);
			while (j < nV) {
				sptv = 0;
				passed = false;
				if (i != j && i != n && j != n) {	
					TotalNumSP = pathfinder(sps_s, j, i, n, passed, &sptv);
					if (TotalNumSP != 0.0) {
						btw += sptv / TotalNumSP;
					}	
				}	
				++j;
			}	
			freeShortestPaths(sps_s);
			++i;
		}
		
		nvs.values[n] = btw;
		++n;
	}	
		
	return nvs;
}

NodeValues betweennessCentralityNormalised(Graph g) {
	int n = GraphNumVertices(g);
	NodeValues nvs = betweennessCentrality(g); 

	
	double normaliser = (n-1)*(n-2);
	
	for (int i = 0; i < n; i++) {
		nvs.values[i] = nvs.values[i] / normaliser;
	}
	return nvs;
}

void showNodeValues(NodeValues nvs) {

}

void freeNodeValues(NodeValues nvs) {
	free(nvs.values);
	return;
}

