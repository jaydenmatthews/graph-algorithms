// Dijkstra API implementation

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"


/***** deleteList *****/
// Recursive freeing of a PredNode list
static void deleteList(PredNode *pred) {
	if (pred == NULL) {
		return;
	} 
	deleteList(pred->next);
	free(pred);
}

/***** insertPredNode *****/
// Insertion of a prednode into a list considering empty and non-empty cases.
static void insertPredNode(Vertex v, Vertex u, PredNode **pred, int dist_same) {	
	PredNode *new_node = malloc(sizeof(PredNode)); 
	new_node->v = u;
	new_node->next = NULL;
	
	if (dist_same == 0) {
		deleteList(pred[v]);
		pred[v] = NULL;
		pred[v] = new_node;
		return;	
	} else if (dist_same == 1 ) {	
		PredNode *curr = pred[v];
		
		while (curr->next != NULL) {
			curr = curr->next;
		}
		
		curr->next = new_node;
		return;	
	}
	return;	

} 


ShortestPaths dijkstra(Graph g, Vertex src) {
	ShortestPaths sps = {0};
	int nV = GraphNumVertices(g);
	
	sps.src = src;		//set source vertex
	sps.numNodes = nV; //set number of vertices
	sps.dist = malloc(sizeof(int) * nV);
	sps.pred = malloc(sizeof(PredNode *) * nV);
	
	
	for (int i = 0; i < nV; i++) {
		sps.dist[i] = INT_MAX;
		sps.pred[i] = NULL;
	}
	
	sps.dist[src] = 0;
	
	PQ pq = PQNew(); //New Priority Queue
	
	PQInsert(pq, src, 0); // Insert source vertex into pq as distance 0
	
	while (!PQIsEmpty(pq)) { // While priority queue isnt empty
		int u = PQDequeue(pq); 
		
		//Get list of adjacent vertices
		AdjList weightlist = GraphOutIncident(g, u); 
		AdjList curr = weightlist;
		
		
		while (curr != NULL) {
			int v = curr->v;
			int weight = curr->weight;
			
			//Case where distance less
			if (sps.dist[v] > sps.dist[u] + weight){ 
				sps.dist[v] = sps.dist[u] + weight;
				PQInsert(pq, v, sps.dist[v]);
				insertPredNode(v, u, sps.pred, 0);
				
			//Case where distance is same	
			} else if (sps.dist[v] == sps.dist[u] + weight) { 
				insertPredNode(v, u, sps.pred, 1);
			}
			curr = curr->next;
			
		}
	}
	PQFree(pq);
	return sps;
}

void showShortestPaths(ShortestPaths sps) {

}


void freeShortestPaths(ShortestPaths sps) {
	//free sps.dist
	free(sps.dist);
	
	int i = 0;
	
	//delete each predecessor list
	while (i < sps.numNodes) { 
		deleteList(sps.pred[i]);
		i++;
	}
	
	free(sps.pred);
}


 

