// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY DBL_MAX
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

static double AdjListWeightFinder(Graph g, Vertex v, Vertex w);
static void DistArrayInitialiser(int size, double dist[size][size]);
static void MinDistanceFinder(int size, double dist[size][size], Vertex *v, Vertex *w);
static int DendContainsVertex(Dendrogram root, Vertex v);
static Dendrogram RecursiveHAC(Graph g, int size, double dist[size][size], Dendrogram *dendA, int method);
static void LanceWilliamsDistance(Graph g, Dendrogram temp, Dendrogram new_cluster, double *min, int method);
static void ClusterCompare(Graph g, Vertex v, Dendrogram new_cluster, double *temp_min, int method);
static void DendContains(Dendrogram new_cluster, Dendrogram temp, int *contains);

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {
	int nV = GraphNumVertices(g);


	double dist[nV][nV];
	DistArrayInitialiser(nV, dist);
	
	// Create our initial distance array
	for (int i = 0; i < nV; i++) {
	 	for (int j = 0; j < nV; j++) {
	 		
	 		if (GraphIsAdjacent(g,i,j)) {
	 			double w1 = AdjListWeightFinder(g, i, j);
	 			double w2 = AdjListWeightFinder(g, j, i);
	 			
	 			double weight = max(w1, w2);
	 			
	 			dist[j][i] = 1 / weight;
	 			dist[i][j] = 1 / weight;
	 			
	 		} else {
	 			if (dist[i][j] == 0.0) {
	 				dist[i][j] = INFINITY;
	 			}
	 		}
	 		 		
	 	}
	}
	
	//create an array of size n dendA where each cell is a pointer to a dendogram nodes
	// set the vertex to the respective v number
	Dendrogram *dendA = malloc(sizeof(Dendrogram) * nV);
	for (int i = 0; i < nV; i++) {
		Dendrogram new_node = malloc(sizeof(DNode));
		new_node->vertex = i;
		new_node->left = NULL;
		new_node->right = NULL;
	 	dendA[i] = new_node;
	}
	
	Dendrogram final = RecursiveHAC(g, nV, dist, dendA, method);
	
	return final;
}

/***** RecursiveHAC *****/
//Recursively builds a dendrogram by freeing and creating a new dend array and dist array
//with each call of size n-1. Removes the need to create new nodes except root nodes of new clusters.
//Free the old dend array each time.
static Dendrogram RecursiveHAC(Graph g, int size, double dist[size][size], Dendrogram *dendA, int method) {
	int prev_size = size;
	int new_size = size - 1;
	
	
	if (new_size == 0) {
		Dendrogram final = dendA[0];
		free(dendA);
		return final;
	}
	// Once MinDistance found between two clusters. Create a new dendogram array
	// of size less by one.
	Vertex v, w;
	MinDistanceFinder(prev_size, dist, &v, &w);

	//The two clusters we will merge
	Dendrogram temp1 = dendA[v]; 
	Dendrogram temp2 = dendA[w];
	
	
	//New dendrogram array
	Dendrogram *newDendA = malloc(sizeof(Dendrogram) * new_size);
	
	//This is the newcluster root node;
	Dendrogram newcluster = malloc(sizeof(DNode));
	newcluster->vertex = -1;
	newcluster->left = temp1;
	newcluster->right = temp2;
	newDendA[new_size-1] = newcluster;
	

	int i = 0;
	int j = 0;
	while (i < prev_size) {
		int contains = 0;
		DendContains(newcluster, dendA[i], &contains);
		if (contains == 0) {
		 	newDendA[j] = dendA[i];
		 	++j;
		}
		
		++i;
	}
	
	free(dendA); //Free the pointers of the old dendA array
	
	//Create a distance array of size less by one, use the previous dist array
	double newDist[new_size][new_size];
	DistArrayInitialiser(new_size, newDist);
	// copying all values of pairs i and j where i and j are not elements of the 
	//new created array
	int posi = 0;
	for (int i = 0; i < prev_size; i++) {
		int posj = 0;
		if (!DendContainsVertex(newcluster, i)) {
			for (int j = 0; j < prev_size; j++) {
				if (!DendContainsVertex(newcluster, j)) {
					newDist[posi][posj] = dist[i][j];
					posj += 1;
				}
			}
			posi += 1;
		}
	}

	//Update dist Array using Lance-Williams Formula (Depends on Method)
	//New array will be of size n-1 x n-1 
	 for (int i = 0; i < new_size - 1; i++) {
	 	Dendrogram temp = newDendA[i];
	 	double min; 	
	 	
	 	if (method == 1) {
	 		min = -1.0;
	 	} else if (method == 2) {
	 		min = INFINITY;
	 	}	

	 	LanceWilliamsDistance(g, temp, newcluster, &min, method);
	 	
	 	if (method == 1 && min == -1.0) {
	 		newDist[i][new_size - 1] = INFINITY;
	 		newDist[new_size - 1][i] = INFINITY;
	 	} else if (method == 2 && min == INFINITY) {
	 		newDist[i][new_size - 1] = INFINITY;
	 		newDist[new_size - 1][i] = INFINITY;
	 	} else {
	 		newDist[i][new_size - 1] = 1 / min;
	 		newDist[new_size - 1][i] = 1 / min;

	 	}
	 }
	 
	 newDist[new_size - 1][new_size - 1] = INFINITY;
	 
	return RecursiveHAC(g, new_size, newDist, newDendA, method);
}
	

/**
 * Frees all memory associated with the given Dendrogram structure.
 */
void freeDendrogram(Dendrogram d) {
	if (d == NULL) {
		return;
	}
	
	freeDendrogram(d->left);
  	freeDendrogram(d->right);
    free(d);
}

/***** AdjListWeightFinder *****/
//Finds the weight of the edge from vertex v to w
static double AdjListWeightFinder(Graph g, Vertex v, Vertex w) {
	AdjList temp = GraphOutIncident(g, v);
	AdjList curr = temp;
	
	while (curr != NULL) {
		if (curr->v == w) {
			return (double)curr->weight;
		}
		curr = curr->next;
	}
	
	return 0.0;
}

/***** DistArrayInitialiser *****/
//Initialises the DistArray with all 0.0
static void DistArrayInitialiser(int size, double dist[size][size]) {
	
	for (int i = 0; i < size; i++) {
	 	for (int j = 0; j < size; j++) {
	 		dist[i][j] = 0.0;
	 	}
	 }
	 
}

/***** MinDistanceFinder *****/
//Finds the minimum distance between two clusters
static void MinDistanceFinder(int size, double dist[size][size], Vertex *v, Vertex *w) {
	double min = INFINITY;
	for (int i = 0; i < size; i++) {
	 	for (int j = i; j < size; j++) {
	 		if (dist[i][j] < min) {
	 			min = dist[i][j];
	 			*v = i;
	 			*w = j;
	 		}
	 	}
	 }
}

/***** DendContainsVertex *****/
//Predorder search to check whether a dendrogram contains a vertex
static int DendContainsVertex(Dendrogram root, Vertex v) {
	if (root == NULL) {
		return 0;
	}
	if (root->vertex == v) {
		return 1;
	}
	
	int res1 = DendContainsVertex(root->left, v);
	
	if (res1) {
		return 1;
	}

	int res2 = DendContainsVertex(root->right, v);
	
	return res2;
}

/***** LanceWilliamsDistance *****/
// Finds the distance between clusters based on the method required.
static void LanceWilliamsDistance(Graph g, Dendrogram temp, Dendrogram new_cluster, double *min, int method) {
	if (temp == NULL) {
		return;
	}
	if (temp->left == NULL && temp->right == NULL) { //is a leaf node
		double weight; 
		
		if (method == 1) {
	 		weight = -1.0;
	 	} else if (method == 2) {
	 		weight = INFINITY;
	 	}		
		//Compare this vertex with the new cluster vertices, finding the weight between them
		ClusterCompare(g, temp->vertex, new_cluster, &weight, method);
		if (method == 1) {
			if (weight > *min) { 
				*min = weight;
			}
		} else if (method == 2) {
			if (weight < *min) { 
				*min = weight;
			}
		}
		//compare the vertex of this node, to all the nodes in new cluster 
	}
	
	if (temp->left != NULL) {
		LanceWilliamsDistance(g, temp->left, new_cluster, min, method);
	}
	if (temp->right != NULL) {
		LanceWilliamsDistance(g, temp->right, new_cluster, min, method);
	}

}

/***** ClusterCompare *****/
// Compares the vertices of a cluster with a specific vertex, and returns 
// (depending on method) the distance.
static void ClusterCompare(Graph g, Vertex v, Dendrogram new_cluster, double *weight, int method) {
	if (new_cluster == NULL) {
		return;
	}
	if (new_cluster->left == NULL && new_cluster->right == NULL) { //is a leaf node
		//printf("Comparing Vertex %d to Vertex %d\n", v, new_cluster->vertex);
		//compare the v to this current leaf node;
		if (GraphIsAdjacent(g, v, new_cluster->vertex) || GraphIsAdjacent(g, new_cluster->vertex, v)) {
			double w1 = AdjListWeightFinder(g, v, new_cluster->vertex);
			double w2 = AdjListWeightFinder(g, new_cluster->vertex, v);
			
			if (method == 1) {
				double w3 = max(w1,w2);
				if (w3 > *weight) {
					*weight = w3;
				}
			} else if (method == 2) {
				if (w1 == 0.0) {
						w1 = INFINITY;
				}
				if (w2 == 0.0) {
						w2 = INFINITY;
				}
				
				double w3 = min(w1,w2);
				
				if (w3 < *weight) {
					*weight = w3;
				}
			}
		} 
	}
	
	if (new_cluster->left != NULL) {
		ClusterCompare(g, v, new_cluster->left, weight, method);
	}
	if (new_cluster->right != NULL) {
		ClusterCompare(g, v, new_cluster->right, weight, method);
	}
}

/***** DendContains *****/
// Checks whether a cluster contains any vertices from a new formed cluster.
static void DendContains(Dendrogram new_cluster, Dendrogram temp, int *contains) {
	if (temp == NULL) {
		return;
	}
	if (temp->left == NULL && temp->right == NULL) { //is a leaf node
		//compare the vertex of this node, to all the nodes in new cluster
		if (DendContainsVertex(new_cluster, temp->vertex)) {
			*contains = 1;
		}
	}
	
	if (temp->left != NULL) {
		DendContains(new_cluster , temp->left, contains);
	}
	if (temp->right != NULL) {
		DendContains(new_cluster, temp->right, contains);
	}
}
