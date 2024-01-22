#include <stdio.h>
#include <stdlib.h>
#include "dna_sim.h"
#include "string_utilities.h"
double treeLength(struct NodeSampledTree* node, char * oldestSample);
double findMRCA (struct NodeSampledTree * node, char * name);
void labelPath (struct NodeSampledTree * node, char * name);

void labelPath (struct NodeSampledTree * node, char * name) {

	if(node->left == NULL) {
		if (strcmp(name, node->name) == 0) {
			
			for (struct NodeSampledTree * pathNode = node; pathNode->parent; pathNode = pathNode->parent) 
				pathNode->label = 1;

		}
		return;
	} else {
		labelPath(node->left, name);
		labelPath(node->right, name);
	}
}


double findMRCA (struct NodeSampledTree * node, char * name) {

	double length = 0;
	if(node->left == NULL) {
		if (strcmp(name, node->name) == 0) {
			for (struct NodeSampledTree * pathNode = node; pathNode->parent ; pathNode = pathNode->parent) {
				if (pathNode->label) 
					break;

				length +=pathNode->branchLength;
			}

		}
		return length;
	} else {
		return findMRCA(node->left, name) + findMRCA(node->right, name);
	}
}

double findCoalescentTime(struct NodeSampledTree * root, char * youngSample, char * oldSample) {

	labelPath(root, oldSample);
	return findMRCA(root, youngSample);
}

double treeLength(struct NodeSampledTree* tip, char * oldestSample) {
	
	double length = 0;
	if(tip->left == NULL) {
		if (strcmp(oldestSample, tip->name) == 0) {
		
			for (struct NodeSampledTree * node = tip; node->parent; node = node->parent) 
				length += node->branchLength;
			
			return length;
	       	} else 
			return  0;
	
       	} else 
		return (treeLength(tip->left, oldestSample) +treeLength(tip->right, oldestSample));
	
	
}

void clearLabel(struct NodeSampledTree * node) {

	node->label = 0;
	if(node->left) {
		clearLabel(node->left);
		clearLabel(node->right);
	}
	return;
}

int main (int argc, char **argv) {
  char* treeFile = argv[1];
  char * nameLongestBranch = argv[2];
  char * otherBranch = argv[3];
  char *treeString = string_from_file(treeFile);
  char * nextTree; 
  
  
  char * branch1 = argv[2];
  char * branch2 = argv[3];

  /*Make a tree */
	nextTree = strtok(treeString, "\n") ;
struct NodeSampledTree* root_ptr;	

  while (nextTree) {
  	root_ptr = newNodeSampledTree(0, "", NULL);
  	int endChar = makeTree(root_ptr, nextTree, 0);

	for(int i = 2; i < argc; i += 2) {
		branch1 = argv[i];
		branch2 = argv[i+1];

  		printf("%s_%s %f\n", branch1, branch2, findCoalescentTime(root_ptr, branch1, branch2));
		clearLabel(root_ptr);
	}
  	//printf("%f\n", treeLength(root_ptr, nameLongestBranch));

  	clearTree(root_ptr); 
	nextTree = strtok(nextTree + strlen(nextTree)+1, "\n") ;
	
  }

}
