#ifndef HIV_HEADER
#define HIV_HEADER
#define _GNU_SOURCE

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <glib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include <unistd.h>
#include <regex.h>
#include <getopt.h>


/* Bifurcating tree */
struct NodeSampledTree {
  char * name;
  double branchLength;
  struct NodeSampledTree* left;  /* NULL if no children*/
  struct NodeSampledTree* right; /* only viruses with no children are replicating */
  struct NodeSampledTree* parent; /* only viruses with no children are replicating */
  int label;

} ;

bool FloatEquals(double a, double b, double threshold);
struct NodeSampledTree* newNodeSampledTree(double branchLength, char * name, struct NodeSampledTree* parent);
int makeTree(struct NodeSampledTree* node_ptr, char* treeString, int currChar);
void readLatent (struct NodeSampledTree* node, FILE* latentFile);
int countTips(struct NodeSampledTree* node, int count);
void printTree(struct NodeSampledTree* root_ptr, char* fileName);
char* NewickStringBranchLength(struct NodeSampledTree* node, char* treeString, double previousTime);
void clearTree(struct NodeSampledTree* node_ptr);

void clearMemory(struct NodeSampledTree* root_ptr, int numTips, int* nodeLatentstate, int* nodeNames, double* nodeSampleTime, char** alignment, char* treeString, gsl_rng *r);
char** allocateAlignmentMem(int numTips, int numBases);

#endif
