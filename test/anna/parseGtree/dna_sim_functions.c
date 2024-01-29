#include "dna_sim.h"
#include "sasprintf.c"
#include <ctype.h>
#include <math.h>

/* Checks if two doubles are within a threshold value of each other */
bool FloatEquals(double a, double b, double threshold) {
  return fabs(a-b) < threshold;
}

/* Checks the input format for the parameters in the rate matrix */
void checkModelInput(int regMatch) {
  if (regMatch != 0) {
    fprintf(stderr, "Model parameters are not in the correct format.\nThe correct format is statFreq_A:statFreq_C:statFreq_G:statFreq-rateParameter_1:rateParameter_2:rateParameter_3:rateParameter_4:rateParameter_5:rateParameter_6\nNumbers must be in decimal format with a digit after the decimal. \nExiting the program. \n");
    exit(1);
  }
}

/* Create a new new node. Allocates memory for the node and the latency states
structure. Assignment values to the variables in the node structure */
struct NodeSampledTree* newNodeSampledTree(double branchLength, char * name, struct NodeSampledTree * parent) {
  /* Allocates memory for a new node */
  struct NodeSampledTree* node_ptr;
  node_ptr  = (struct NodeSampledTree*)malloc(sizeof(struct NodeSampledTree));

  if (node_ptr == NULL) {
    fprintf(stderr, "Insufficient Memory\n");
    exit(1);
  }

  node_ptr->branchLength = branchLength;
  node_ptr->name = name;

  /* Sets pointers to daughter nodes */
  node_ptr->left = NULL;
  node_ptr->right = NULL;
  node_ptr->parent = parent;
  node_ptr->label = 0;
  return(node_ptr);
}

// Read in files

/* The function parses a string into a tree structure. Every time there's an open
parenthesis, two new nodes are created and the function is called twice. The
function keeps track of which character in the treeString (currChar) it is looking
at. */
int makeTree(struct NodeSampledTree* node_ptr, char* treeString, int currChar) {
	char* nodeName;
  /* Note: You should never actually reach the end of string character. The last
  time you return is after the last closed parenthesis */
  while( treeString[currChar] != '\0') {

    /* If the current character is an open parenthesis, make new node, call
    makeTree recursively. currChar is updated along the way */
    if (treeString[currChar] == '(') {
      node_ptr->left = newNodeSampledTree(-1, "", node_ptr);
      currChar = makeTree(node_ptr->left, treeString, currChar + 1);
      node_ptr->right = newNodeSampledTree(-1, "", node_ptr);
      currChar = makeTree(node_ptr->right, treeString, currChar);


    } else {
      /* First, find the node name. Then, find the branch length. Return the
      current character */
      char* ptr; /* Needed for the strtod function */
      char* residualTreeString = malloc(strlen(treeString) + 1); /* Used so the original string isn't modified */
      strcpy(residualTreeString, treeString + currChar);

      regex_t regDecimal, regInt;
      int regCompiled, regMatch;
      /* Regex of decimal number */
      regCompiled = regcomp(&regDecimal, "^([0-9]+)((\\.)([0-9]+))$" , REG_EXTENDED);
      if (regCompiled == 1) {
        fprintf(stderr, "Regular expression did not compile.\n");
        exit(1);
      }
      /* Regex of integer number */
      regCompiled = regcomp(&regInt, "^([0-9,A-Z,a-z,^]+)$" , REG_EXTENDED);
      if (regCompiled == 1) {
        fprintf(stderr, "Regular expression did not compile.\n");
        exit(1);
      }

      /* Finds the nodeName by looking for the next ":". Convert it into an
      integer, save it in the node structure. Update currChar. */
      if (residualTreeString[0] !=  ':') {
      	nodeName = strtok(residualTreeString, ":"); /*Note this function modified residualTreeString */

      	regMatch = regexec(&regInt, nodeName, 0, NULL, 0);
      	if (regMatch != 0) {
      		fprintf(stderr, "Problem reading in tree file. Regular expression does not match a integers.\n");
      		exit(1);
      	}
      	node_ptr->name = strdup(nodeName);
      } else {
      		node_ptr->name = strdup("");
		nodeName = "";
      }

      
      //node_ptr->name = strdup(nodeName);
      currChar = currChar + strlen(nodeName) + 1;

      residualTreeString = strcpy(residualTreeString, treeString + currChar);

      /* Finds the branch length, converts it to a double, saves it in the node
      structure */
      char* branchLength = strtok(residualTreeString, ",);");

      regMatch = regexec(&regDecimal, branchLength, 0, NULL, 0);
      if (regMatch != 0) {

	      printf("%s %s\n", nodeName, branchLength);
        fprintf(stderr, "Problem reading in tree file. Regular expression does not match a decimal number.\n");
        exit(1);
      }
      node_ptr->branchLength = strtod(branchLength, &ptr);
      currChar = currChar + strlen(branchLength) + 1 ;

      free(residualTreeString);
      regfree(&regDecimal);
      regfree(&regInt);

      /* Returns the updated current character */
      return(currChar);
    }
  }
  return(currChar);
}


/* Counts the number of tips on the tree structure.*/
int countTips(struct NodeSampledTree* node, int count) {
  if (node->left == NULL && node->right == NULL) {
    return(count + 1);

  } else {
    count = countTips(node->left, count);
    count = countTips(node->right, count);
    return count;
  }
}



/* This function prints the tree structure in newick format.
This can be used to check that the file is identical to the file that was
read in. A different function is required than the function originally used to
print the tree since the birth age is actually a branch length when it is read
n.*/
void printTree(struct NodeSampledTree* root_ptr, char* fileName) {
  char *treeStringCheck = strdup("");
  treeStringCheck = NewickStringBranchLength(root_ptr, treeStringCheck, 0);

  FILE *treeFile;
  treeFile = fopen(fileName, "w");
  fputs(treeStringCheck, treeFile);
  fputs("\n", treeFile);
  fclose(treeFile);
  free(treeStringCheck);
}

/* Used to print the tree to test that it matches the tree that was read in. See
printTree for more detail. */
char* NewickStringBranchLength(struct NodeSampledTree* node, char* treeString, double previousTime) {

  if (node->left == NULL && node->right == NULL) {
    Sasprintf(treeString, "%s%s:%f", treeString, node->name, node->branchLength);
    return (treeString);

  } else {
    Sasprintf(treeString, "%s(", treeString);
    treeString = NewickStringBranchLength(node->left, treeString, node->branchLength);
    Sasprintf(treeString, "%s,", treeString);
    treeString = NewickStringBranchLength(node->right, treeString, node->branchLength);
    Sasprintf(treeString, "%s)%s:%f", treeString, node->name, node->branchLength);

  }
  return (treeString);
}

/* Clears memory for the tree. Stem must be cleared separately */
void clearTree(struct NodeSampledTree* node_ptr) {
  if (node_ptr->left != NULL) {
    clearTree(node_ptr->left);
    clearTree(node_ptr->right);
  }

  free(node_ptr->name);
  /* Clear node */
  free(node_ptr);
  return;
}

/* Frees some of the memory in the program. (There's still a leak)*/
void clearMemory(struct NodeSampledTree* root_ptr, int numTips, int* nodeLatentstate, int* nodeNames, double* nodeSampleTime, char** alignment, char* treeString, gsl_rng *r) {

  clearTree(root_ptr);

  free(treeString);

  return;
}

