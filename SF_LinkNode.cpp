#include "SF_LinkNode.h"

SF_Link::SF_Link(void) { }

SF_Link::SF_Link(int _length, int _node[2], int _seg[2]) {
	linklength=_length;
	linkseg[0]=_seg[0];	
	linkseg[1]=_seg[1];	
	linknode[0]=_node[0];	
	linknode[1]=_node[1];
}

SF_Link::~SF_Link() { }

int SF_Link::node(int number) { 
  if(number>0&&number<3) return linknode[number-1];
  else return -1;
}

/* old implementation of seg()
int SF_Link::seg(int number) { 
  if(number>0&&number<3) return linkseg[number-1];
  else return -1;
} */

int SF_Link::seg(int number) { 
	if(number==linknode[0]) return 1;
	else if (number==linknode[1]) return linklength;
	else return -1; // error code
}

int SF_Link::length(void) { return linklength;}

void SF_Link::dump(void) {
  printf("%2d %2d %2d %2d %2d\n",linkseg[1],linkseg[0],linklength,linknode[0],linknode[1]);
}


/* 
******************************
here the node begins 
******************************** 
*/

SF_Node::SF_Node(void) {
  node_n_links=-1; // just construct the bare object, do not allocate anything
}

SF_Node::SF_Node(int _n_links, Array<int> _link, Array<int> _seg) {
  int i;
  node_n_links=_n_links;
  nodeseg=(int*)malloc(node_n_links*sizeof(int));
  nodelink=(int*)malloc(node_n_links*sizeof(int));
  for(i=0;i<node_n_links;i++) {
    nodeseg[i]=_seg[i];
    nodelink[i]=_link[i]; 
  }
}

SF_Node::~SF_Node() { }

int SF_Node::n_links() { return node_n_links; }

int SF_Node::link(int number) { 
  if(number>0&&number<(node_n_links+1)) return nodelink[number-1];
  return -1;
}

int SF_Node::seg(int number) { 
  if(number>0&&number<(node_n_links+1)) return nodeseg[number-1];
  return -1;
}

void SF_Node::dump(int n_dump) {
  int i;
  if(node_n_links<0) { 
    fprintf(stderr,"Error: trying to dump a non-allocated SF_Node\n");
    return;
  }
/*  for(i=n_dump;i>0;i--) {
    if(i>node_n_links) printf("-  ");
    else printf("%2d ",nodeseg[i]); 
  }
  */
  printf("  %2d   ", node_n_links);
  for(i=0;i<n_dump;i++) { 
   if(i<node_n_links) printf(" %2d",nodelink[i]); 
   else printf("  -");
  }
//  for(i=0;i<MAXLINKS;i++) printf(" %d",nodelink[i+1]); 
  printf("\n"); fflush(stdout);
}

void SF_LinkNodeTree::Dump(void) {
	// dump the nodes
	int i; 
	int numlinks=n_links;
	printf("Links:\n");
	for (i=0;i<numlinks;) LinkList[++i].dump();
	numlinks++;
	printf("Nodes:\n");
	printf ("n_links, links ->\n");
	for (i=0;i<numlinks;) { 
		NodeList[i+1].dump(numlinks);
		i++;
	}
	fflush(stdout);
	return;	
}
