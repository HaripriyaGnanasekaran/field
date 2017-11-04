#ifndef SF_LINKNODExH
#define SF_LINKNODExH

#define MAXLINKS 3

#include<stdio.h>
#include<stdlib.h>
#include <fenk.h>

/*! \brief A Link class for branched polymers.

A link goes from node1 to node2; seg1 belongs both to link.node1 
and this link, seg1 belongs both to link.node1 and this link.

All int functions in this class return -1 on error
*/
class SF_Link {
  public:
  	//! Constructs a link wiht all data initialised to -1.
	SF_Link(void);
	//! Constructs a link with the data determined by the input parameters.
	/*! Constructs a link with the data determined by the input parameters.
	\param _length is the length of the link.
	\param _node[2] provides the IDs of the start and the end node.
	\param _seg[2] provides the IDs of the corresponding segments that are parts of the nodes
	\sa SF_Link(void);
	*/
	SF_Link(int _length, int _node[2], int _seg[2]); 
	//! A destructor.
	~SF_Link();
	//! Returns link length.
	int length(); 
	//! Returns the corresponding node number.
	/*! Returns the corresponding node number.
	\param number can be 1 for the start node or 2 for the end node
	*/
	int node(int number); 
	//! Returns the corresponding segment number by which the link is connected to the node number.
	/*!  \param number can be 1 for the start segment or 2 for the end segment
	*/
	int seg(int number); 
	//! Dumps the link to stdout.
	/*! Obsolete! Made just for testing purposes. 
	Dumps the link in format as stored in Golitah: 
	seg2 seg1 length node1 node2
	*/
	void dump(void);
  private: 
  	int linkseg[2]; /*! IDs ot the segments at the ends of the link */
	int linknode[2]; /*! IDs ot the nodes to which the link is connected */
	int linklength; /*! The length of the link */
};

/*! \brief A Node class for branched polymers.

Each Node contains the numbers of links connected to it as well as the numbers 
of segments by which they are connected.
The sgements are considered all to be parts of the node.

All int functions in this class return -1 on error
*/
class SF_Node {
  public:
  	//! Constructs a node wihtout any data
	SF_Node(void);
  	//! Constructs a node with the data determined by the input parameters.
	/*! Constructs a node with the data determined by the input parameters.a
	IMPORATNT: assumes that _link and _seg are allocated to the size of _n_links
	\param _n_links is the length of the number of links of the node
	\param _link[] provides the IDs of the links connected to the node
	\param _seg[] provides the IDs of the corresponding segments by which the links are connected to the node
	\sa SF_Node(void);
	*/
	SF_Node(int _n_links, Array<int> _link, Array<int> _seg); 

	//! A destructor
	~SF_Node();
	//! Returns the number of links of the node
	int n_links(void);
	//! Returns the corresponding link ID of the link link_num
	/*! \param link_num is the number of the link as stored inside the node, i.e. between 1 and n_links */
	int link(int link_num);
	//! Returns the corresponding segment number of the link link_num
        /*! \param link_num is the number of the link as stored inside the node, i.e. between 1 and n_links */
	int seg(int link_num);
	//! Dumps the Node to stdout. 
	/*! \param n_dump defines the number of columns that should
	bed dumped, the columns that contain no data are dumped as "-"
	*/
	void dump(int n_dump);
  private:
  	int *nodelink; /*! IDs ot the links conected to the node */
	// FIXME *nodeseg is probably not going to be used anymore
	int *nodeseg; /*! numberes ot the segments conected to the node */
	int node_n_links; /*! number of links connected to the node */
};

/*
typedef struct SF_LinkNodeTree{
	Array<SF_Link> LinkList;
	Array<SF_Node> NodeList;

private:
	void SF_LinkNodeTree:Dump(const SF_LinkNodeTree* tree) const;
};
*/
//! The data object that stores the topology of the branched polymer
/*! NOTE: it is the responsibility of the one who uses the data structure that
proper space is allocated for the link[] and node[] arrays

*/

class SF_LinkNodeTree{
public:
	Array<SF_Link> LinkList; //! the array of SF_Link's
	Array<SF_Node> NodeList; //! array of SF_Node's
	int n_links; //! the number of links
	//! dump the tree to stdout
	void Dump(void);
};

#endif
