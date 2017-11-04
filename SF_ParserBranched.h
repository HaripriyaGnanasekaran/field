#ifndef SF_PARSERBRANCHEDxH
#define SF_PARSERBRANCHEDxH
#include <fenk.h>
#include "SF_LinkNode.h"
#include "SF_MolStructure.h"


//forward declaration
class SF_MolStructure;

//!	A helper class, reads branched compositions, builds LinkNodeTree and SegBlockQ.
class SF_ParserBranched {
	
public:

//!	Constructs the parser, processes the composition, sets values of SF_MolStructure members.
/*! \param composition_ is the composition to be parsed
	\param MolStructure_ is a pointer to a SF_MolStructure object to be initialised. 
		SF_ParserBranched directly sets SF_MolStructure::LinkNodeTree,  SF_MolStructure::SegBlockQ,
		SF_MolStructure::repeatQ and SF_MolStructure::numGenerations.

	What it does, in words:
		1)	Breaks the composition to "blocks". Blocks are substrings of the initial composition separated by 
			square brackets. Does not distinguish '[' and ']'.
		2)	Assigns ranks to the blocks in the order of their appearance in the initial composition string (left to right). 
				The first block has rank 1. 
				A block is assigned with the rank of the preceding block +1 if they are separated by '['.
				A block is assigned with the rank of the preceding block -1 if they are separated by ']'.
				A block following '][' (multiple branches) is assigned with	the rank of the preceding block, 
				the ranks of the preceding blocks are then increased by 1 in reverse order (right to left) 
				untill a block with a rank which is lower than that of the block which triggered the elevation
				of ranks is found. 
				Combinations '[[', ']]', '[]' result in a syntax error.
		3)	Constructs LinkNodeTree by connecting the first segment of each block to the last segment of a preceding
			block with a rank which does not exceed its own rank.
		4)	Sets 
				SF_MolStructure::LinkNodeTree - Pointer to the constructed tree 
				SF_MolStructure::SegBlockQ	- New SegBlock for each branch
				SF_MolStructure::repeatQ - 1
				SF_MolStructure::numGenerations - number of blocks

*/
	SF_ParserBranched(Text composition_, SF_MolStructure* MolStructure_, Input* MyInput_);
	~SF_ParserBranched();  

private:
	SF_MolStructure* Struct;
	Text Composition;
	int numBrackets;
	int numBlocks;
	Array<int> BracketPos;  
	Array<int> BracketType;
	Array<int> BlockRanks;
	Array<Text> Blocks;
	Input* MyInput;

	///Expands brackets containing branched fragments, eg ((A)5[(B)1](C)5)2 -> (A)5[(B)1](C)5(A)5[(B)1](C)5
	void ExpandBrackets();
	///Reads square brackets positions and types, stores them in BracketPos and BracketType (types: '['=1; ']'=-1)
	void ParseBrackets();
	///Calls ParseBrackets(), splits Compositon into linear fragments, stores them in Blocks, assigns BlockRanks.
	void ParseBlocks();
	///Calls FillLinks() and FillNodes()
	void FillLinksAndNodes();
	///Fills Links in the LinkNodeTree of a MolStructure object
	void FillLinks();
	///Fills Nodes in the LinkNodeTree of a MolStructure object
	void FillNodes();
	///Fills SegBlockQ of a MolStructure object with blocks corresponding to linear fragments of the composition
	void FillSegBlockQ();
	///Returns the position of the first unmatched '(' to the left from the position pos_ in Composition
	int FindOpen(int pos_); 
	///Returns the position of the first unmatched ')' to the right from the position pos_ in Composition
	int FindClose(int pos_); 
	///Returns text containing only digits to the right from the position pos_ in Composition
	Text ReadRepeat(int pos_); 
#if DEBUG
	void DumpBrackets();
	void DumpBlocks();
	void DumpLinkNodeTree();
#endif
};


#endif
