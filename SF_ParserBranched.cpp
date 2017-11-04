#include "SF_ParserBranched.h"

SF_ParserBranched::SF_ParserBranched(Text composition_, SF_MolStructure* MolStructure_, Input* MyInput_) {
	Struct = MolStructure_;
	MyInput = MyInput_;
	Composition = Copy(composition_);
	char test;
	ExpandBrackets();
	#if DEBUG
		Message(debug, "SF_ParserBranched expanded composition: " + Composition);
	#endif
	Composition.Setpos(1);
	test=Composition.Getchar();
	if (test!='[') {
		Composition = "[" + Composition + "]";
	}
	ParseBlocks();
	FillSegBlockQ();
	FillLinksAndNodes();
	Struct->numGenerations = numBlocks;
//	cout << "numGenerations =" << Struct->numGenerations << endl;
	#if DEBUG
		DumpLinkNodeTree();
	#endif
}

SF_ParserBranched::~SF_ParserBranched() {

}
void
SF_ParserBranched::FillSegBlockQ(){
	int blnum;
	Struct->SegBlockQ.Dim(1,numBlocks);
	Struct->repeatQ.Dim(1,numBlocks);
	for (blnum=1; blnum<=numBlocks; blnum++) {
		Struct->SegBlockQ[blnum] = new SF_SegmentBlock(Blocks[blnum], Struct->molName, Struct->SegQ, MyInput, "", "");
		Struct->repeatQ[blnum] = 1;
	}
	if (Struct->SegBlockQ[1]->GetMaxLength() == 1)
		Message(fatal, MyInput, "Error in composition '" + Struct->composition + "'\nThe first linear fragment contains only one monomer.");
}

void
SF_ParserBranched::FillLinks(){
	int blnum, i, length;
	int nodes[2] = {1,2}, dummy[2] = {0,0};
	Struct->LinkNodeTree->LinkList.Dim(1,numBlocks);
	Struct->LinkNodeTree->n_links = numBlocks;
	length = Struct->SegBlockQ[1]->GetMaxLength();
	Struct->LinkNodeTree->LinkList[1] = SF_Link(length,nodes,dummy);
	for (blnum=2; blnum<=numBlocks; blnum++){
		length = Struct->SegBlockQ[blnum]->GetMaxLength()+1;
		nodes[1] += 1;
		i = blnum - 1;
		while ((BlockRanks[i] > BlockRanks[blnum]) && (i>0)) i--;
		nodes[0] = Struct->LinkNodeTree->LinkList[i].node(2);
		Struct->LinkNodeTree->LinkList[blnum] = SF_Link(length,nodes,dummy);
	}
}

void
SF_ParserBranched::FillNodes(){

	int blnum, i, j;
	Array<int> n_links(numBlocks+1), dumbo(2);
	Struct->LinkNodeTree->NodeList.Dim(1,numBlocks+1);
	for (i=0; i<=numBlocks; i++) n_links[i] = 0;
	for (blnum=1; blnum<=numBlocks; blnum++){ //count links connected to each node
		n_links[Struct->LinkNodeTree->LinkList[blnum].node(1)-1] += 1;
		n_links[Struct->LinkNodeTree->LinkList[blnum].node(2)-1] += 1;
	}
	Array<int> links(1);
	links[0] = 1;	//parent link
	Struct->LinkNodeTree->NodeList[1] = SF_Node(1, links, dumbo);
	for (i=1; i<=numBlocks; i++){
		Array<int> links(n_links[i]);
		j=1;
		links[0] = i;	//parent link
		for (blnum=i+1;((blnum<=numBlocks)&&(j<n_links[i]));blnum++)
			if (Struct->LinkNodeTree->LinkList[blnum].node(1) == i+1) links[j++] = blnum;
		dumbo.Dim(n_links[i]);
		Struct->LinkNodeTree->NodeList[i+1] = SF_Node(n_links[i], links, dumbo);
	}
}
void
SF_ParserBranched::FillLinksAndNodes(){
	Struct->LinkNodeTree = new SF_LinkNodeTree;
	FillLinks();
	FillNodes();
}


void
SF_ParserBranched::ExpandBrackets(){
	int i,bra,ket,count;
	int repeat, repeatlength;
	Text s, s_repeat, temp, t;
	t = Composition;
	char test;
	t.Setpos(1);
	count = 0;
	while (t.More()) {
		test=t.Getchar();
		if (test == '(') count++;
		if (test == ')') count--;
		if ((test == '[') && (count>0)) {
			bra = FindOpen(t.Pos()-1);
			ket = FindClose(t.Pos()-1);
			s = t.Sub(bra+1, ket-bra-1);
			s_repeat = ReadRepeat(ket);
			repeatlength = s_repeat.Length();
			repeat = repeatlength==0 ? 1 : s_repeat.Getint();
			temp = ( bra>1 ? t.Sub(1, bra-1) : "");
			for (i=1; i<=repeat; i++) temp = temp + s;
			temp = temp + (ket+repeatlength+1 < t.Length() ? t.Sub(ket+repeatlength+1, t.Length()-ket-repeatlength) : "");
			Composition = Copy(temp);
			ExpandBrackets();
			break;
		}
	}

}
int
SF_ParserBranched::FindOpen(int pos_){
	int pos = pos_, count=-1;
	Text t;
	char test;
	t = Composition;
	t.Setpos(pos);
	while ((count<0) && (pos>1)) {
		pos--;
		t.Setpos(pos);
		test=t.Getchar();
		if (test=='(') count++;
		if (test==')') count--;
	}
	if (count==0) return pos;
	Message(fatal, MyInput, "Parser has problems expanding brackets in composition '" + Composition + "'");
	return 0; //never get here

}
int
SF_ParserBranched::FindClose(int pos_){
	int pos = pos_, count=1;
	Text t;
	char test;
	t = Composition;
	t.Setpos(pos);
	while ((count>0) && (pos<t.Length())) {
		pos++;
		t.Setpos(pos);
		test=t.Getchar();
		if (test=='(') count++;
		if (test==')') count--;
	}
	if (count==0) return pos;
	Message(fatal, MyInput, "Parser has problems expanding brackets in composition '" + Composition + "'");
	return 0; //never get here
}
Text
SF_ParserBranched::ReadRepeat(int pos_){
	char test;
	Text t = Blanks(Composition.Length());
	Composition.Setpos(pos_);
	test=Composition.Getchar();
	while (Composition.More()) {
		test=Composition.Getchar();
		if (isdigit(test)) t.Putchar(test);
		else break;
	}
	return t.Strip();
}

void
SF_ParserBranched::ParseBrackets() {
	char test;
	int brnum = Composition.Length();
	Array<int> brpos(brnum), brtype(brnum);	//these local arrays overestimate the number of brackets
	int diff = 0;
	brnum = 0;
	Composition.Setpos(1);
	while (Composition.More()) {
		test=Composition.Getchar();
		if (test == '[') {
			diff++; brnum++;
			brtype[brnum] = 1;
			brpos[brnum] = Composition.Pos() -1;
		}
		if (test == ']') {
			diff--;	brnum++;
			brtype[brnum] = -1;
			brpos[brnum] = Composition.Pos() -1;
		}
	}
	if (diff != 0) {
		Message(fatal, MyInput, "Error in composition '" + Struct->composition + "'\nSquare brackets do not match");
	}
	numBrackets = brnum;
	BracketPos.Dim(1, numBrackets);
	BracketType.Dim(1,numBrackets);
	for (brnum=1; brnum<=numBrackets; brnum++){
		BracketPos[brnum] = brpos[brnum];
		BracketType[brnum] = brtype[brnum];
	}
}

void
SF_ParserBranched::ParseBlocks(){
	int k, brnum, blnum, bllength;
	ParseBrackets();
	blnum = numBrackets;
	Array<int> blranks(blnum);		//these local arrays overestimate the number of blocks
	Array<Text> blocks(blnum);		//
	Text s;
	blnum = 0;
	blranks[0] = 0;
	Boolean multibranch = false;
	for (brnum=1; brnum < numBrackets; brnum++) {
		bllength = BracketPos[brnum+1]-BracketPos[brnum]-1; //block length excluding brackets
		if (bllength > 0) {
			blnum++;
			blocks[blnum] = Composition.Sub(BracketPos[brnum]+1,bllength);
			if (!multibranch) {
				blranks[blnum] = blranks[blnum-1] + BracketType[brnum];
			}
			else {
				blranks[blnum] = blranks[blnum-1];
				k = blnum - 1;
				while (blranks[k]>=blranks[blnum]){
					blranks[k] += 1;
					k--;
				}
				multibranch = false;
			}
		}
		else {
			s = Copy(Composition.Sub(BracketPos[brnum],2));
			if (*s == *Copy("][")) multibranch = true;
			else Message(fatal, MyInput, "Unnecessary square brackets in composition : " + Struct->composition );
		}

	}
	numBlocks = blnum;
	BlockRanks.Dim(1, numBlocks);
	Blocks.Dim(1, numBlocks);
	for (blnum=1; blnum<=numBlocks; blnum++){
		BlockRanks[blnum] = blranks[blnum];
		Blocks[blnum] = blocks[blnum];
	}
}

#if DEBUG

void
SF_ParserBranched::DumpBrackets() {
	int brnum;
	Text s = Copy("Brackets: "), tpos, ttype;
	for (brnum=1;brnum<=numBrackets;brnum++){
		tpos = Blanks(3);
		ttype = Blanks(3);
		tpos.Putint(BracketPos[brnum]);
		ttype.Putint(BracketType[brnum]);
		s = s+" {"+tpos+ttype+"} ";
	}
	Message(debug,"Comp="+Composition+"\n" + s);
}

void
SF_ParserBranched::DumpBlocks() {
	int blnum;
	Text s = Copy("Blocks: "), trank, block;
	for (blnum=1;blnum<=numBlocks;blnum++){
		block = Copy(Blocks[blnum]);
		trank = Blanks(3);
		trank.Putint(BlockRanks[blnum]);
		s = s+" {"+block+trank+"} ";
	}
	Message(debug,"Comp="+Composition+"\n" + s);
}

void
SF_ParserBranched::DumpLinkNodeTree(){
	int blnum,i;
	for (blnum=1;blnum<=numBlocks;blnum++) {
		printf("%d: [%d,%d] l=%d\n",blnum,Struct->LinkNodeTree->LinkList[blnum].node(1),Struct->LinkNodeTree->LinkList[blnum].node(2), Struct->LinkNodeTree->LinkList[blnum].length());
		fflush(stdout);
	}
	for (i=1;i<=numBlocks+1;i++) {
		printf("\n Node %d has %d links: ",i,Struct->LinkNodeTree->NodeList[i].n_links());
		for (blnum=1;blnum<=Struct->LinkNodeTree->NodeList[i].n_links();blnum++) printf("%d ",Struct->LinkNodeTree->NodeList[i].link(blnum));
		fflush(stdout);
	}
}

#endif
