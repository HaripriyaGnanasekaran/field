#include "SF_SegmentBlock.h"

SF_SegmentBlock::SF_SegmentBlock(Text t, 
								 Text molName,
								 SF_SegmentList* SegQ_,
								 Input* MyInput_,
								 Text before,
								 Text after) {
	blockName = Copy(t);
	nameBlockBefore = Copy(before);
	nameBlockAfter = Copy(after);
	SegQ = SegQ_;
	MyInput = MyInput_;
	char test;
	repeat = "";
	numBlocks = CountNumBlocks(t,molName);
	if (numBlocks == 0) {
		t.Setpos(2);
		test = t.Getchar();
		Text segName = Blanks(t.Length());
		while (test != ')') {
			segName.Putchar(test);
			test = t.Getchar();
		}
		segName = Copy(segName.Strip());
		repeated = Copy(segName);
		repeat = ReadRepeat(t,molName);
		numDiffSegments = 1;
		SF_MolSegment* MolSeg = SegQ->NewMolSegment(segName,molName);
		DiffSegmentQ.Dim(1,1);
		DiffSegmentQ[1] = MolSeg;
	} else if (numBlocks == 1) {
		BlockQ.Dim(1,1);
		AssignOneBlock(t,molName);
		AssignDiffSegmentQ();
	} else {
		BlockQ.Dim(1,numBlocks);
		AssignMultiBlock(t,molName);
		AssignDiffSegmentQ();
	}
	if (repeat.More()) {
		test = repeat.Getchar();
		if (isdigit(test)) {
			polydisperse = false;
			--repeat;
			avRepeat = minRepeat = maxRepeat = repeat.Getint();
		}
		/*
		else {
			polydisperse = true;
			Polydisp = NewPolydispersity(repeat,molName);
			avRepeat = Polydisp->GetAvRepeat();
			maxRepeat = Polydisp->GetMaxRepeat();
			minRepeat = Polydisp->GetMinRepeat();
		}
		*/
	} else {
		polydisperse = false;
		avRepeat = minRepeat = maxRepeat = 1;
	}
	if (numBlocks == 0) {
		avLength = avRepeat;
		minLength = minRepeat;
		maxLength = maxRepeat;
	} else {
		avLength = minLength =  maxLength = 0;
		for (int i=1; i<=numBlocks; i++) {
			avLength += (BlockQ[i])->GetAvLength()*avRepeat;
			minLength += (BlockQ[i])->GetMinLength()*minRepeat;
			maxLength += (BlockQ[i])->GetMaxLength()*maxRepeat;
		}
	}
}
SF_SegmentBlock::~SF_SegmentBlock() {
	SF_SegmentBlock *bl;
	for (int i=1; i<= numBlocks; i++) {
		bl = BlockQ[i];
		delete bl;
	}
}
Text
SF_SegmentBlock::GetBlockName() const {
	return blockName;
}
Text
SF_SegmentBlock::GetBlockRepeat() const {
	return repeat;
}
double
SF_SegmentBlock::GetAvLength() const {
	return avLength;
}
int
SF_SegmentBlock::GetMinLength() const {
	return minLength;
}
int
SF_SegmentBlock::GetMaxLength() const {
	return maxLength;
}
double
SF_SegmentBlock::GetAvRepeat() const {
	return avRepeat;
}
int
SF_SegmentBlock::GetMinRepeat() const {
	return minRepeat;
}
int
SF_SegmentBlock::GetMaxRepeat() const {
	return maxRepeat;
}
SF_MolSegment*
SF_SegmentBlock::GetSegment(int s) const {
	SF_SegmentBlock* bl;
	int i;
	if (numBlocks == 1) {
		bl = BlockQ[1];
		int length = bl->GetMaxLength();
		for (i=1; i<=maxRepeat; i++) {
			s-=length;
			if (s<=0) {
				return bl->GetSegment(s+length);
			}
		}
	}
	for (i=1; i<= numBlocks; i++) {
		bl = BlockQ[i];
		s -= bl->GetMaxLength();
		if (s <= 0) {
			s += bl->GetMaxLength();
			return bl->GetSegment(s);
		}
	}
	return DiffSegmentQ[1];
}
int
SF_SegmentBlock::GetNumBlocks() const {
	return numBlocks;
}
SF_SegmentBlock*
SF_SegmentBlock::GetBlock(int number) const {
	return BlockQ[number];
}
int
SF_SegmentBlock::GetNumDiffSegments() const {
	return numDiffSegments;
}
SF_MolSegment*
SF_SegmentBlock::GetDiffSegment(int number) const {
	return DiffSegmentQ[number];
}
int
SF_SegmentBlock::GetAvNumSegments(const SF_MolSegment* MolSeg) const {
	int number = 0;//verander voor polydisp
	if (numBlocks == 0 && MolSeg == DiffSegmentQ[1]) {
		return maxLength;
	}
	for (int i=1; i<=numBlocks; i++) {
		number += maxRepeat*(BlockQ[i])->GetAvNumSegments(MolSeg);
	}
	return number;
}
Boolean 
SF_SegmentBlock::Symmetric() const {
	for (int i=1; i<=maxLength; i++) {
		if (GetSegment(i) != GetSegment(maxLength - i +1)) {
			return false;
		}
	}
	return true;
}
Boolean 
SF_SegmentBlock::PolydisperseChain() const {
	if (numBlocks == 0) {
		return polydisperse;
	}
	Boolean value = false;
	for (int i=1; i<=numBlocks; i++) {
		if ((BlockQ[i])->PolydisperseChain()) {
			value = true;
		}
	}
	return value;
}
Boolean 
SF_SegmentBlock::Polydisperse() const {
	return polydisperse;
}
/*
double 
SF_SegmentBlock::GetFraction(int s) const {
	return Polydisp->GetFraction(s);
}
*/
Array<Text>
SF_SegmentBlock::GetFractionNames() const {
	Array<Text> Names;
	Array<Text> Temp;
	int numNames = 0;
	if (numBlocks == 0) {
		numNames = maxRepeat - minRepeat + 1;
		Names.Dim(1,numNames);
		for (int i=1; i<=numNames; i++) {
			Names[i] = "(" + repeated + ")";
			Text temp = Blanks(100);
			temp.Putint(minRepeat + i - 1);
			Names[i] = Names[i] + Copy(temp.Frontstrip());
		}
	} else if (numBlocks == 1) {
		SF_SegmentBlock* Block = BlockQ[1];
		Temp = Block->GetFractionNames();
		int numFrBlock = Temp.Upperbound();
		numNames = (maxRepeat - minRepeat + 1)*numFrBlock;
		Names.Dim(1,numNames);
		int count = 0;
		for (int i=1; i<=maxRepeat - minRepeat + 1; i++) {
			for (int j=1; j<=numFrBlock; j++) {
				count++;
				Names[count] = "(" + Temp[j] + ")";
				Text temp = Blanks(100);
				temp.Putint(minRepeat + i - 1);
				Names[count] = Names[count] + Copy(temp.Frontstrip());
			}
		}
	} else {
		int i,j,k;
		Boolean polydisp = false;
		for (i=1; i<=numBlocks; i++) {
			Temp = (BlockQ[i])->GetFractionNames();
			if ((BlockQ[i])->PolydisperseChain()) {
				numNames += Temp.Upperbound();
				polydisp = true;
			}
		}
		if (!polydisp) {
			Names.Dim(1,1);
			Names[1] = "";
			for (i=1; i<=numBlocks; i++) {
				Names[1] = Names[1] + (BlockQ[i])->GetBlockName();
			}
			return Names;
		}
		Names.Dim(1,numNames);
		numNames = 0;
		int count = 0;
		for (i=1; i<=numBlocks; i++) {
			if (!(BlockQ[i])->PolydisperseChain()) {
				continue;
			}
			Temp = (BlockQ[i])->GetFractionNames();
			numNames += count;
			count = Temp.Upperbound();
			for (j=1; j<=count; j++) {
				Names[numNames+j] = "";
			}
			for (int j=1; j<=numBlocks; j++) {
				if (i==j) {
					for (k=1; k<=count; k++) {
						Names[numNames+k] = Names[numNames+k] + Copy(Temp[k]);
					}
				} else {
					for (k=1; k<=count; k++) {
						Names[numNames+k] = Names[numNames+k] + (BlockQ[j])->GetBlockName();
					}
				}	
			}
		}
	}
	return Names;
}
int
SF_SegmentBlock::CountNumBlocks(Text& t, Text molName) const {
	int value = 0;
	t.Setpos(1);
	Boolean oneSegment = true;
	char test;
	while(t.More()) {
		test = t.Getchar();
		if (test == '(') {
			value++;
			int count = 1;
			while (count > 0) {
				if (!t.More()) {
					Message(fatal, MyInput, "Syntax error reading block '" 
						+ t + "' in composition for molecule '" 
						+ molName + "'");
				}
				test = t.Getchar();
				if (test == '(') {
					oneSegment = false;
					count++;
				}
				if (test == ')') {
					count--;
				}
			}
			Text dummy = Copy(ReadRepeat(t,molName));
		} else {
			if (SegQ->SegmentDefined(t)) {
				t = "(" + t + ")1";
				return 0;
			}
			if (t.Pos() < 3) {
				Message(fatal, MyInput, "Syntax error reading block '" 
					+ t + "' in composition for molecule '" 
					+ molName + "'\nProbably you try to use the old style for"
					" defining compositions.\nThe old style a10b10 should now "
					"be written as (a)10(b)10");
			} else {
				Message(fatal, MyInput, "error reading composition '" 
					+ t + "' for molecule '" 
					+ molName + "'");
			}
		}
	}
	if (value > 1) {
		return value;
	} else {
		if (oneSegment) return 0;
		else return 1;
	}
}
void
SF_SegmentBlock::AssignOneBlock(Text t, Text molName) {
	Text block(Blanks(t.Length()));
	t.Setpos(1);
	char test = t.Getchar();
	while(t.More()) {
		if (test == '(') {
			int count = 1;
			while (count > 0) {
				test = t.Getchar();
				if (test == '(') count++;
				if (test == ')') count--;
				if (count > 0) block.Putchar(test);
			}
			repeated = Copy(block.Strip());
			repeat = ReadRepeat(t,molName);
			BlockQ[1] = new SF_SegmentBlock(repeated,molName,SegQ,MyInput,
				nameBlockBefore+"(",")"+repeat+nameBlockAfter);
		}
	}
	SF_SegmentBlock* Block = BlockQ[1];
	avLength = Block->GetAvLength()*Block->GetAvRepeat();
	minLength = Block->GetMinLength()*Block->GetMinRepeat();
	maxLength = Block->GetMaxLength()*Block->GetMaxRepeat();
}
void
SF_SegmentBlock::AssignMultiBlock(Text t, Text molName) {
	int number = 0;
	t.Setpos(1);
	char test;
	Text before = "";
	Text after;
	while(t.More()) {
		test = t.Getchar();
		if (test == '(') {
			number++;
			int count = 1;
			Text block(Blanks(t.Length()));
			block.Putchar(test);
			while (count > 0) {
				test = t.Getchar();
				if (test == '(') count++;
				if (test == ')') count--;
				block.Putchar(test);
			}
			block = Copy(block.Strip()) + ReadRepeat(t,molName);
			after = Copy(t.Sub(t.Pos(),t.Length()-t.Pos()+1));
			BlockQ[number] = new SF_SegmentBlock(block,molName,SegQ,MyInput,
				nameBlockBefore+before,after+nameBlockAfter);
			before = before + block;
				
		}
	}
	SF_SegmentBlock* Block;
	avLength = 0;
	minLength = 0;
	maxLength = 0;
	for (int i=1; i<=numBlocks; i++) {
		Block = BlockQ[i];
		avLength += Block->GetAvLength()*Block->GetAvRepeat();
		minLength += Block->GetMinLength()*Block->GetMinRepeat();
		maxLength += Block->GetMaxLength()*Block->GetMaxRepeat();
	}		
}
void
SF_SegmentBlock::AssignDiffSegmentQ() {
	if (numBlocks == 0) {
		return;
	}
	SF_SegmentBlock* Block;
	SF_MolSegment* MolSeg;
	SF_MolSegment* MolSeg2;
	int i,j,count=0;
	numDiffSegments = 0;
	for (i=1; i<=numBlocks; i++) {
		Block = BlockQ[i];
		int numDiffSeg = Block->GetNumDiffSegments();
		for (j=1; j<=numDiffSeg; j++) {
			numDiffSegments++;
		}
	}
	DiffSegmentQ.Dim(1,numDiffSegments);
	for (i=1; i<=numBlocks; i++) {
		Block = BlockQ[i];
		int numDiffSeg = Block->GetNumDiffSegments();
		for (j=1; j<=numDiffSeg; j++) {
			DiffSegmentQ[++count] = Block->GetDiffSegment(j);
		}
	}
	for (i=1; i<=numDiffSegments; i++) {
		MolSeg = (SF_MolSegment*) DiffSegmentQ[i];
		for (int j=i+1; j<=numDiffSegments; j++) {
			MolSeg2 = DiffSegmentQ[j];
			if (MolSeg == MolSeg2) {
				Array<SF_MolSegment*> DiffSegmentQTemp(1,numDiffSegments-1);
				for (int k=1; k<=numDiffSegments; k++) {
					if (k<j) DiffSegmentQTemp[k] = DiffSegmentQ[k];
					if (k>j) DiffSegmentQTemp[k-1] = DiffSegmentQ[k];
				}
				numDiffSegments--;
				j--;
				DiffSegmentQ = DiffSegmentQTemp;
			}
		}
	}
}
// ? Reads how many times a monomer is repeated: (A)10 is repeated 10 times
Text
SF_SegmentBlock::ReadRepeat(Text& t, Text molName) const {
	Text value = Blanks(t.Length());
	char test;
	if (t.More()) {
		test = t.Getchar();
	} else {
		return Copy(value.Strip());
	}
	if (test == '(') {
		--t;
		return Copy(value.Strip());
	} else if (isalpha(test)) {
		value.Putchar(test);
		if (!t.More()) {
			Message(fatal,"error reading composition for molecule '" 
			+ molName + "' (block: " + t + ")");
		}
		test = t.Getchar();
		value.Putchar(test);
		while(test != '(') {
			if (!t.More()) {
				Message(fatal,"error reading composition for molecule '" 
					+ molName + "' (block: " + t + ")");
			}
			test = t.Getchar();
			value.Putchar(test);
		}
		if (!t.More()) {
			Message(fatal,"error reading composition for molecule '" 
				+ molName + "' (block: " + t + ")");
		}
		test = t.Getchar();
		value.Putchar(test);
		while(test != ')') {
			if (!t.More()) {
				Message(fatal,"error reading composition for molecule '" 
					+ molName + "' (block: " + t + ")");
			}
			test = t.Getchar();
			value.Putchar(test);
		}
		return Copy(value.Strip());
	} else if (isdigit(test)) {
		while(isdigit(test)) {
			value.Putchar(test);
			if (t.More()) {
				test = t.Getchar();
			} else {
				break;
			}
		}
		if (!isdigit(test)) {
			--t;
		}
		return Copy(value.Strip());
	}
	Message(fatal,"programming error in SF_SegmentBlock::ReadRepeat");
	return Copy(""); // never get here
}
/*
Polydispersity*
SF_SegmentBlock::NewPolydispersity(Text t, Text molName) {
	t.Setpos(1);
	char test = t.Getchar();
	Text type = Blanks(t.Length());
	Text parameters = Blanks(t.Length());
	while (test != '(') {
		type.Putchar(test);
		if (t.More()) {
			test = t.Getchar();
		} else {
			Message(fatal,"invalid polydispersity '" + t 
				+ "' in composition for molecule '" + molName + "'");
		}
	}
	if (t.More()) {
		test = t.Getchar();
	} else {
		Message(fatal,"invalid polydispersity '" + t 
			+ "' in composition for molecule '" + molName + "'");
	}
	while (test != ')') {
		parameters.Putchar(test);
		if (t.More()) {
			test = t.Getchar();
		} else {
			Message(fatal,"invalid polydispersity '" + t 
				+ "' in composition for molecule '" + molName + "'");
		}
	}
	Sysout().Outtext("param : " + parameters);
	Sysout().Outimage();
	Sysout().Outtext("molName : " + molName);
	Sysout().Outimage();
	if (*Copy(type.Strip()) == *Copy("Flory")) {
		return new PolydispFlory(parameters,molName);
	} else {
		Message(fatal,"invalid polydispersity '" + t 
			+ "' in composition for molecule '" + molName + "'");
		return NULL; //never get here
	}
}
*/
Boolean
SF_SegmentBlock::operator==(const SF_SegmentBlock& b) const {
	if (numBlocks == 0 && b.numBlocks == 0) {	//verander voor polydisp
		if (DiffSegmentQ[1] == b.DiffSegmentQ[1] && avLength == b.avLength) {
			return true;
		}
	} 
	if (numBlocks != b.numBlocks) {
		return false;
	} else {
		int i;
		SF_SegmentBlock *bl1, *bl2;
		for (i=1; i<= numBlocks; i++) {
			bl1 = BlockQ[i];
			bl2 = b.BlockQ[i];
			if (*bl1 != *bl2) {
				return false;
			}
		}
		return true;
	}
}
Boolean
SF_SegmentBlock::operator!=(const SF_SegmentBlock& b) const {
	return !(*this == b);
}
