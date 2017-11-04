#include "SF_MolStructure.h"

SF_MolStructure::SF_MolStructure(Text composition_,
								 Text molName_,
								 SF_SegmentList* SegQ_,
								 Lattice* Lat_,
								 Input* MyInput) {

	molName = molName_;
	composition = composition_;
	SegQ = SegQ_;
	Lat = Lat_;

	if (MyInput->GetNumNames("alias") > 0) {
		numAlias = MyInput->GetNumNames("alias");
		AliasNames = MyInput->GetNames("alias");
		AliasComps.Dim(1,numAlias);
		for (int i=1; i<=numAlias; i++) {
			AliasComps[i]=MyInput->GetText("alias", AliasNames[i], "value", molName);
		}
	}

	Text t = Copy(composition);

	char test;
	side_type = 0;
	lin=1,den=2,a_den=3,bra=4;

	int count = 0, bracount = 0, commacount = 0, periodcount=0, semicoloncount=0;
	//! remove ' ' and '-' from the composition; check bracket pairs;
	//! check intersecting [ and ( brackets; count square brackets
	bool substituted=true;
	//bool found;
	t.Setpos(1);
	test = t.Getchar();
	if (test == '@') {
		t.Setpos(1);
		Text s = Copy(t.Scanto('('));
		if (*s != *Copy("dend")) t.Setpos(2);
	}
	int pb;
	int nloops=0;
	while (substituted && nloops < 25) {
		substituted=false; found=false; nloops++;
		pb=0;t.Setpos(2);
		while(t.More()&& nloops < 25) {
			test=t.Getchar();
			if (test=='@') {
				pb = t.Pos()-2;
				Text s = t.Scanto('(').Strip();
				int spos=t.Pos();
				t.Setpos(pb+2);
				Text s1= t.Scanto('@').Strip();
				int s1pos=t.Pos();
				t.Setpos(pb+2);
				Text s2=t.Scanto(')').Strip();
				int s2pos=t.Pos();
				t.Setpos(pb+2);
				Text s3=t.Scanto(',').Strip();
				int s3pos=t.Pos();
				t.Setpos(pb+2);
				Text s4=t.Scanto('.').Strip();
				int s4pos=t.Pos();
				t.Setpos(pb+2);
				Text s5=t.Scanto(';').Strip();
				int s5pos=t.Pos();
				if (*s == *Copy("dend") || *s1 == *Copy("dend") || *s2 == *Copy("dend") || *s3 == *Copy("dend") || *s4 == *Copy("dend") || *s5 == *Copy("dend"))
				{
					Message(fatal,MyInput,"@dend is reserved keyword. Do not use 'dend' as alias.name. In the case you try to get a dendrimer as side group in comb, put only the dendrimer arguments in the 'third' field");}
				else {
					int I=0;
					for (int i=1; i<=numAlias; i++) {
						if (*s == *AliasNames[i]) I=i; else
						if (*s1 == *AliasNames[i]) {spos=s1pos; I=i;} else
						if (*s2 == *AliasNames[i]) {spos=s2pos; I=i;} else
						if (*s3 == *AliasNames[i]) {spos=s3pos; I=i;} else
						if (*s4 == *AliasNames[i]) {spos=s4pos; I=i;} else
						if (*s5 == *AliasNames[i]) {spos=s5pos; I=i;};
					}
					if (I==0) Message(fatal,MyInput,"Alias substitution failed after: " + t.Sub(1,pb));
					int length= AliasNames[I].Length();
					substituted = true; found=true;
					Text tnew = Blanks(t.Length()+AliasComps[I].Length());
					if (spos>t.Length()+1)
						tnew=t.Sub(1,pb)+AliasComps[I];
						else tnew=t.Sub(1,pb)+AliasComps[I]+t.Sub(pb+AliasNames[I].Length()+2,t.Length()-(pb+length+1));
					t=Copy(tnew); t.Strip();
					found=false;
					nloops++;
					t.Setpos(2);
					if (t.Length() >1e6) Message(fatal,MyInput,"Alias substitutions grow out of hand. Substitutions aborted. Possible you have an infinite loop!");
				}

			}
		}
	}
	if (nloops>24) Message(fatal,MyInput,"Nesting limit (25) reached. Substitutions aborted. Check your Aliases: Possible you have defined an infinite loop!");
	composition =Copy(t);
	Text t2 = Blanks(t.Length());
	t.Setpos(1);
	while (t.More()) {
		test=t.Getchar();
		if ((!isspace(test)) && (test!='-')) t2.Putchar(test);
		if (test == '(') count++;
		if (test == ')') count--;
		if ((test == '[') || (test == ']')) bracount++;
		if (test == ',') commacount++;
		if (test == '.') periodcount++;
		if (test == ';') semicoloncount++;
	}
	t = Copy(t2.Frontstrip().Strip());
	if (count < 0) Message(fatal,MyInput,"Error in composition '" + t + "'\nMore brackets closed than open");
	if (count > 0) Message(fatal,MyInput,"Error in composition '" + t + "'\nMore brackets open than closed");
	t.Setpos(1);
	test = t.Getchar();
	if (test == '@') {
		//mmm dendrimer only
		Text s = Copy(t.Scanto('('));
		if (*s != *Copy("dend")) {
			Message(fatal,"syntax error in composition for mol : " + molName);
		}
		s = Copy(t.Scanto(','));
		s = Copy(s.Frontstrip().Strip());
		if (MyInput->CheckInt(s)) {
			MolType = dendrimer;
			numGenerations = s.Getint();
			if (bracount > 0) {
				Message(fatal,"Error in composition for mol : '" + molName +
				"'. Square brackets are not allowed in symmetric dendrimer composition.");
			}
			if (periodcount>0) {
				Message(fatal,"Error in composition for mol : '" + molName +
				"'. No 'periods' allowed.");
			}
			if (semicoloncount>0) {
				Message(fatal,"Error in composition for mol : '" + molName +
				"'. No 'semicolon' allowed.");
			}
			if (commacount%2 == 1) {
				Message(fatal,"Error in composition for mol : '" + molName +
				"'. Number of ',' should be even.");
			}
		} else {
			t.Setpos(1);
			t.Scanto('(');
			s=Copy(t.Scanto(';')); //signature for asymmetric dendrimers.
			s = Copy(s.Frontstrip().Strip());
			if (MyInput->CheckInt(s)) {
				MolType = asymmetric_dendrimer;
				numGenerations = s.Getint();
				if (bracount > 0) {
					Message(fatal,"Error in composition for mol : '" + molName +
					"'. Square brackets are not allowed in asymmetic dendrimer composition.");
				}
				if (periodcount>0) {
					Message(fatal,"Error in composition for mol : '" + molName +
					"'. No 'periods' allowed.");
				}
				if (semicoloncount>1) {
					Message(fatal,"Error in composition for mol : '" + molName +
					"'. Only one 'semicolon' allowed.");
				}
				if (commacount%2 == 0) {
					Message(fatal,"Error in composition for mol : '" + molName +
					"'. Number of ',' should be odd.");
				}
			} else {

				t.Setpos(1);
				t.Scanto('(');
				s=Copy(t.Scanto('.')); //signature for comb.
				s = Copy(s.Frontstrip().Strip());
				MolType = comb;
				if (MyInput->CheckInt(s)) numArms = s.Getint();
				else {
					Message(literal,"Help @dend: \n"
					"There are three options for 'dendrimers'. \n"
					"Symmetric and Asymmetric and Combs. Use the following: \n "
					"Sym dendrimer: @dend(G,T0, f0, T1, f1, ....TG, fG) \n "
					"Aysm dendrimer: @dend(G; S11, S12 ; S21 , S22; S31, S32; ...., SG1, SG2) \n"
					"Comb: @dend(NA. Main1, Spacer, Side, Main2) equivalent to 'Main1-(Spacer-[Side])NA-Main2' \n"
					"with G = generations, NA = repeats of (Spacer,Side), T, S*, M* = valid (linear) composition, f= number. \n"
					"     Sx1 and Sx2 are valid compositions \n"
					"Note: first segments Sx1(1)=Sx2(1) ! \n"
					"TAKE NOTE of the use of ',', '.', and ';' in these definitions.");
					Message(fatal, "try again");
				}
				if (periodcount>1) {
					Message(fatal,"Error in composition for mol : '" + molName +
					"'. Only one 'period' allowed, after the NA, i.e. the number of arms, argument");
				}
				if (semicoloncount>1) {
					Message(fatal,"Error in composition for mol : '" + molName +
					"'. Only one 'semicolon' allowed (in third argument; in case of asymmetric dendrimers in side).");
				}
				if (semicoloncount==0) {
					if (commacount%2 ==0) {
						Message(fatal,"Error in composition for mol : '" + molName +
						"'. Number of ',' should be odd.");
					}
				} else {
					if (commacount%2 ==1) {
						Message(fatal,"Error in composition for mol : '" + molName +
						"'. Number of ',' should be even.");
					}
				}
				SegBlock_combQ.Dim(1,4);
				if (commacount ==3) {side_type = lin; } else {
					t.Scanto(',');
					t.Scanto(',');
					if (semicoloncount ==1) {
						side_type=a_den;
						s=Copy(t.Scanto(';')); //signature for asymmetric dendrimers.
						s = Copy(s.Frontstrip().Strip());
						if (MyInput->CheckInt(s)) {
							numGenerations = s.Getint();
							SegBlockQ.Dim(1,1);
							SegBlock_1Q.Dim(1,numGenerations);
							SegBlock_2Q.Dim(1,numGenerations);
							repeatQ.Dim(1,numGenerations);
						} else 	Message(fatal,"Error in composition for mol : '" + molName +
							"'. Number of generations not found for asymmetric dendrimer side chains in comb.");
					} else {
						side_type=den;
						s=Copy(t.Scanto(',')); //signature for dendrimers.
						s = Copy(s.Frontstrip().Strip());

						if (MyInput->CheckInt(s)) {
							numGenerations=s.Getint();
							SegBlockQ.Dim(1,numGenerations);
							repeatQ.Dim(1,numGenerations); //mmm repeatQ (integers) gives the functionality of the branching points for each generation
						} else 	Message(fatal,"Error in composition for mol : '" + molName +
							"'. Number of generations not found for dendrimer side chains in comb.");
					}
					t.Setpos(1); t.Scanto('('); t.Scanto('.');
				}
			}
		}


		switch (MolType) {
			case comb:
				s = Copy(t.Scanto(','));
				s = Copy(s.Frontstrip().Strip());

				s.Setpos(1);
				while (s.More()) {
					test=s.Getchar();
					if ((test == '[') || (test == ']'))
					Message(fatal,"Error in composition for mol : '" + molName +
					"'. Square brackets are not allowed in argument 1:" + s);
				}
				s.Setpos(1);
				SegBlock_combQ[1] = new SF_SegmentBlock(s,molName,SegQ,MyInput,"","");

				s = Copy(t.Scanto(','));
				s = Copy(s.Frontstrip().Strip());

				s.Setpos(1);
				while (s.More()) {
					test=s.Getchar();
					if ((test == '[') || (test == ']'))
					Message(fatal,"Error in composition for mol : '" + molName +
					"'. Square brackets are not allowed in argument 2:" + s);
				}
				s.Setpos(1);
				SegBlock_combQ[2] = new SF_SegmentBlock(s,molName,SegQ,MyInput,"","");

				if (commacount ==3) {
					s = Copy(t.Scanto(','));
					s = Copy(s.Frontstrip().Strip());
					if (bracount>0) {
						side_type=bra;
						SF_ParserBranched(s, this, MyInput);
						SegBlock_combQ[3]=SegBlock_combQ[1];
					} else
					SegBlock_combQ[3] = new SF_SegmentBlock(s,molName,SegQ,MyInput,"","");
				} else {
					SegBlock_combQ[3]=SegBlock_combQ[1];
					if (semicoloncount ==1) {
						s=Copy(t.Scanto(';')); //signature for asymmetric dendrimers.
						s = Copy(s.Frontstrip().Strip());
						for (int i=1; i<=numGenerations; i++) {
							s = Copy(t.Scanto(','));
							s = Copy(s.Frontstrip().Strip());
							SegBlock_1Q[i] = new SF_SegmentBlock(s,molName,SegQ,MyInput,"","");
							s = Copy(t.Scanto(','));
							s = Copy(s.Frontstrip().Strip());
							s.Setpos(s.Length());
							test = s.Getchar();
							if (test==')') s = Copy(s.Sub(1,s.Length()-1));
							SegBlock_2Q[i] = new SF_SegmentBlock(s,molName,SegQ,MyInput,"","");
							if (i==1) {repeatQ[i] = 1;} else {repeatQ[i]=2;}
							if (!(SegBlock_1Q[i]->GetSegment(1)==SegBlock_2Q[i]->GetSegment(1)))
							Message(fatal,"In Asymmetric dendrimers: first segments of long and short blocks need to be the same. This appears not to be the case in side chain asymmetric dendrimer arguments");
						}
					} else {
						s=Copy(t.Scanto(',')); //signature for dendrimers.

						s = Copy(s.Frontstrip().Strip());
						for (int i=1; i<=numGenerations; i++) {
							s = Copy(t.Scanto(','));
							s = Copy(s.Frontstrip().Strip());
							SegBlockQ[i] = new SF_SegmentBlock(s,molName,SegQ,MyInput,"","");
							s = Copy(t.Scanto(','));
							s = Copy(s.Frontstrip().Strip());
							repeatQ[i] = s.Getint();
						}
					}
				}
				s = Copy(t.Scanto(','));
				s = Copy(s.Frontstrip().Strip());
				s.Setpos(s.Length());
				test = s.Getchar();
				if (test==')') s = Copy(s.Sub(1,s.Length()-1));

				s.Setpos(1);
				while (s.More()) {
					test=s.Getchar();
					if ((test == '[') || (test == ']'))
					Message(fatal,"Error in composition for mol : '" + molName +
					"'. Square brackets are not allowed in argument 4:" + s);
				}
				s.Setpos(1);
				SegBlock_combQ[4] = new SF_SegmentBlock(s,molName,SegQ,MyInput,"","");

			break;
			case asymmetric_dendrimer:
				SegBlockQ.Dim(1,1);
				SegBlock_1Q.Dim(1,numGenerations);
				SegBlock_2Q.Dim(1,numGenerations);
				repeatQ.Dim(1,numGenerations);
				for (int i=1; i<=numGenerations; i++) {
					s = Copy(t.Scanto(','));
					s = Copy(s.Frontstrip().Strip());
					SegBlock_1Q[i] = new SF_SegmentBlock(s,molName,SegQ,MyInput,"","");
					s = Copy(t.Scanto(','));
					s = Copy(s.Frontstrip().Strip());
					s.Setpos(s.Length());
					test = s.Getchar();
					if (test==')') s = Copy(s.Sub(1,s.Length()-1));
					SegBlock_2Q[i] = new SF_SegmentBlock(s,molName,SegQ,MyInput,"","");
					if (i==1) {repeatQ[i] = 1;} else {repeatQ[i]=2;}
					if (!(SegBlock_1Q[i]->GetSegment(1)==SegBlock_2Q[i]->GetSegment(1)))
					Message(fatal,"In Asymmetric dendrimers: first segments of long and short blocks need to be the same. This appears not to be the case...");
				}
			break;
			case dendrimer:
				SegBlockQ.Dim(1,numGenerations);
				repeatQ.Dim(1,numGenerations); //mmm repeatQ (integers) gives the functionality of the branching points for each generation
				for (int i=1; i<=numGenerations; i++) {
					s = Copy(t.Scanto(','));
					s = Copy(s.Frontstrip().Strip());
					SegBlockQ[i] = new SF_SegmentBlock(s,molName,SegQ,MyInput,"","");
					s = Copy(t.Scanto(','));
					s = Copy(s.Frontstrip().Strip());
					repeatQ[i] = s.Getint();
				}
			break;

			default:
				Message(fatal,"Error in analyse dendrimer.");
			break;
		}
		SetMolLength();
		SetDiffSegmentQ();
		SetPhiRef();
	}
	else {
		if (bracount > 0) {
			SF_ParserBranched(t, this, MyInput);  //Initialize LinkNodeTree, SegBlockQ, repeatQ, numGenerations
			MolType = branched;
			SetMolLength();	//watch out - dendrimer-specific
			SetDiffSegmentQ();
			SetPhiRef();
		}
		else {
		//mmm if not a dendrimer or a branched molecule
			t.Setpos(1);
			SF_SegmentBlock* SegBlock;
			SegBlock = new SF_SegmentBlock(t,molName,SegQ,MyInput,"","");
			SegBlockQ.Dim(1,1);
			repeatQ.Dim(1,1);
			repeatQ[1] = 1;
			SegBlockQ[1] = SegBlock;
			numDiffSegments = SegBlock->GetNumDiffSegments();
			molLength = SegBlock->GetAvLength();
			DiffSegmentQ.Dim(1,numDiffSegments);
			numGenerations = 1;
			int i;
			SF_MolSegment* MolSeg;
			double number;
//			equivalent to SetPhiRef(), but does not use SegBlockQ
			for (i=1; i<=numDiffSegments; i++) {
				MolSeg = DiffSegmentQ[i] = SegBlock->GetDiffSegment(i);
				number = SegBlock->GetAvNumSegments(MolSeg);
				MolSeg->SetPhiRef(number/molLength);
			}
			if (molLength == 1) {
				MolType = monomer;
			}
			else {
				if (numDiffSegments == 1) {
					MolType = homopolymer;
				}
				else {
					MolType = copolymer;
				}
			}
		}
	}
	chemIntRef = 0;
	int i,j,numStates;
	SF_MolSegment* MolSeg;
	SF_MolState* MolState;
	for (i=1; i<=numDiffSegments; i++) {
		MolSeg = DiffSegmentQ[i];
		numStates = MolSeg->GetNumStates();
		for (j=1; j<=numStates; j++) {
			MolState = MolSeg->GetState(j);
			chemIntRef += MolState->GetPhiRef()*SegQ->ChemIntRef(MolState);
		}
	}
	chemIntRef /= 2;
	someSegmentsPinned = false;
	for (i=1; i<=numDiffSegments; i++) {
		MolSeg = GetDiffSegment(i);
		if (MolSeg->GetFreedom() == pinned) {
			someSegmentsPinned = true;
		}
	}
	someSegmentsMC = false;
	someSegmentsLD = false;
	if (someSegmentsPinned) {
		for (i=1; i<=numDiffSegments; i++) {
			MolSeg = GetDiffSegment(i);
			if (MolSeg->GetMC()) {
				someSegmentsMC = true;
			}
			if (MolSeg->GetLD()) {
				someSegmentsLD = true;
			}
		}
	}
	someSegmentsGrafted = false;
	for (i=1; i<=numDiffSegments; i++) {
		MolSeg = GetDiffSegment(i);
		if (MolSeg->GetFreedom() == grafted && !someSegmentsGrafted) {
			someSegmentsGrafted = true;
			GraftedSegment = MolSeg;
		} else if (MolSeg->GetFreedom() == grafted && someSegmentsGrafted) {
			Message(warning,"Thermodynamics is only right "
			"for one segment per molecule which has freedom : grafted. "
			"You probably want to set the freedom to 'pinned'.");
		}

	}

}
SF_MolStructure::~SF_MolStructure() {
	int i;
	for (i=1; i<=numDiffSegments; i++) {
		delete DiffSegmentQ[i];
	}
	for (i=1; i<=SegBlockQ.Upperbound(); i++) {
		delete SegBlockQ[i];
	}
//	if (MolType = branched) delete LinkNodeTree;

}
double
SF_MolStructure::GetAvLength() const {
	return molLength;
}
int
SF_MolStructure::GetMaxLength() const {
	return int(molLength);
}
SF_MolSegment*
SF_MolStructure::GetSegment(int number) const {
	return SegBlockQ[1]->GetSegment(number);
}

SF_MolSegment*
SF_MolStructure::GetSegment(int linknum, int r) const {
	#if DEBUG
	//!range check
		if ((linknum<1)||(linknum>numGenerations)) {
			printf("Programming error in SF_MolStructure::GetSegment(linknum=%d,segnum=%d)\n",linknum,r);
			fflush(stdout);
			Message(fatal, "linknum out of range");
		}
		else if((r > LinkNodeTree->LinkList[linknum].length()) || (r<1) )
		{
			printf("Programming error in SF_MolStructure::GetSegment(linknum=%d,segnum=%d)\n",linknum,r);
			fflush(stdout);
			Message(fatal, "segmentnum out of range");
		};

	#endif

	if (linknum == 1) return (SegBlockQ[linknum])->GetSegment(r);
	if (r == 1){
		int parentnode, parentlink;
		parentnode = LinkNodeTree->LinkList[linknum].node(1);
		parentlink = LinkNodeTree->NodeList[parentnode].link(1);
		return (SegBlockQ[parentlink])->GetSegment(SegBlockQ[parentlink]->GetMaxLength());
	}
	return (SegBlockQ[linknum])->GetSegment(r-1);
}
Text
SF_MolStructure::GetComposition() const {
	return composition;
}
int
SF_MolStructure::GetNumDiffSegments() const {
	return numDiffSegments;
}
SF_MolSegment*
SF_MolStructure::GetDiffSegment(int number) const {
	return DiffSegmentQ[number];
}
int
SF_MolStructure::GetNumSegmentBlocks() const {
	return SegBlockQ.Upperbound();
}
SF_SegmentBlock*
SF_MolStructure::GetSegmentBlock(int i) const {
	return SegBlockQ[i];
}
SF_SegmentBlock*
SF_MolStructure::GetSegmentBlockComb(int i) const {
	return SegBlock_combQ[i];
}

SF_SegmentBlock*
SF_MolStructure::GetSegmentBlock(int i, int j) const {

	int length1 = SegBlock_1Q[i]->GetAvLength();
	int length2 = SegBlock_2Q[i]->GetAvLength();

	if (j==1) {
		if (length1<length2) return SegBlock_1Q[i]; else return SegBlock_2Q[i];
	}
	if (j==2) {
		if (length1<length2) return SegBlock_2Q[i]; else return SegBlock_1Q[i];
	}
	return NULL;
}

double
SF_MolStructure::GetAvNumSegments(const SF_MolSegment* MolSeg) const {

	double number = 0;

	switch (MolType) {
		case comb:
			if (side_type==lin) {
				number=SegBlock_combQ[3]->GetAvNumSegments(MolSeg);
			} else if (side_type == a_den) {
				for (int j=numGenerations; j>=1; j--) {
					if (MolSeg == SegBlock_1Q[j]->GetSegment(1)) {
						number = (number + SegBlock_1Q[j]->GetAvNumSegments(MolSeg)+ SegBlock_2Q[j]->GetAvNumSegments(MolSeg) - 1)*repeatQ[j];
					} else number = (number + SegBlock_1Q[j]->GetAvNumSegments(MolSeg)+SegBlock_2Q[j]->GetAvNumSegments(MolSeg))*repeatQ[j];
				}
			} else {
				for (int j=numGenerations; j>=1; j--) {
					if (MolSeg == SegBlockQ[j]->GetSegment(1)) {
						number = (number + SegBlockQ[j]->GetAvNumSegments(MolSeg) - 1)*repeatQ[j]+1;
					} else {
						number = (number + SegBlockQ[j]->GetAvNumSegments(MolSeg))*repeatQ[j];
					}
				}
			}
			number = SegBlock_combQ[1]->GetAvNumSegments(MolSeg) +
		             SegBlock_combQ[4]->GetAvNumSegments(MolSeg) +
		             (SegBlock_combQ[2]->GetAvNumSegments(MolSeg)+number)*numArms;

		break;
		case asymmetric_dendrimer:
			for (int j=numGenerations; j>=1; j--) {
				if (MolSeg == SegBlock_1Q[j]->GetSegment(1)) {
					number = (number + SegBlock_1Q[j]->GetAvNumSegments(MolSeg)+ SegBlock_2Q[j]->GetAvNumSegments(MolSeg) - 1)*repeatQ[j];
				} else number = (number + SegBlock_1Q[j]->GetAvNumSegments(MolSeg)+SegBlock_2Q[j]->GetAvNumSegments(MolSeg))*repeatQ[j];
			}
		break;
		default:
			for (int j=numGenerations; j>=1; j--) {
				if (MolSeg == SegBlockQ[j]->GetSegment(1)) {
					number = (number + SegBlockQ[j]->GetAvNumSegments(MolSeg) - 1)*repeatQ[j]+1;
				} else {
					number = (number + SegBlockQ[j]->GetAvNumSegments(MolSeg))*repeatQ[j];
				}
			}
		break;
	}
	return number;
}

MoleculeType
SF_MolStructure::GetMoleculeType() const {
	return MolType;
}
void
SF_MolStructure::SetPhiBulk(double phiBulk) {
	SF_MolSegment* MolSeg;
	for (int i=1; i<=numDiffSegments; i++) {
		MolSeg = GetDiffSegment(i);
		double number = GetAvNumSegments(MolSeg);
		MolSeg->SetPhiBulk((phiBulk*number)/molLength);
	}
}

void
SF_MolStructure::SetPhiBulkBoundaries(Vector phiBulk) {
	SF_MolSegment* MolSeg;
	for (int i=1; i<=numDiffSegments; i++) {
		MolSeg = GetDiffSegment(i);
		double number = GetAvNumSegments(MolSeg);
		int numBounds = 2*Lat->GetNumGradients();
		Vector phiBulkBounds(1,numBounds);
		for (int j=1; j<=numBounds; j++) {
			phiBulkBounds[j] = phiBulk[j]*number/molLength;
		}
		MolSeg->SetPhiBulkBoundaries(phiBulkBounds);
	}
}
double
SF_MolStructure::ChemIntRef() const {
	return chemIntRef;
}
Boolean
SF_MolStructure::Symmetric() const {
	int N1,N2,N4,N;
	SF_MolSegment* MolSeg1;
	SF_MolSegment* MolSeg2;
	switch (MolType) {
		case comb:
			N1=SegBlock_combQ[1]->GetAvLength();
			N2=SegBlock_combQ[2]->GetAvLength();
		    N4=SegBlock_combQ[4]->GetAvLength();
		    if (N4 != N1+N2-1) return false;
		    N=N1+N2+N4;
			for (int s=1; s<N/2; s++) {
				if (s<=N1) MolSeg1=SegBlock_combQ[1]->GetSegment(s); else
				MolSeg1=SegBlock_combQ[2]->GetSegment(s-N1);
				MolSeg2=SegBlock_combQ[4]->GetSegment(N4-s+1);
				if (!(MolSeg1==MolSeg2)) {return false;}
			}
			if (numArms>1) {
				for (int s=1; s<N2; s++) {
					MolSeg1=SegBlock_combQ[2]->GetSegment(s);
					MolSeg2=SegBlock_combQ[2]->GetSegment(N2-s);
					if (!(MolSeg1==MolSeg2)) {return false;}
				}
			}
			return true;
			break;

		default:
			if (numGenerations == 1) {
				return (SegBlockQ[1])->Symmetric();
			} else return false;
			break;
	}
	return false;
}
Boolean
SF_MolStructure::SomeSegmentsPinned() const {
	return someSegmentsPinned;
}

Boolean
SF_MolStructure::SomeSegmentsMC() const {
	return someSegmentsMC;
}

Boolean
SF_MolStructure::SomeSegmentsLD() const {
	return someSegmentsLD;
}

Boolean
SF_MolStructure::SomeSegmentsGrafted() const {
	return someSegmentsGrafted;
}
SF_MolSegment*
SF_MolStructure::GetGraftedSegment() const {
	return GraftedSegment;
}
Boolean
SF_MolStructure::AllInteractionsEqual() const {
	SF_MolSegment* MolSeg;
	SF_MolSegment* MolSeg2;
	SF_State* State;
	SF_State* State2;
	for (int i=1; i<=numDiffSegments; i++) {
		MolSeg = GetDiffSegment(i);
		for (int j=i+1; j<=numDiffSegments; j++) {
			MolSeg2 = GetDiffSegment(j);
			if (MolSeg->GetNumStates() != MolSeg2->GetNumStates()) {
				return false;
			}
			for (int k=1; k<=MolSeg->GetNumStates(); k++) {
				State = SegQ->GetState(MolSeg->GetState(k)->GetName());
				State2 = SegQ->GetState(MolSeg2->GetState(k)->GetName());
				if (!SegQ->ChemIntStatesEqual(State,State2)) {
					return false;
				}
				if (State->GetValence() != State2->GetValence()) {
					return false;
				}
				if (State->GetEpsilon() != State2->GetEpsilon()) {
					return false;
				}
			}
		}
	}
	return true;
}
int
SF_MolStructure::GetNumGenerations() const {
	return numGenerations;
}
int
SF_MolStructure::GetNumRepeat(int i) const {
	return repeatQ[i];
}
int
SF_MolStructure::GetNumArms() const {
	return numArms;
}

int
SF_MolStructure::GetSideType() const {
	//1 = linear sides
	//2 = dendrimer
	//3 = asymmetric_dendrimer
	//4 = branched
	return side_type;
}

void
SF_MolStructure::SetDiffSegmentQ(){
	numDiffSegments = 0;
	switch (MolType) {
		case comb:
			for (int i=1; i<5; i++){
				numDiffSegments += SegBlock_combQ[i]->GetNumDiffSegments();
			}

			if (side_type == a_den) {
				for (int i=1; i<=numGenerations; i++) {
					numDiffSegments += SegBlock_1Q[i]->GetNumDiffSegments()+SegBlock_2Q[i]->GetNumDiffSegments();
				}

			} else if (side_type == den || side_type == bra){
				for (int i=1; i<=numGenerations; i++) {
					numDiffSegments += SegBlockQ[i]->GetNumDiffSegments();
				}

			}

			DiffSegmentQ.Dim(1,numDiffSegments);
			numDiffSegments = 0;
			for (int i=1; i<5; i++){
				for (int j=1; j<=SegBlock_combQ[i]->GetNumDiffSegments(); j++) {
					DiffSegmentQ[++numDiffSegments] = SegBlock_combQ[i]->GetDiffSegment(j);
				}
			}
			if (side_type == a_den) {
				for (int i=1; i<=numGenerations; i++) {
					for (int j=1; j<=SegBlock_1Q[i]->GetNumDiffSegments(); j++) {
						DiffSegmentQ[++numDiffSegments] = SegBlock_1Q[i]->GetDiffSegment(j);
					}
					for (int j=1; j<=SegBlock_2Q[i]->GetNumDiffSegments(); j++) {
						DiffSegmentQ[++numDiffSegments] = SegBlock_2Q[i]->GetDiffSegment(j);
					}
				}
			} else if (side_type==den || side_type ==bra) {
				for (int i=1; i<=numGenerations; i++) {
					for (int j=1; j<=SegBlockQ[i]->GetNumDiffSegments(); j++) {
						DiffSegmentQ[++numDiffSegments] = SegBlockQ[i]->GetDiffSegment(j);
					}
				}

			}

		break;
		case asymmetric_dendrimer:
			for (int i=1; i<=numGenerations; i++) {
				numDiffSegments += SegBlock_1Q[i]->GetNumDiffSegments()+SegBlock_2Q[i]->GetNumDiffSegments();
			}
			DiffSegmentQ.Dim(1,numDiffSegments);
			numDiffSegments = 0;
			for (int i=1; i<=numGenerations; i++) {
				for (int j=1; j<=SegBlock_1Q[i]->GetNumDiffSegments(); j++) {
					DiffSegmentQ[++numDiffSegments] = SegBlock_1Q[i]->GetDiffSegment(j);
				}
				for (int j=1; j<=SegBlock_2Q[i]->GetNumDiffSegments(); j++) {
					DiffSegmentQ[++numDiffSegments] = SegBlock_2Q[i]->GetDiffSegment(j);
				}
			}
		break;
		default:
			for (int i=1; i<=numGenerations; i++) {
				numDiffSegments += SegBlockQ[i]->GetNumDiffSegments();
			}
			DiffSegmentQ.Dim(1,numDiffSegments);
			numDiffSegments = 0;
			for (int i=1; i<=numGenerations; i++) {
				for (int j=1; j<=SegBlockQ[i]->GetNumDiffSegments(); j++) {
					DiffSegmentQ[++numDiffSegments] = SegBlockQ[i]->GetDiffSegment(j);
				}
			}
		break;
	}

	SF_MolSegment* MolSeg1;
	SF_MolSegment* MolSeg2;
	// remove duplicate entries from the DiffSegmentQ
	for (int i=1; i<=numDiffSegments; i++) {
		MolSeg1 = DiffSegmentQ[i];
		for (int j=i+1; j<=numDiffSegments; j++) {
			MolSeg2 = DiffSegmentQ[j];
			if (MolSeg1 == MolSeg2) {
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
	for (int i=1; i<=numDiffSegments; i++) {
		MolSeg1 = DiffSegmentQ[i];
		double number = GetAvNumSegments(MolSeg1);
		MolSeg1->SetPhiRef(number/molLength);
	}
}


void
SF_MolStructure::SetMolLength(){
	int i;
	molLength = 0;
	switch (MolType) {
		case comb:
			if (side_type==lin) molLength = SegBlock_combQ[3]->GetAvLength(); else
			if (side_type == a_den) {
				for (i=numGenerations; i>=1; i--) {
					int REPEAT=2;
					molLength = molLength*REPEAT + SegBlock_1Q[i]->GetAvLength() + SegBlock_2Q[i]->GetAvLength()-1;
				}
			} else {
				for (i=numGenerations; i>=1; i--)
				molLength = (molLength + SegBlockQ[i]->GetAvLength() - 1)*repeatQ[i] + 1;
			}

			molLength = (molLength+SegBlock_combQ[2]->GetAvLength())*numArms +
						SegBlock_combQ[1]->GetAvLength() +
					    SegBlock_combQ[4]->GetAvLength();
		break;
		case asymmetric_dendrimer:
			for (i=numGenerations; i>=1; i--) {
			int REPEAT=2;
			molLength = molLength*REPEAT + SegBlock_1Q[i]->GetAvLength() + SegBlock_2Q[i]->GetAvLength()-1;
		}
		break;
		default:
			for (i=numGenerations; i>=1; i--) {
			molLength = (molLength + SegBlockQ[i]->GetAvLength() - 1)*repeatQ[i] + 1;
		}
		break;
	}
}

void
SF_MolStructure::SetPhiRef(){
	int i;
	SF_MolSegment* MolSeg1;
	for (i=1; i<=numDiffSegments; i++) {
		MolSeg1 = DiffSegmentQ[i];
		double number = GetAvNumSegments(MolSeg1);
		MolSeg1->SetPhiRef(number/molLength);
	}
}
SF_LinkNodeTree*
SF_MolStructure::GetLinkNodeTree() const {
	return LinkNodeTree;
}

