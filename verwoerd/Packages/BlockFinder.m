(* Mathematica Package *)

(* COPYRIGHT (C) 2010 Wynand Verwoerd 
under the terms of the GNU General Public License V3
This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
*)

BeginPackage["BlockFinder`",{"NetSplitter`"}]
(* Exported symbols added here with SymbolName::usage *)
Blacken::usage="Binarizes and fills up row and column gaps to display blocks"  
BandFinder::usage="Finds bands that are zero except for a block in a binary matrix "
BlockFinder::usage="Locates blocks in a matrix and arranges them on the diagonal"

Begin["`Private`"] (* Begin Private Context *) 

Cols[A_,i_]:=Module[{Wite=0.,Blacks},
	Blacks=Position[A[[i]],x_/;x!=Wite];
	If[Length[Blacks]>0,Flatten[{First[Blacks],Last[Blacks]}],{0,0}]]
	
Merge[A_,B_]:={Min[A[[1]],B[[1]]],Max[A[[2]],B[[2]]]}

BandFinder[A_]:=Module[{Col={},Ro={},i,starti=1,colrange,newcols,Newblock},
	colrange=Cols[A,starti];
	For[i=2,i<= First[Dimensions[A]],i++,newcols=Cols[A,i];
		colrange=If[colrange[[2]]<newcols[[1]]||newcols[[2]]<colrange[[1]],
			Newblock=True;AppendTo[Col,colrange];newcols,
			Newblock=False;Merge[newcols,colrange]];
		If[Newblock,AppendTo[Ro,{starti,i-1}];starti=i;,Null]];
	AppendTo[Col,colrange];
	AppendTo[Ro,{starti,i-1}];{Ro,Col}]

Blockdim[BlockMat_,i_]:={BlockMat[[1,i,2]]-BlockMat[[1,i,1]]+1,BlockMat[[2,i,2]]-BlockMat[[2,i,1]]+1}

BlockMerge[BlockMat_,i_]:=Module[{short},
	short=Delete[BlockMat,{{1,i+1},{2,i+1}}];
	ReplacePart[short,{
		{1,i}->{Min[BlockMat[[1,i,1]],BlockMat[[1,i+1,1]]],Max[BlockMat[[1,i,2]],BlockMat[[1,i+1,2]]]},
		{2,i}->{Min[BlockMat[[2,i,1]],BlockMat[[2,i+1,1]]],Max[BlockMat[[2,i,2]],BlockMat[[2,i+1,2]]]}}]]

NoDuplicates[seq_]:=Module[{lis,checklist,i},
	lis=seq;checklist=Range[Max[lis]];
	For[i=1,i<= Length[lis],i++,
		If[Position[checklist,lis[[i]]]!={},
			checklist=Delete[checklist,Position[checklist,lis[[i]]]],
			lis=Delete[lis,i];i--]];
	lis]

RowBlacken[rowvec_] := 
 Module[{nonzero, row}, row = rowvec; 
  nonzero = Flatten@Position[row, a_ /; a != 0]; 
  row[[First[nonzero] ;; Last[nonzero]]] = 1; row]
  
Blacken[Mat_] := Module[{BlockedA},
  BlockedA = Normal@Chop[Mat];
  BlockedA = Map[RowBlacken, Normal[BlockedA]];
  BlockedA = Transpose[Map[RowBlacken, Transpose[BlockedA]]]]
  
BlockFinder[A_]:=Module[{ComBlockCount=0,LastComBlockCount=1,VBlocks1,VBlocks2,VBlocks,
	HBlocks1,HBlocks2,HBlocks,ord,CommonBlocks,Remainder,ComBlx,RemBlx,BlockedA,
	CrossColseq,CrossRowseq,NewColseq,NewRowseq,RemBlock,CrossBlocks,ColCount,RowCount},
	BlockedA=Blacken[A];
(*	BlockedA=Map[RowBlacken, Normal[BlockedA]]; 
	BlockedA=Transpose[Map[RowBlacken, Transpose[BlockedA]]];*)
	{RowCount,ColCount}=Dimensions[BlockedA];
	NewColseq=Range[ColCount];NewRowseq=Range[RowCount];
	While[ComBlockCount!=LastComBlockCount,LastComBlockCount=ComBlockCount;
		VBlocks1=BandFinder[Transpose[BlockedA]];
		VBlocks2=BandFinder[Reverse[Transpose[BlockedA]]];
		VBlocks2=Reverse[VBlocks2,2];VBlocks2=ReplacePart[VBlocks2,1->Reverse[ColCount-VBlocks2[[1]]+1,2]];
		VBlocks=If[Length[VBlocks1[[1]]]>Length[VBlocks2[[1]]],VBlocks1,VBlocks2];
		HBlocks1=BandFinder[BlockedA];
		HBlocks2=BandFinder[Reverse[BlockedA]];
		HBlocks2=Reverse[HBlocks2,2];HBlocks2=ReplacePart[HBlocks2,1->Reverse[RowCount-HBlocks2[[1]]+1,2]];
		HBlocks=If[Length[HBlocks1[[1]]]>Length[HBlocks2[[1]]],Reverse[HBlocks1],Reverse[HBlocks2]];
		ord=Ordering[HBlocks[[1]]]; HBlocks={HBlocks[[1]][[ord]],HBlocks[[2]][[ord]]};
		Module[{BlockDif,i=1,rows,gap},
			While[HBlocks!= VBlocks,BlockDif=Blockdim[VBlocks,i]-Blockdim[HBlocks,i];
				Which[BlockDif== {0,0},i++,
					Max[BlockDif]>0,If[i==Length[HBlocks[[1]]],Break[]];
						rows=Sort[{HBlocks[[2,i]],HBlocks[[2,i+1]]}];
						gap=rows[[2,1]]-rows[[1,2]];
						If[gap<2,HBlocks=BlockMerge[HBlocks,i],Break[]],
					Min[BlockDif]<0,If[i==Length[VBlocks[[1]]],Break[]];
						VBlocks=BlockMerge[VBlocks,i]]]];
		ComBlx=Transpose[VBlocks]\[Intersection]Transpose[HBlocks];
		CommonBlocks=If[ComBlx== {},{{}},Transpose[ComBlx]];
		ComBlockCount=Length[CommonBlocks[[1]]];
		Remainder=If[Length[VBlocks[[1]]]>Length[HBlocks[[1]]],
			RemBlx=Complement[Transpose[VBlocks],ComBlx];If[RemBlx== {},{{}},Transpose[RemBlx]],
			RemBlx=Complement[Transpose[HBlocks],ComBlx];If[RemBlx== {},{{}},Transpose[RemBlx]]];
		CrossBlocks=Join[CommonBlocks,Remainder,2];
		{CrossColseq,CrossRowseq}=Map[NoDuplicates,Apply[Join,Apply[Range,CrossBlocks,{2}],{1}],2];
		BlockedA=BlockedA[[CrossRowseq, CrossColseq]];
		NewColseq=NewColseq[[CrossColseq]];
		NewRowseq=NewRowseq[[CrossRowseq]];
		];
	BlockedA=A[[NewRowseq, NewColseq]];
	RemBlock=If[Remainder== {{}},{{}},
				{{{Min[Remainder[[1]]],ColCount}},{{Min[Remainder[[2]]],RowCount}}}];
	CrossBlocks=Join[CommonBlocks,RemBlock,2];

	{BlockedA,NewRowseq,NewColseq,CrossBlocks}]

End[] (* End Private Context *)

EndPackage[]