(* Mathematica Package *)

(* COPYRIGHT (C) 2010 Wynand Verwoerd 
under the terms of the GNU General Public License V3
This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
*)

BeginPackage["RandomWalkBlocker`",{"NetSplitter`","MatrixDisplay`","BlockFinder`","FileReadWrite`","Printout`"}]
(* Exported symbols added here with SymbolName::usage *)  
FillOut::usage="Performs one-shot blocking with all externals, does incorporation and collection of supersinks "
BlockMaker::usage="Control routine that calculates DAG, hierarchical grouping and scoring of blocks"
Combine::usage="Function to merge selected blocks and orphans"
ExtendedBlock::usage="Function to incorporate non-bridging externals, and orphans into blocks"
HierarchGroup::usage="Rearranges a matrix with rows and columns hierarchically grouped"
Concertina::usage="Utility to spread out duplicates in a numerical list so Bincounts don't eliminate them"

Begin["`Private`"] (* Begin Private Context *) 

Probab[A_]:=N[Thread[Quiet[Normalize[Normal[A],Total[#]&]]]]

Expanse[A_]:=A.A

Residual[A_]:=Norm[A,"Frobenius"]

NodeZapper[nodelist_]:=(CR=Delete[CR,Split[nodelist]];RC=Transpose[Delete[Transpose[RC],Split[nodelist]]];InternalComplist=Delete[InternalComplist,Split[nodelist]];InternalCNodes=InternalCNodes-Length[nodelist];)

Concertina[folds_, delta_] := 
 Module[{i = 1, itop, j, fold}, itop = Length[folds]; fold = folds; 
  While[i < itop, j = 1; 
   While[i + j <= itop && fold[[i + j]] == fold[[i]], 
    fold[[i + j]] += j*delta; j++]; i = i + j]; fold]

FillOut[SuckIn_]:=Module[{Internals,CRP,RCP,Prob,count,residue,change=1,CB1,CB2,LIB,TotBlocks,blockno,
	domina,domination,i,bot,top},
Status="Weeding out surplus externals";
FinalBlockCount=Length[FinalBlocks];
Internals=Flatten[FinalBlocks];
Internalrows=Flatten[Table[Position[Complist,Internals[[i]]],{i,1,Length[Internals]}]];

(* externals gathered interactively are now applied to all the blocks simultaneously*)
If[SuckIn,
	Do[FinalBlocks[[i]]=ExtendedBlock[FinalBlocks[[i]]];
		FinalBlocks[[FinalBlockCount]]=IsolationBlock;
		Internals=Flatten[FinalBlocks];
		Internalrows=Flatten[Table[Position[Complist,Internals[[i]]],{i,1,Length[Internals]}]],{i,1,FinalBlockCount-1}];
		If[ Iterationcount>0, (*Repeat blocking procedure with final set of externals*)
			Status="Apply collected externals to all blocks ";IsolationBlock={};
			Externals=Complement[Complist,Internals];
			Externals=Union[Externals,BasicExternals];
			Externalrows=Flatten[Table[Position[Complist,Externals[[i]]],{i,1,Length[Externals]}]];
			CR=Abs[Sign[S]]-Sign[S];RC=Transpose[Abs[Sign[S]]+Sign[S]];
			If[ReverseThem,
				RC=Join[RC,Transpose@CR[[All,Reversibles]]];
				CR=Transpose@Join[Transpose[CR],RC[[Reversibles]]]];
			InternalComplist=Complist;InternalCNodes=CNodes;
			BlockMaker[Externalrows];
			IRB=Transpose[MapAt[#+Sinks&,CrossBlocks,2]];
			IRBlocks=Apply[Range,IRB,{2}];IRBlocks=Map[Flatten,IRBlocks];
			IRBlocks=Table[Nameseq[[IRBlocks[[i]]]],{i,1,Length[IRBlocks]}];
			FinalBlocks=Join[IRBlocks,{IsolationBlock}];
			FinalBlockCount=Length[FinalBlocks];
			Internals=Flatten[FinalBlocks];
			Internalrows=Flatten[Table[Position[Complist,Internals[[i]]],{i,1,Length[Internals]}]];
			Do[FinalBlocks[[i]]=ExtendedBlock[FinalBlocks[[i]]];
				FinalBlocks[[FinalBlockCount]]=IsolationBlock;
				Internals=Flatten[FinalBlocks];
				Internalrows=Flatten[Table[Position[Complist,Internals[[i]]],{i,1,Length[Internals]}]],{i,1,FinalBlockCount-1}];
		];
	];

(* recalculate the DAG according to new block allocation*)
Externals=Complement[Complist,Internals];
Externalrows=Complement[Range[CNodes],Internalrows];
CR=Abs[Sign[S]]-Sign[S];RC=Transpose[Abs[Sign[S]]+Sign[S]];
If[ReverseThem,
	RC=Join[RC,Transpose@CR[[All,Reversibles]]];CR=Transpose@Join[Transpose[CR],RC[[Reversibles]]]];
InternalComplist=Complist;InternalCNodes=CNodes;
NodeZapper[Externalrows];
CRP=Probab[CR]; RCP=Probab[RC];
Prob=Probab[(CRP.RCP)+0.25*IdentityMatrix[InternalCNodes]];
DAG=Prob; count=1;residue=Residual[Prob];change=1;
While[Abs[change]>10^(-10),DAG=Chop[DAG.DAG];oldres=residue; residue=Residual[DAG];change=residue-oldres;count++];
Blockseq=Flatten@Table[Position[InternalComplist,Internals[[i]]],{i,1,Length[Internals]}];
BlockDAG=Transpose[Transpose[DAG][[Blockseq]]][[Blockseq]];
Transbin=Transpose[Sign[BlockDAG]];
Sources=Count[Transbin,Table[0,{InternalCNodes}]];
Sinks=InternalCNodes-Sources;
Collectseq=Join[Nonzerows[Transbin],Complement[Range[InternalCNodes],Nonzerows[Transbin]]];
Shuffleseq=Blockseq[[Collectseq]];Nameseq=InternalComplist[[Shuffleseq]];
Sinkseq=Nameseq[[1;;Sinks]];
Sourceseq=Nameseq[[-Sources;;]];
ShuffleDAG=Transpose[Transpose[DAG][[Shuffleseq]]][[Shuffleseq]][[All,1;;Sinks]];
BlockCount=FinalBlockCount-1;
BlockSinks=Table[Length[Sinkseq\[Intersection]FinalBlocks[[i]]],{i,1,BlockCount}];

CB1=Module[{topp=0},
	Reap[Do[Sow[{topp+1,topp+=BlockSinks[[i]]}],{i,1,BlockCount}]][[2,1]]];
CB2=Module[{topp=0,incr},
	Reap[Do[incr=Length[FinalBlocks[[i]]]-BlockSinks[[i]];
		Sow[{topp+If[incr>0,1,0],If[incr>0,topp+=incr,topp-1]}],{i,1,BlockCount}]][[2,1]]];
CrossBlocks={CB1,CB2};
top=CB1[[BlockCount,2]];
LIB=Length[IsolationBlock];
If[LIB>0,AppendTo[CB1,{top+1,top+LIB}];
	TotBlocks=FinalBlockCount,TotBlocks=FinalBlockCount-1];

(*Shuffle the sinks so that supersinks in each block are grouped together*)
Shuffleseq=
Join[Flatten@Table[bot=CB1[[blockno,1]];top=CB1[[blockno,2]];
	Which[top-bot>0,SinkblockBin=Sign[ShuffleDAG[[bot;;top,bot;;top]]];
		Grouper[SinkblockBin];GroupColseq+bot-1,top-bot==0,{bot},top-bot<0,{}],{blockno,1,TotBlocks}],
	Range[Sinks+1,Sinks+Sources]];
Nameseq=Nameseq[[Shuffleseq]];
Sinkseq=Nameseq[[1;;Sinks]];
ShuffleDAG=Transpose[Transpose[ShuffleDAG][[Shuffleseq[[1;;Sinks]]]]][[Shuffleseq]];

(* Next, orphan block supersinks are reclassified as small blocks.  *)
If[LIB>0,
OrphanSinklimits=BandFinder[ShuffleDAG[[Sinks-LIB+1;;Sinks,Sinks-LIB+1;;Sinks]]][[1]];
SuperFirst=Map[Apply[Range,#]&,Sort[OrphanSinklimits,#2[[2]]-#2[[1]]<#1[[2]]-#1[[1]]&]];
If[Length[SuperFirst[[1]]]>1,Shuffleseq=Join[Range[Sinks-LIB],
	Sinks-LIB+Flatten[SuperFirst],Range[Sinks+1,Sinks+Sources]];
	Nameseq=Nameseq[[Shuffleseq]];Sinkseq=Nameseq[[1;;Sinks]];
	ShuffleDAG=Transpose[Transpose[ShuffleDAG][[Shuffleseq[[1;;Sinks]]]]][[Shuffleseq]];
	SuperSinks=Length/@Select[SuperFirst,Length[#]>1&];LSuper=Length[SuperSinks];
	CB1Super=Sinks-LIB+
	Module[{topp=0},
		Reap[Do[Sow[{topp+1,topp+=SuperSinks[[i]]}],{i,1,Length[SuperSinks]}]][[2,1]]];
	CB2Super=Array[{Sources,Sources-1}&,LSuper];
	CB1=Join[Most[CB1],CB1Super];CB2=Join[CB2,CB2Super];
	CrossBlocks={CB1,CB2};
	LIB=LIB-Total[SuperSinks];
	If[LIB>0,AppendTo[CB1,{Sinks-LIB+1,Sinks}]];
	FinalBlockCount=Length[CB1];BlockCount=Length[CB2]]];

FinalBlocks=Table[Join[Nameseq[[CrossBlocks[[1,blockno,1]];;CrossBlocks[[1,blockno,2]]]],
						Nameseq[[CrossBlocks[[2,blockno,1]]+Sinks;;CrossBlocks[[2,blockno,2]]+Sinks]]],
					{blockno,1,BlockCount}];

IsolationBlock=If[LIB>0,Sinkseq[[-LIB;;]],{}];
AppendTo[FinalBlocks,IsolationBlock];FinalBlockCount=Length[FinalBlocks];
IntrowsperBlock=Table[Flatten@Table[Position[Complist,FinalBlocks[[i,j]]],
					{j,1,Length[FinalBlocks[[i]]]}],{i,1,FinalBlockCount}];
ColsperBlock=Quiet[Table[Nonzerows[Normal[Transpose[S[[IntrowsperBlock[[i]]]]]]], {i, 1, 
     FinalBlockCount}], {Transpose::"nmtx"}];
     
(* Identify the dominating compartment for each block *)
domina = Map[Commonest[Compartlist[[#]]] &, IntrowsperBlock];
domination = MapThread[If[#1 != {}, 
	Count[Compartlist[[#2]], #1[[1]]]/Length[#2],0] &, {domina, IntrowsperBlock}];
BlockCompart= MapThread[
	Which[#2 > 0.6666, #1[[1]], #2 == 0, "None", True,"Mixed"] &, {domina, domination}];
BlockHue=ConstantArray[0,FinalBlockCount]; BlockBands[Sources]; 

]

Get["HierarchicalClustering`"];
Off[Agglomerate::"ties"];

CFlat[clust_,count_]:=Replace[Thread[ClusterFlatten[ClusterSplit[clust,count]]],x_/;Not[VectorQ[x,IntegerQ]]:>{x[[1]]},1]

DistanceList[Clust_]:=Sort[Extract[Clust,Map[Append[#,3]&,Position[Clust,_Cluster]]],Greater]

HierarchGroup[Matrix_, DistChoice_: SokalSneathDissimilarity, LinkChoice_: "Single"] := 
 Module[{Rowdim, Coldim, BlockRowseq, BlockColseq, Mat}, 
   {Rowdim, Coldim} = Dimensions[Matrix]; Mat=Normal[Matrix];
  BlockRowseq = 
    ClusterFlatten[Agglomerate[Mat -> Range[Rowdim],
    	DistanceFunction -> DistChoice,Linkage -> LinkChoice]];
  BlockColseq = 
    ClusterFlatten[Agglomerate[Transpose[Mat] -> Range[Coldim], 
      	DistanceFunction -> DistChoice, Linkage -> LinkChoice]];
   {Matrix[[BlockRowseq, BlockColseq]],BlockRowseq, BlockColseq}]

Grouper[ThinBin_]:=
  (* Specialized version of HierarchGroup, takes only a binary matrix ThinBin, 
     hierarchy parameters taken from global variables and returns only the matrix itself *)
Module[{Rowdim,Coldim},
	{Rowdim,Coldim}=Dimensions[ThinBin];
	(RowHierarchy=Agglomerate[ThinBin->Range[Rowdim],DistanceFunction->DistChoice,Linkage->LinkChoice];
	GroupRowseq=ClusterFlatten[RowHierarchy];
	ColHierarchy=Agglomerate[Transpose[ThinBin]->Range[Coldim],DistanceFunction->DistChoice,Linkage->LinkChoice];
	GroupColseq=ClusterFlatten[ColHierarchy];
	Transpose[Transpose[ThinBin][[GroupColseq]]][[GroupRowseq]])]

Blocker[GroupMat_,Rob_,Cob_]:=Module[{RowGroups,ColGroups,RowTops,ColTops,Rowtot,Coltot,Blocks,Blocksize},
	RowGroups=CFlat[RowHierarchy,Rob];ColGroups=CFlat[ColHierarchy,Cob];
	RowTops=Map[Length,RowGroups];PrependTo[RowTops,0];
	Do[RowTops[[i]]=RowTops[[i]]+RowTops[[i-1]],{i,2,Length[RowTops]}];
	ColTops=Map[Length,ColGroups];PrependTo[ColTops,0];
	Do[ColTops[[i]]=ColTops[[i]]+ColTops[[i-1]],{i,2,Length[ColTops]}];
	Blocks=Length[ColTops];
	Rowtot=Replace[Total[Transpose[GroupMat]],0->1,1];
	Blocksize=Table[ColTops[[k+1]]-ColTops[[k]],{k,1,Blocks-1}];
	RowBlocking=Table[Flatten[Table[ConstantArray[N[Total[GroupMat[[i,ColTops[[k]]+1;;ColTops[[k+1]]]]]/Rowtot[[i]]],Blocksize[[k]]],
		{k,1,Blocks-1}]],{i,1,Sources}];
	Blocks=Length[RowTops];
	Coltot=Replace[Total[GroupMat],0->1,1];Blocksize=Table[RowTops[[k+1]]-RowTops[[k]],{k,1,Blocks-1}];
	ColBlocking=Transpose[Table[Flatten[Table[ConstantArray[N[Total[GroupMat[[RowTops[[k]]+1;;RowTops[[k+1]],i]]]/Coltot[[i]]],Blocksize[[k]]],
		{k,1,Blocks-1}]],{i,1,Sinks}]];
	CombinedBlocking=(RowBlocking+ColBlocking)/2;
	whites=Count[CombinedBlocking,x_/;x==0,2];

((Total[CombinedBlocking,2])/(Sinks*Sources-whites))^2*whites/(Sinks*Sources)
]

BlackScore[GLevel_]:=
(RowBlocks=Count[DistanceList[RowHierarchy],x_/;x>= GLevel,2]+1;
 ColBlocks=Count[DistanceList[ColHierarchy],x_/;x>= GLevel,2]+1;
Blocker[GroupBin,RowBlocks,ColBlocks])

BlockMaker[Externalrows_]:=Module[{CRP,RCP,Prob,count,residue,change,oldres,Transbin,Collectseq,Collectcols,Empties,ThinBin,
									DistanceValues,TrialLevels,BlockRowseq,BlockColseq,RowBlocking2,ColBlocking2,CrossColseq,CrossRowseq},
NodeZapper[Externalrows];
If[InternalCNodes>0,CRP=Probab[CR]; RCP=Probab[RC];
	Prob=Probab[(CRP.RCP)+0.25*IdentityMatrix[InternalCNodes]];DAG=Prob; count=1;
	residue=Residual[Prob];change=1;
	While[Abs[change]>10^(-10),DAG=Chop[DAG.DAG];oldres=residue; residue=Residual[DAG];change=residue-oldres;count++],
	DAG={{1.0}}];
Transbin=Transpose[Sign[DAG]];Sources=Count[Transbin,Table[0,{InternalCNodes}]];
Sinks=InternalCNodes-Sources;
Collectseq=Join[Nonzerows[Transbin],Complement[Range[InternalCNodes],Nonzerows[Transbin]]];

If[Sinks<2 || Sources<2,
	EmitSound[Crisis];
	ChoiceDialog["No blocks to process, as the DAG has only one row or column",{" OK "->1},WindowFloating->True,ButtonBoxOptions->{ImageMargins->6}];
	FinalBlocks={InternalComplist,{}};BlockCount=1;FinalBlockCount=2;Shuffleseq=Collectseq;
	ShuffleDAG=Transpose[Transpose[DAG][[Shuffleseq]]][[Shuffleseq]][[All,1;;Sinks]];
	Nameseq=InternalComplist[[Shuffleseq]];BlockCompart={DefaultCompartment,DefaultCompartment};
	CrossBlocks={{{1,Sinks}},{{1,Sources}}};Return[False]];

Collectcols=Transpose[Transbin[[Collectseq[[1;;Sinks]]]]][[Collectseq]];
Empties=Position[Total[Collectcols[[-Sources;;]]],0];
NewExternals=InternalComplist[[Collectseq[[Flatten[Empties]]]]];(*Orphans=NewExternals;*)
IsolationBlock=Union[IsolationBlock,NewExternals];
Collectseq=Delete[Collectseq,Empties];
Sinks=Sinks-Length[Empties];
If[Sinks<2,
	EmitSound[Crisis];
	ChoiceDialog["No blocks to process, as the DAG has only one row or column",{" OK "->1},WindowFloating->True,ButtonBoxOptions->{ImageMargins->6}];
	FinalBlocks={InternalComplist,{}};BlockCount=1;FinalBlockCount=2;Shuffleseq=Collectseq;
	ShuffleDAG=Transpose[Transpose[DAG][[Shuffleseq]]][[Shuffleseq]][[All,1;;Sinks]];
	Nameseq=InternalComplist[[Shuffleseq]];BlockCompart={DefaultCompartment,DefaultCompartment};
	CrossBlocks={{{1,Sinks}},{{1,Sources}}};Return[False]];
	
Collectcols=Transpose[Transbin[[Collectseq[[1;;Sinks]]]]][[Collectseq]];
ThinBin=Collectcols[[-Sources;;,1;;Sinks]];
Status="Setting up hierarchical grouping";
GroupBin=Grouper[ThinBin];DistanceValues=Union[N[DistanceList[RowHierarchy]],N[DistanceList[ColHierarchy]]];
TrialLevels=Cases[(Most[DistanceValues]+Rest[DistanceValues])/2,x_/;GroupLevelRange[[1]]<x<GroupLevelRange[[2]]];
(*Even in a nontrivial DAG, it can happen that all rows/columns are so similar that no groups could be distinguished; in this case abort the blocking attempt*)
If[TrialLevels== {},FinalBlocks={InternalComplist,{}};BlockCount=1;FinalBlockCount=2;Shuffleseq=Collectseq;
	ShuffleDAG=Transpose[Transpose[DAG][[Shuffleseq]]][[Shuffleseq]][[All,1;;Sinks]];
	Nameseq=InternalComplist[[Shuffleseq]];CrossBlocks={{{1,Sinks}},{{1,Sources}}};Return[False]];

GroupLevel=TrialLevels[[Ordering[Map[BlackScore,TrialLevels],-1]]][[1]];
BlackScore[GroupLevel];Status="Set up blocking and rearrange accordingly";
BlockRowseq=ClusterFlatten[Agglomerate[CombinedBlocking->Range[Sources],DistanceFunction->BrayCurtisDistance,Linkage->LinkChoice]];
BlockColseq=ClusterFlatten[Agglomerate[Transpose[CombinedBlocking]->Range[Sinks],DistanceFunction->BrayCurtisDistance,Linkage->LinkChoice]];
RowBlocking2=Transpose[Transpose[RowBlocking][[BlockColseq]]][[BlockRowseq]];
ColBlocking2=Transpose[Transpose[ColBlocking][[BlockColseq]]][[BlockRowseq]];
ConsolBlocking=(RowBlocking2+ColBlocking2)/2;
Status="Scan for blocks and sort them diagonally";
{ConsolBlocking,CrossRowseq,CrossColseq,CrossBlocks}=BlockFinder[ConsolBlocking];
BlockCount=Length[CrossBlocks[[1]]];
BlockColseq=BlockColseq[[CrossColseq]];BlockRowseq=BlockRowseq[[CrossRowseq]];
ConsolBin=Transpose[Transpose[GroupBin][[BlockColseq]]][[BlockRowseq]];Blockseq=Join[GroupColseq[[BlockColseq]],Sinks+GroupRowseq[[BlockRowseq]]];
Shuffleseq=Collectseq[[Blockseq]];
ShuffleDAG=Transpose[Transpose[DAG][[Shuffleseq]]][[Shuffleseq]][[All,1;;Sinks]];Nameseq=InternalComplist[[Shuffleseq]];
Sinkseq=Nameseq[[1;;Sinks]];Sourceseq=Nameseq[[-Sources;;]];True]


MergeDialog[step_]:=(Module[{Chosen={},Choice,compnos={},BlockHeaders,Tinge,Orphanstyles,Textarea,wide,aspect,
							Vetter,Matview,Targets,ChosenCompart,CompartSetter,MergeControls,Carryon},
	SelectedBlocks={};BlockHeaders=Table["Block "<>ToString[i]<>"\n"<>BlockCompart[[i]],{i,1,FinalBlockCount}];
	BlockHeaders=Pick[BlockHeaders,Map[Length[#]>0&,FinalBlocks]];
	BlockHeaders=If[IsolationBlock=={},{BlockHeaders},{Most[BlockHeaders]}];
	If[Length[IsolationBlock]>0,compnos=Map[Position[Nameseq,#][[1,1]]&,IsolationBlock,{1}];
		Orphanstyles=MapThread[(If[MemberQ[TargetCompounds,#2],Tinge=Targetinge,Tinge=None];
			Style[ToString[#1]<>"  "<>#2,Background->Tinge])&,{compnos,IsolationBlock}],
		Orphanstyles={}];
	Targets=Sort@Flatten[Table[Position[Nameseq,TargetCompounds[[i]]],{i,1,Length[TargetCompounds]}]];
	Textarea = Rasterize[TextCell[Style[TableForm[PrintBlocks[SelectedBlocks], 
      TableSpacing -> {4, 5, 0}, TableHeadings -> BlockHeaders],FontSize -> 10, Bold], "Output"]];
	wide = ImageDimensions[Textarea][[1]]; aspect = Max[ImageAspectRatio[Textarea], 1];
	ChosenCompart = "Not Selected"; 
	CompartSetter = RadioButtonBar[Dynamic[ChosenCompart, 
		(SelectedBlocks = Flatten@Position[BlockCompart[[1;;BlockCount]], #];
		Chosen = Choice[[Flatten@Position[Compartlist[[IntrowsperBlock[[-1]]]], #]]];
		ChosenCompart = #) &],
		Prepend[CompartIDs, "Deselect all"], Appearance -> "Horizontal"];
	(*Activate sections referring to Choice to present isolated compounds as tick boxes. 
	  This can become slow with large networks that have many orphans*)
	Choice=Flatten[Table[Position[Nameseq,IsolationBlock[[i]]],{i,1,Length[IsolationBlock]}]];
	If[Choice!= {},Vetter=CheckboxBar[Dynamic[Chosen],
		Table[Choice[[i]]->ToString[Choice[[i]]]~~"  "~~Nameseq[[Choice]][[i]]~~"\t",{i,1,Length[Choice]}],
		Appearance->"Row",Spacings->2]];
		
  MergeControls = Column[{
	 Row[{
		Panel[Grid[{{Tooltip[Style["Selected Blocks \n", "Subsubsection"], 
			"Select by mouseclick, or type below in the format {2,4,5}"]},
			{InputField[Dynamic[SelectedBlocks,(SelectedBlocks = Select[#, 0 < # <= BlockCount &]) &], FieldSize -> 40],
			SpanFromLeft}}]], 
		Panel[Grid[{
      	{"Efficacy (" <> ToString[Eff[Effpower, {}]]<>"%)",	ProgressIndicator[Eff[Effpower, {}]/100]},
      	{"Merge and display ", DefaultButton["Continue", DialogReturn[True]]},
      	{"Merge and exit", CancelButton["Finish", DialogReturn[False]]}}],Background->Lighter[Magenta,0.75]]}, 
     Spacer[30]],
     If[Length[CompartIDs] > 1, 
		Panel[Grid[{{Tooltip[Style["Compartment Selection \n", "Subsubsection"],
       		"Selects all blocks and orphans associated with a compartment"]}, {CompartSetter, SpanFromLeft}}]]],
     If[Choice != {},
          Panel[Grid[{{Tooltip[Style["Selected Orphans \n", "Subsubsection"],
     	"Select individually to merge with chosen blocks"]},{Vetter, SpanFromLeft}}],ImageSize -> wide]]
     	
     	(*Activate the next If statement instead of previous to a give text box where orphans are selected by number*)
	(*If[compnos!={},Panel[Grid[{{Tooltip[Style["Selected orphans \n","Subsubsection"],
    	"Type orphan numbers in the format {2,4,5} to merge individually with chosen blocks"]},
    	{InputField[Dynamic[Chosen,(Chosen=Select[#,MemberQ[compnos,#]&])&],FieldSize->30],SpanFromLeft}}]]],*)
    	}, Spacings->2];

  Matview = Module[{xbin, ybin, ClickBlock},
		xbin = Table[CrossBlocks[[1, i, 2]], {i, 1, BlockCount}]; PrependTo[xbin, 0];
   		xbin = Concertina[xbin, 0.0001];
   		ybin = Table[CrossBlocks[[2, i, 2]], {i, 1, BlockCount}]; PrependTo[ybin, 0];
   		Do[If[ybin[[i]] < ybin[[i - 1]], ybin[[i]] += 1; ybin[[i - 1]] -= 0.5;], {i, 2, Length[ybin]}];
		ybin = Concertina[ybin, 0.0001];
	Panel@ClickPane[Dynamic@Refresh[Show[{
		MatrixPlot[BlockBands[Sinks],ColorFunction -> CompartHue, ColorFunctionScaling -> False, 
          	ImageSize -> wide, AspectRatio -> aspect, FrameTicks -> TickSpec[{}, 0]], 
		MatrixPlot[PinkBlue[Sinks, Targets, SelectedBlocks], 
			ColorRules -> {0 -> None, 2 -> Targetinge, 4 -> Foundtinge, 6 -> Chosentinge}, 
			PlotLabel -> Style["Blocked DAG matrix", "Subsection"]], 
		MatrixPlot[Map[If[# == 0, None, #] &, ShuffleDAG, {2}], ColorRules -> {0 -> None}]}, 
	  Background -> White], 
      TrackedSymbols -> SelectedBlocks],
      	(ClickBlock = Flatten@Position[Diagonal@BinCounts[{{#[[1]], Sources - #[[2]]}}, {xbin}, {ybin}], 1];
	  SelectedBlocks = BlockSelection[ClickBlock]; 
      ChoosO = (SelectedBlocks != {})) &]];
       
	CalculationTime+=TimeUsed[]-StartTime;EmitSound[Attention];
	Carryon=DialogInput[TabView[{
		"Selector" -> MergeControls,
		"Matrix View" -> Matview,
		"Block listing" -> Dynamic[Column[{
			TextCell["Internal metabolites in each block ", "Subtitle"], 
			TextCell[Style[TableForm[PrintBlocks[SelectedBlocks], TableSpacing -> {4, 5, 0}, TableHeadings -> BlockHeaders], 
  					FontSize -> 10, Bold], "Output"]}, 
  			Background -> LightYellow]], 
  		"Orphans" -> Dynamic[Column[{
  			TextCell["Isolated sinks", "Subtitle"], 
  			TextCell["---------------------"], 
  			TextCell[Style[TableForm[Partition[Orphanstyles, 3, 3, 1, {}], 
  				TableSpacing -> {1, 5}], FontSize -> 10, Bold], "Output"]}, 
  			Background -> LightPink]], 
  		"Targets" -> Dynamic[Column[{
  			TextCell["Compounds of special interest","Subtitle"], 
  			TextCell["-----------------------------------------------"], 
  			TextCell[Style[TableForm[Partition[TargetCompounds, 3, 3, 1, {}]],FontSize -> 10, Bold], "Output"]}, 
		Background -> Targetinge]]}],
		WindowTitle->"Blocking display and merging: Step "<>ToString[step],
		WindowElements->{"HorizontalScrollBar","VerticalScrollBar","MagnificationPopUp"},
		WindowSize->{All,All},WindowFrameElements->{"ZoomBox"},
		WindowMargins->{{Automatic,Automatic},{10,0}}];
	StartTime=TimeUsed[];
	{Chosen,Carryon}])

Combine[LongBlock_]:=(Module[{MergeOrphanrows,Carryon=True,MergeOrphans,MergeBlock,itemcount,newblock,selblock},
	If[Length[LongBlock]>1,
		selblock=LongBlock[[1]];MergeOrphanrows=LongBlock[[2]];
		If[RepeatMerge && (Complement[MergeOrphanrows,Range[Sinks-Length[IsolationBlock]+1,Sinks]]!= {}||Complement[selblock,Range[1,BlockCount]]!= {}),
			DialogInput[DialogNotebook[{TextCell["Invalid merge data read for step "<>ToString[Stepcount+1]<>"\nFurther imported merges abandoned"],Button["OK",DialogReturn[1]]}]];RepeatMerge=False;];
		If[RepeatMerge,Status="Implement merge history loaded from Externals file";
			MergeOrphans=Nameseq[[MergeOrphanrows]];MergeBlock=Union@@Append[FinalBlocks[[selblock]],MergeOrphans];
			If[Length[MergeBlock]>1,IsolationBlock=Complement[IsolationBlock,MergeOrphans];
				Do[FinalBlocks=ReplacePart[FinalBlocks,selblock[[i]]->{}],{i,1,Length[selblock]}];
				newblock=Length[FinalBlocks];FinalBlocks=ReplacePart[FinalBlocks,newblock->ExtendedBlock[MergeBlock]];
				AppendTo[FinalBlocks,IsolationBlock];
				AppendTo[MergeTrace,{++Stepcount,If[selblock== {},"None",selblock],If[MergeOrphanrows== {},"None",MergeOrphanrows],newblock}];
				FillOut[False];];PostMergeEff = Eff[Effpower, {}];Return[]],
		While[Carryon,Status="Waiting for merge instructions";
			{MergeOrphanrows,Carryon}=MergeDialog[Stepcount+1];
			itemcount=Length[SelectedBlocks]+Length[MergeOrphanrows];
			If[Carryon&&(itemcount<2),Beep[];
				DialogInput[DialogNotebook[{TextCell["At least two items must be selected for merging."],DefaultButton[]},WindowTitle->"Merging Canceled"]]];
			If[itemcount>1,	
				MergeOrphans=Nameseq[[MergeOrphanrows]];
				MergeBlock=Union@@Append[FinalBlocks[[SelectedBlocks]],MergeOrphans];
				If[Length[MergeBlock]>1,IsolationBlock=Complement[IsolationBlock,MergeOrphans];
					Do[FinalBlocks=ReplacePart[FinalBlocks,SelectedBlocks[[i]]->{}],{i,1,Length[SelectedBlocks]}];
					newblock=Length[FinalBlocks];
					FinalBlocks=ReplacePart[FinalBlocks,newblock->ExtendedBlock[MergeBlock]];
					AppendTo[FinalBlocks,IsolationBlock];
					AppendTo[MergeTrace,{++Stepcount,If[SelectedBlocks== {},"None",SelectedBlocks],If[MergeOrphanrows== {},"None",MergeOrphanrows],newblock}];
					FillOut[False];];
				SelectedBlocks={}];
			];
			PostMergeEff = Eff[Effpower, {}];
		]])

ExtendedBlock[Block_]:=Module[{WorkBlock=Block,BlockIntRows,BlockIntCols,BlockRows,BlockExtRows,BlockRowShortlist,OrphanRows,NonBlockIntRows,pp,ShortlistCols,Harvest,NewBlockIntRows,NewInts=10,ElectedFromShort,a,b,i,j,LinkingCols,NoLinks,Rehab,OrphansNomore},
	While[NewInts>0,
		BlockIntRows=Flatten[Table[Position[Complist,WorkBlock[[i]]],{i,1,Length[WorkBlock]}]];
		BlockIntCols=Nonzerows[Normal[Transpose[S[[BlockIntRows]]]]];
		BlockRows=Nonzerows[Normal[S[[All,BlockIntCols]]]];
		BlockExtRows=Complement[BlockRows,BlockIntRows];
		BlockRowShortlist=Complement[BlockExtRows,StructuralExternalrows\[Union]StExrows];
		OrphanRows=Flatten[Table[Position[Complist,IsolationBlock[[i]]],{i,1,Length[IsolationBlock]}]];
		NonBlockIntRows=Complement[Internalrows,BlockIntRows\[Union]OrphanRows];
		pp=Split[Position[Normal@S[[BlockRowShortlist,All]],x_/;x!=0],#1[[1]]== #2[[1]]&];
		ShortlistCols=Table[pp[[i,j,2]],{i,1,Length[pp]},{j,1,Length[pp[[i]]]}];
		Harvest=If[Length[NonBlockIntRows]>0&&Length[ShortlistCols]>0,
			Reap@Do[If[MatrixQ[Transpose@Normal@S[[NonBlockIntRows,ShortlistCols[[i]]]],#==0&],Sow[i,a];Sow[BlockRowShortlist[[i]]],b];,{i,1,Length[BlockRowShortlist]}],{}];
		ElectedFromShort=If[Depth[Harvest]== 4,Harvest[[2,1]],{}];
		NewBlockIntRows=If[Depth[Harvest]== 4,Harvest[[2,2]],{}];
		NewInts=Length[NewBlockIntRows];
		If[IsolationBlock!= {},OrphanRows=Flatten[Table[Position[Complist,IsolationBlock[[i]]],{i,1,Length[IsolationBlock]}]];
			LinkingCols=ShortlistCols[[ElectedFromShort]];
			NoLinks=Length[LinkingCols];
			Rehab=If[NoLinks>0,Total[Abs@Table[Transpose@Normal@S[[OrphanRows,LinkingCols[[i]]]],{i,1,NoLinks}],2],{}];
			OrphansNomore=If[Total[Rehab]>0,Pick[OrphanRows,Rehab,x_/;x!= 0],{}];
			IsolationBlock=Complement[IsolationBlock,Complist[[OrphansNomore]]];
			Internalrows=Complement[Internalrows,OrphansNomore]];
		WorkBlock=Join[WorkBlock,Complist[[NewBlockIntRows]]]];
		WorkBlock]


End[] (* End Private Context *)

EndPackage[]