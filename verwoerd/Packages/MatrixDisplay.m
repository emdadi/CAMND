(* Mathematica Package *)

(* COPYRIGHT (C) 2010 Wynand Verwoerd 
under the terms of the GNU General Public License V3
This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
*)

BeginPackage["MatrixDisplay`",{"NetSplitter`"}]
(* Exported symbols added here with SymbolName::usage *)  
BlockSelection::usage="Routine to handle interactive selection of a block by clicking on the matrix plot"
PinkBlue::usage="Produces coloured background blocks on matrix plots"
Candidates::usage="Finds candidate metabolites to be made external"
BlockBands::usage="Displays a background band for each block, coloured according to cellular compartment"
CompartHue::usage="Colorfunction that asigns the compartment colour"
TickSpec::usage="Specify ticks on matrix plots"
Eff::usage="Calculates efficacy of a partitioning into blocks and orphans"

Begin["`Private`"] (* Begin Private Context *) 

BlockSelection[CBlock_]:=Complement[SelectedBlocks\[Union]CBlock,SelectedBlocks\[Intersection]CBlock]

BlockPlot[Offset_,ChooseBlocks_]:=Module[{BBlocks=ConstantArray[0,{Sources+Offset, Sinks}],Tinge,i},
	For[i=1,i<= BlockCount,i++,
		Tinge=If[MemberQ[ChooseBlocks,i],6,4];
		BBlocks[[CrossBlocks[[2,i,1]]+Offset;;Max[CrossBlocks[[2,i,2]],CrossBlocks[[2,i,1]]]+Offset,CrossBlocks[[1,i,1]];;CrossBlocks[[1,i,2]]]]=Tinge;];
BBlocks]

Pinkblock[Offset_,Targets_]:=Module[{rows,cols,i,TT=ConstantArray[0,{Sources+Offset, Sinks}]},
	If[Targets!= {},
		{cols,rows}=Which[First[Targets]>Sinks,{{},Targets},
			Last[Targets]<= Sinks,{Targets,{}},
			True,Split[Targets,!(#1<= Sinks&&#2>Sinks)&]];
			Do[TT[[1;;,cols[[i]]]]=2,{i,1,Length[cols]}];
			Do[TT[[rows[[i]]-Sinks+Offset,1;;]]=2,{i,1,Length[rows]}]];
TT];

PinkBlue[Ofset_,Targets_,ChooseBlocks_]:=
	If[Sinks>0 && Sources>0,MapThread[Max,{BlockPlot[Ofset,ChooseBlocks],Pinkblock[Ofset,Targets]},2],{{0}}];

Candidates[Chosen_]:=Module[{BB=ConstantArray[0.,{Sources,Sinks}],CompNo,i},
For[i=1,i<= Length[Chosen],i++,CompNo=Chosen[[i]];
	If[CompNo<= Sinks,BB[[All,CompNo]]=ConsolBlocking[[All,CompNo]],BB[[CompNo-Sinks]]=ConsolBlocking[[CompNo-Sinks]]]];
BB]

CompartHue[a_] := Hue[a, 0.15, 1]

BlockBands[Offset_] := 
 Module[{BBands = ConstantArray[None, {Sources + Offset, Sinks}], 
   sortedcomparts, colors, colorseq, middle, colorgap, loro, hiro, blockclist, commons, maincolor, meancolor},
  commons = Map [Count[Compartlist, #] &, CompartIDs];
  sortedcomparts = CompartIDs[[Ordering[commons, All, Greater]]];
  colors = Max[3, Length[CompartIDs]];
  middle = Floor[colors/2]; colorgap = 1/(2. middle);
  colorseq = Riffle[Table[middle - i, {i, 0, colors/2}],Table[middle + i, {i, 1, colors/2}]]*colorgap;
  For[i = 1, i <= BlockCount, i++,
   hiro = CrossBlocks[[2, i, 2]]; loro = CrossBlocks[[2, i, 1]]; 
   If[hiro<loro,Continue[]];
   blockclist=Map[Compartlist[[#]] &, IntrowsperBlock[[i]]];
   maincolor = Position[sortedcomparts, BlockCompart[[i]]]; 
   maincolor = If[maincolor == {}, colorseq[[1]], colorseq[[maincolor[[1, 1]]]]];
   meancolor = Mean[Map[colorseq[[Position[sortedcomparts, #][[1, 1]]]] &, blockclist]];
   meancolor = Clip[meancolor, {maincolor - colorgap/3, maincolor + colorgap/3}];
   BlockHue[[i]]=meancolor;
   BBands[[loro + Offset ;; hiro + Offset, ;;]] = meancolor];
  BBands]
  
TickSpec[Chosen_,Ofset_]:=Module[{HT,VT,LefTicks,RiTicks,BotTicks,TopTicks,i},
	rowrange=Sinks+Sources-Ofset;
	LefTicks=Table[{i,i+Ofset},{i,1,rowrange,Max[5,Floor[rowrange/20,10]]}];
	BotTicks=Table[i,{i,1,Sinks,Max[Floor[Sinks/7,5],5]}];
	HT=Select[Chosen,#<= Sinks&];PrependTo[HT,1];
	TopTicks=AppendTo[HT,Sinks];
	VT=Select[Chosen,#>Sinks&]-Sinks;PrependTo[VT,1];AppendTo[VT,rowrange];
	RiTicks=Map[{#,#+Ofset}&,VT];
{{LefTicks,RiTicks},{BotTicks,TopTicks}}]

HorizontalTickSpec[Chosen_,Ofset_]:=Module[{HT,VT,LefTicks,RiTicks,BotTicks,TopTicks,i,rowrange},
	rowrange=Sinks+Sources-Ofset;
	LefTicks=Table[{i,i+Ofset},{i,1,rowrange,Max[5,Floor[rowrange/20,10]]}];
	BotTicks=Table[i,{i,1,Sinks,Max[Floor[Sinks/7,5],5]}];
	HT=Select[Chosen,#<= Sinks&];PrependTo[HT,1];
	TopTicks=AppendTo[HT,Sinks];
	VT=Select[Chosen,#>Sinks&]-Sinks;PrependTo[VT,1];AppendTo[VT,rowrange];
	RiTicks=Map[{#,#+Ofset}&,VT];
{{BotTicks,TopTicks},{LefTicks,RiTicks}}]

Eff[p_, FoundBlocks_] := 
 Module[{NodeCount, InternalCounts, Blockcount, OrphanCount, k}, 
  OrphanCount = Length[IsolationBlock]; 
  InternalCounts = 
   Join[Total[
     Map[Flatten@Differences[#, {0, 1}] + 1 &, CrossBlocks, 1]], 
    Length /@ FoundBlocks]; 
  InternalCounts = DeleteCases[InternalCounts, 0]; 
  (*  Remove dummy placeholder blocks*)
  Blockcount = Length[InternalCounts]; k = Blockcount + OrphanCount;
  NodeCount = Total[InternalCounts] + OrphanCount;
  Round[100*(
    Log[NodeCount^p] - 
     If[k > 1, Log[k^p + 1/k (Total[InternalCounts^p] + OrphanCount)],
       Log[NodeCount^p]])/(Log[NodeCount^p] - Log[2 NodeCount^(p/2)])]
  ]
  
End[] (* End Private Context *)

EndPackage[]