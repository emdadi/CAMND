(* Mathematica Package *)

(* COPYRIGHT (C) 2010 Wynand Verwoerd 
under the terms of the GNU General Public License V3
This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
*)

BeginPackage["Printout`",{"NetSplitter`","FileReadWrite`","MatrixDisplay`"}]
(* Exported symbols added here with SymbolName::usage *)  
Printout::usage="Produces a printout on screen and sends this to printer or file"
PrintBlocks::usage="Formats block metabolite listings in 3 columns"

Begin["`Private`"] (* Begin Private Context *) 

PrintBlocks[ChooseBlocks_]:=Module[{pickem,picked,compnos,NumberBlocks,PB,Tinge},
	PB=If[Length[FinalBlocks]>0,compnos=Map[Position[Nameseq,#][[1,1]]&,FinalBlocks,{2}];NumberBlocks=compnos[[1;;BlockCount]];
		If[ChooseBlocks!= {},pickem=Transpose[MapAt[#+Sinks&,CrossBlocks,2]][[ChooseBlocks]];pickem=Flatten@Transpose@Apply[Range,pickem,{2}],
			pickem={0}];picked=Nameseq[[pickem]];
			Do[Which[MemberQ[TargetCompounds,FinalBlocks[[i,j]]],Tinge=Targetinge,MemberQ[picked,FinalBlocks[[i,j]]],Tinge=Chosentinge,True,Tinge=None];
				NumberBlocks[[i,j]]=Style[ToString[compnos[[i,j]]]<>"  "<>FinalBlocks[[i,j]],Background->Tinge],{i,1,BlockCount},{j,1,Length[compnos[[i]]]}];
				Map[Partition[#,3,3,1," "]&,NumberBlocks]," "];
				PB=Select[PB,Length[#]>0&];Map[Transpose,PB,{1}]]

ShowSubnet[Blockno_]:=Module[{BlockIntRows,BlockIntCols,BlockRows,BlockExtRows,Inflows,Outflows,Crossflows,lines,Bname,Detailstyles},
	If[Length[FinalBlocks]>0 && Length[FinalBlocks[[Blockno]]]>0 ,
		BlockIntRows=Flatten[Table[Position[Complist,FinalBlocks[[Blockno,i]]],{i,1,Length[FinalBlocks[[Blockno]]]}]];
		BlockIntCols=Nonzerows[Normal[Transpose[S[[BlockIntRows]]]]];
		BlockRows=Nonzerows[Normal[S[[All,BlockIntCols]]]];
		BlockExtRows=Complement[BlockRows,BlockIntRows];
		Inflows=BlockExtRows\[Intersection]InputExternalrows;
		Outflows=BlockExtRows\[Intersection]OutputExternalrows;
		Crossflows=Complement[BlockExtRows,Inflows\[Union]Outflows];
		lines =Max[{Length[Inflows],Length[Crossflows],Length[Outflows]}];
		Inflows=PadRight[Complist[[Inflows]],lines," "];
		Crossflows=PadRight[Complist[[Crossflows]],lines," "];
		Outflows=PadRight[Complist[[Outflows]],lines," "];
		If[lines>0,Detailstyles=Map[Style[#,If[MemberQ[TargetCompounds,#],Background->Targetinge," "]]&,{ConstantArray[" ",lines],Inflows,Crossflows,Outflows},{2}],Detailstyles={" "}];
		Bname=If[Blockno<FinalBlockCount,"Block "<>ToString[Blockno],"Orphan Block "];
		TextCell[Style[TableForm[Detailstyles,TableDirections->Row,TableSpacing->{1,0},TableHeadings->{{Bname,"Inflows","Crossflows","Outflows"}}],FontFamily->"Arial Narrow",FontSize->7,Bold],"Output"],
		Null]
]

Printout[]:=(Module[{Inpath,Modeldir,compnos,NumberBlocks,NontrivialBlockCount,PrintBlocks,BlockHeaders,Tinge,i,j,Orphanstyles,OrphanRows,Harvest,Blocklength,BlockIntRows,BlockIntCols,BlockExtRows,BlockRows,Overlaprows,Overlaplist,Otherblocks,OverlapReport,NonStructx,picked,Shortlist,Survivors,OrphanToExt,Targets,Drop2,Drop3,nb,dest},
	Inpath=DirectoryName[Infile];Modeldir=StringTake[Inpath,Transpose[Take[StringPosition[Inpath,"\\"],-2]][[1]]+{1,-1}];
	Status="Preparing the printout";
	If[Length[FinalBlocks]>0,
		NontrivialBlockCount=Count[FinalBlocks,x_/;Length[x]>0]-1;
		compnos=Map[Position[Nameseq,#][[1,1]]&,FinalBlocks,{2}];
		NumberBlocks=compnos[[1;;BlockCount]];
		Do[If[MemberQ[TargetCompounds,FinalBlocks[[i,j]]],Tinge=Targetinge,Tinge=None];NumberBlocks[[i,j]]=Style[ToString[compnos[[i,j]]]<>"  "<>FinalBlocks[[i,j]],Background->Tinge],{i,1,BlockCount},{j,1,Length[compnos[[i]]]}];
		PrintBlocks=Map[Partition[#,PrintColumns,PrintColumns,1," "]&,NumberBlocks];,PrintBlocks=" "];PrintBlocks=Select[PrintBlocks,Length[#]>0&];
		BlockHeaders=Table["Block "<>ToString[i]<>"\n"<>BlockCompart[[i]],{i,1,FinalBlockCount}];
		BlockHeaders=Pick[BlockHeaders,Map[Length[#]>0&,FinalBlocks]];
		BlockHeaders=If[IsolationBlock=={},{BlockHeaders},{Most[BlockHeaders]}];
		NonStructx=Complement[Externalrows,StructuralExternalrows];
		picked=Extract[ConnectionExternalrows,Position[ConnectionExternalrows,x_/;MemberQ[NonStructx,x],{1},Heads->False]];
		Shortlist=Transpose[{Connex[[picked]],Complist[[picked]]}];
		Survivors=ChosenExternals\[Intersection]Externals;
		OrphanToExt=Complement[Complist[[Complement[NonStructx,ConnectionExternalrows]]],StoichioExternals,Extras,Survivors];
		If[Length[IsolationBlock]>0,compnos=Map[Position[Nameseq,#][[1,1]]&,IsolationBlock,{1}];Orphanstyles=MapThread[(If[MemberQ[TargetCompounds,#2],Tinge=Targetinge,Tinge=None];Style[ToString[#1]<>"  "<>#2,Background->Tinge])&,{compnos,IsolationBlock}],compnos={};Orphanstyles={}];
		OrphanRows=Flatten[Table[Position[Complist,IsolationBlock[[i]]],{i,1,Length[IsolationBlock]}]];
		Harvest=Reap[Do[Blocklength=Length[FinalBlocks[[Blockno]]];
			If[Blocklength>0,BlockIntRows=Flatten[Table[Position[Complist,FinalBlocks[[Blockno,i]]],{i,1,Blocklength}]];
				BlockIntCols=Nonzerows[Normal[Transpose[S[[BlockIntRows]]]]];
				BlockRows=Nonzerows[Normal[S[[All,BlockIntCols]]]];
				BlockExtRows=Complement[BlockRows,BlockIntRows];
				Overlaprows=BlockExtRows\[Intersection]OrphanRows;
				Overlaplist=Complist[[Overlaprows]];
				Otherblocks={};Do[Otherblocks=Join[Otherblocks,If[BlockExtRows\[Intersection]IntrowsperBlock[[i]]!= {},{i},{}]],{i,Blockno+1,BlockCount}];
				Which[Otherblocks!= {}&&Overlaprows!= {},Sow[{Blockno,Otherblocks,Overlaplist}],Otherblocks!= {},Sow[{Blockno,Otherblocks,"No overlaps"}],Overlaprows!= {},Sow[{Blockno,"No overlaps",Overlaplist}]];];
				,{Blockno,1,BlockCount}]];
				OverlapReport=If[Depth[Harvest]== 6,Harvest[[2,1]],"No block overlaps found"];
				Targets=Sort@Flatten[Table[Position[Nameseq,TargetCompounds[[i]]],{i,1,Length[TargetCompounds]}]];
				Drop2=If[DropRXN!= {},Style[TableForm[Partition[RemovedRXN,PrintColumns,PrintColumns,1,{}],TableSpacing->{0,8},TableHeadings->{{"Reactions "}}],FontFamily->"Arial Narrow",FontSize->7,Bold],"No reactions removed"];
				Drop3=If[DropComps!= {},Style[TableForm[Partition[RemovedComp,PrintColumns,PrintColumns,1,{}],TableSpacing->{0,8},TableHeadings->{{"Metabolites "}}],FontFamily->"Arial Narrow",FontSize->7,Bold],"No metabolitess removed"];
				nb=CreateDocument[{TextCell[Row[{Style["Netsplitter results: "<>Modeldir,"Subtitle"],
					Style["\nStoichiometry input file: "<> Infile],
					Style["\n\nNotes on the model: "<>ModelNote],
					Style["\nA total of "<>ToString[RNodes]<>" reactions were imported and "<>ToString[Reversecount]<>" of these specified as reversible."],
					Style["\n A total of "<>ToString[CNodes]<>" metabolites were listed, of which "<>ToString[Length[StoichioExternals]]<>" were specified as external in the stoichiometry file."],
					Style[Externalmess],Style[Fluxmess],Style[Targetmess],
					Style["\nCalculated at "<>DateString[]],
					Style["\nA total of "<>ToString[Iterationcount]<>" selection rounds were performed "<>If[BlockCount>1," resolving "<>ToString[NontrivialBlockCount]<> " non-trivial subnet blocks and "<>ToString[Length[IsolationBlock]]<>" orphans.","but no block structure could be resolved."]],
					Style["\n\tSplit Efficacy = " <> ToString[StartEff] <> "% (input externals only); " <> ToString[PreMergeEff] <> "% (after selection and reincorporation); and " <> 
						If[Stepcount > 0, ToString[PostMergeEff] <> "% (after merging).", "no merging was done."]],
					"\nRunning time for this calculation = "<>ToString[CalculationTime]<>" seconds."}]],
					TextCell["Optional parameters used","Subsection"],
					TextCell[Row[{"    Reactions omitted with ID containing:\t"<>If[Dropstring== "","No omissions",Dropstring],"\n    Max connectivity kept as internal:\t"<>ToString[ConnexMax],"\n    Binary vector similarity measure:\t"<>ToString[DistChoice],"\n    Intercluster distance:           \t"<>ToString[LinkChoice],"\n    Interblock overlap cutoff value =\t"<>ToString[NumberForm[GreyLimit,4]],"\n    Grouplevel search interval      =\t"<>ToString[GroupLevelRange],"\n    Optimised Grouplevel value      =\t"<>ToString[NumberForm[GroupLevel,4]],If[ReverseThem,"\n    Reversible reactions were expanded by duplication in the reverse direction.","\n    Reversibility ignored - only explicit reaction directions included. "]}]],
					Grid[{{Show[{
						MatrixPlot[BlockBands[Sinks], ColorFunction -> CompartHue,ColorFunctionScaling -> False, ImageSize -> 400, FrameTicks -> None],
						MatrixPlot[PinkBlue[Sinks,Targets,SelectedBlocks],ImageSize->400,ColorRules->{0->None,2->Targetinge,4->Foundtinge,6->Chosentinge},PlotLabel->Style["Rearranged DAG matrix","Section"]],
						MatrixPlot[Map[If[#== 0,None,#]&,ShuffleDAG,{2}],ImageSize->400,ColorRules->{0->None}]}]}},
						Spacings->Scaled[0.05]] ,
					TextCell["\n"<>ToString[InternalCNodes]<>" Internal Metabolites","Subsubtitle"],TextCell["Metabolites in each block","Subsubsection"],TextCell["Numbering as in the plot; target metabolites are highlighted.","Text"],
					TextCell[Style[TableForm[Map[Transpose,PrintBlocks,{1}],TableSpacing->{4,2,0},TableHeadings->BlockHeaders],FontFamily->"Arial Narrow",FontSize->7],"Output",PageBreakWithin->True],
					TextCell[ToString@Length[IsolationBlock]<>" Orphans (subnets with a single internal metabolite)","Subsubsection"],
					TextCell[Style[TableForm[Partition[Orphanstyles,PrintColumns,PrintColumns,1,{}],TableSpacing->{0,2}],FontFamily->"Arial Narrow",FontSize->7,Bold],"Output"],
					TextCell["Externals/internals overlap between blocks","Subsubsection"],
					TextCell[Style[If[Length[OverlapReport]== 0,OverlapReport,Grid[Prepend[OverlapReport,{"Block","overlaps with blocks","and with isolated sinks"}],Dividers->{Center,{False,True}}]],FontFamily->"Arial",FontSize->7],"Output"],
					TextCell[ToString[CNodes-InternalCNodes]<>" External Metabolites","Subsubtitle"],
					TextCell[ToString@Length[Refusals]<>" Metabolites refused status as externals ","Subsubsection"],
					TextCell[Style[TableForm[Partition[Refusals,PrintColumns,PrintColumns,1,{}],TableSpacing->{0,2}],FontFamily->"Arial Narrow",FontSize->7],"Output"],
					TextCell[ToString@Length[StructuralExternalrows]<>" Structural externals","Subsubsection"],
					TextCell["These are metabolites that only ever appear as either an inflow or an outflow ","Output",FontSize->8],
					TextCell[ToString@Length[Shortlist]<>" High connectivity metabolites (with connectivity counts)","Subsubsection"],
					TextCell[Style[TableForm[Partition[Map[ToString[#[[1]]]<>": "<>#[[2]]&,Shortlist],PrintColumns,PrintColumns,1,{}],TableSpacing->{0,2}],FontFamily->"Arial Narrow",FontSize->7,Bold],"Output"],
					TextCell[ToString@Length[Complement[StoichioExternals,BasicExternals]]<>" Metabolites specified in the stoichiometry specification","Subsubsection"],
					TextCell[Style[TableForm[Partition[Complement[StoichioExternals,BasicExternals],PrintColumns,PrintColumns,1,{}],TableSpacing->{0,2}],FontFamily->"Arial Narrow",FontSize->7,Bold],"Output"],
					TextCell[ToString@Length[Extras\[Intersection]Externals]<>" Metabolites proposed in the externals file ","Subsubsection"],
					TextCell[Style[TableForm[Partition[Extras\[Intersection]Externals,PrintColumns,PrintColumns,1,{}],TableSpacing->{0,2}],FontFamily->"Arial Narrow",FontSize->7,Bold],"Output"],
					TextCell[ToString@Length[Survivors]<>" Metabolites chosen interactively ","Subsubsection"],
					TextCell[Style[TableForm[Partition[Survivors,PrintColumns,PrintColumns,1,{}],TableSpacing->{0,2}],FontFamily->"Arial Narrow",FontSize->7,Bold],"Output"],
					TextCell[ToString@Length[OrphanToExt]<>" Ex-orphans becoming external due to reincorporation ","Subsubsection"],
					TextCell[Style[TableForm[Partition[OrphanToExt,PrintColumns,PrintColumns,1,{}],TableSpacing->{0,2}],FontFamily->"Arial Narrow",FontSize->7,Bold],"Output"],
					TextCell["Blockwise details for external metabolites","Subsubtitle"],
					Sequence@@Cases[Table[ShowSubnet[i], {i, 1, FinalBlockCount}], Except[Null]],
					TextCell["Items deleted from Stoichiometry","Subsubsection"],
					TextCell[Drop2,"Output"],TextCell[Drop3,"Output"],
					TextCell["Record of interactive choices","Subsubsection"],
					TextCell[Style[Grid[Prepend[DecisionTrace,{"Selection Round","Restrict to blocks","Refused externals"}],Dividers->{Center,{False,True}}],FontFamily->"Arial",FontSize->10],"Output"],
					TextCell[Style[Grid[Prepend[MergeTrace,{"Merge step","Merged Blocks","Included Orphans","New Block"}],Dividers->{Center,{False,True}}],FontFamily->"Arial",FontSize->10],"Output"]
					},WindowTitle->"Subnetwork Output",WindowMargins->{{Automatic,Automatic},{0,Automatic}},ShowCellBracket->False];
				SetOptions[nb,PrintingOptions->{"PrintingMargins"->{{36,36},{36,36}}}];
				SelectionMove[nb,Next,Cell,9,AutoScroll->False];FrontEndExecute[FrontEndToken[nb,"SelectionConvert","PostScript"]];(*SetOptions[nb,PrintingOptions->{"PrintMultipleHorizontalPages"->True}];*)
				SetOptions[nb,ShowPageBreaks->True];
				dest=ChoiceDialog[" ",{"Printer"->"printer","PDF File"->ToFileName[Inpath,"NSPL_Printout.pdf"],"Postscript file"->ToFileName[Inpath,"NSPL_Printout.ps"],"Mathematica Notebook"->"notebook","Keep on Screen"->"open","Cancel"->"exit"},WindowTitle->"Where do you want this ouput to go?",WindowFloating->True,WindowSize->{72 7.0,All},WindowMargins->{{0,Automatic},{Automatic,0}},ButtonBoxOptions->{ImageMargins->4}];
				If[dest=="printer",FrontEndExecute[FrontEndToken["SystemPrintOptionsDialog"]]];
				Which[dest== "exit",Null,dest=="open",Null,dest== "notebook",NotebookSave[nb,ToFileName[Inpath,"NSPL_Printout.nb"]],dest== "printer",FrontEndExecute[FrontEndToken[nb,"PrintDialog"]],True,NotebookPrint[nb,dest]];
				If[dest!= "open",NotebookClose[nb]];
])


End[] (* End Private Context *)

EndPackage[]