(* Mathematica Package *)

(* COPYRIGHT (C) 2010 Wynand Verwoerd 
under the terms of the GNU General Public License V3
This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
*)

BeginPackage["SelectExternals`",{"NetSplitter`","MatrixDisplay`","Printout`","RandomWalkBlocker`"}]
(* Exported symbols added here with SymbolName::usage *)  
ExternalSearch::usage="Finds candidate externals and administers user choices made in the selection dialog"

Begin["`Private`"] (* Begin Private Context *) 

SelectionDialog[Choice_,GreyLimit_]:=(Module[{
	Chosen=Choice,GreyChoice=GreyLimit,ChoiceOnly=True,Carryon,ChoosO,PrintBlocks,StyleBlock,Orphanstyles,Refusalstyles,
	Vetter,Ontheright,Targets,Ontheleft,Greytop,Greybot,i,hiro,loro,hicol,locol,blockclist,domina,domination},
Greybot=Ceiling[First[Greys],0.001];Greytop=Floor[Last[Greys],0.001];
If[Length[FoundBlocks]>0,PrintBlocks=Map[Partition[#,3,3,1," "]&,FoundBlocks];
						StyleBlock=Map[Style[#,If[MemberQ[TargetCompounds,#],Background->Targetinge," "]]&,PrintBlocks,{3}],
						StyleBlock=" "];
If[Length[IsolationBlock]>0,Orphanstyles=Map[Style[#,If[MemberQ[TargetCompounds,#],Background->Targetinge," "]]&,IsolationBlock],
							Orphanstyles={" "}];
If[Length[Refusals]>0,Refusalstyles=Map[Style[#,If[MemberQ[TargetCompounds,#],Background->Targetinge," "]]&,Refusals],
						Refusalstyles={" "}];

Vetter=CheckboxBar[Dynamic[Chosen],Table[Choice[[i]]->ToString[Choice[[i]]]~~"  "~~Nameseq[[Choice]][[i]],{i,1,Length[Choice]}],Appearance->"Vertical"];
Ontheright=Column[{ Panel[Grid[{{Tooltip[Style["Overlap compounds \n","Subsubsection"],"Select to approve as an external compound"]},
						{Vetter,SpanFromLeft},{"   "},{"Overlap Cutoff",InputField[Dynamic[GreyChoice,(GreyChoice=Max[Greybot,Min[Greytop,#]])&],FieldSize->5]}}]],
					Panel[Grid[{{Tooltip[Style["Selected Blocks \n","Subsubsection"],"Select by mouseclick, or type below in the format {2,4,5}"]},
						{InputField[Dynamic[SelectedBlocks,(SelectedBlocks=Select[#,0<#<= BlockCount&];ChoosO=(#!= {}))&],FieldSize->24],SpanFromLeft},
						{"Restrict next round to selected",Checkbox[Dynamic[ChoiceOnly],Enabled->Dynamic[ChoosO]]}}]],
					Panel[Grid[{
						{"Efficacy (" <> ToString[Eff[Effpower, FoundBlocks]] <> "%)",ProgressIndicator[Eff[Effpower, FoundBlocks]/100]},
						{"Next round ",DefaultButton["Continue",DialogReturn[True]]},
						{"Complete all rounds ",Button["Automatic",DialogReturn[2]]},
						{"Do reincorporation ",CancelButton["Finish",DialogReturn[False]]},
						{"Discard blocks ",Button["Retry",DialogReturn[3]]}
							}],Background->Lighter[Magenta,0.75]]
				},Spacings->2];
Targets=Sort@Flatten[Table[Position[Nameseq,TargetCompounds[[i]]],{i,1,Length[TargetCompounds]}]];

(* Identify the dominating compartment for each block *)
Colseq = Flatten[Table[Position[Complist, Sinkseq[[i]]], {i, 1, Length[Sinkseq]}]];
Rowseq = Flatten[Table[Position[Complist, Sourceseq[[i]]], {i, 1, Length[Sourceseq]}]];
RowCompartlist = Compartlist[[Rowseq]]; ColCompartlist = Compartlist[[Colseq]];
BlockCompart = ConstantArray["None", BlockCount]; IntrowsperBlock = ConstantArray[0, BlockCount];
BlockHue=ConstantArray[0,BlockCount];
For[i = 1, i <= BlockCount, i++,
	hiro = CrossBlocks[[2, i, 2]]; loro = CrossBlocks[[2, i, 1]];
	If[hiro < loro, Continue[]];
	hicol = CrossBlocks[[1, i, 2]]; locol = CrossBlocks[[1, i, 1]];
	IntrowsperBlock[[i]] = Join[Rowseq[[loro ;; hiro]], Colseq[[locol ;; hicol]]];
	blockclist = Join[RowCompartlist[[loro ;; hiro]], ColCompartlist[[locol ;; hicol]]];
	domina = Commonest[blockclist];
	domination = If[domina != {}, Count[blockclist, domina[[1]]]/Length[blockclist],0];
	BlockCompart[[i]] = 
		Which[domination > 0.6666, domina[[1]], domination == 0, "None", True, "Mixed"]];
 

Ontheleft=Module[{xbin,ybin,ClickBlock},
	xbin=Table[CrossBlocks[[1,i,1]]-1,{i,1,BlockCount}];AppendTo[xbin,CrossBlocks[[1,BlockCount,2]]];
	xbin = Concertina[xbin, 0.0001];
	ybin=Table[CrossBlocks[[2,i,1]]-1,{i,1,BlockCount}];AppendTo[ybin,CrossBlocks[[2,BlockCount,2]]];
	ybin = Concertina[ybin, 0.0001];
	Panel@ClickPane[Dynamic@Refresh[Show[{
		MatrixPlot[BlockBands[0], ColorFunction -> CompartHue,
			ColorFunctionScaling -> False, ImageSize -> 400, FrameTicks -> TickSpec[Chosen,Sinks]],
		MatrixPlot[PinkBlue[0, Targets, SelectedBlocks], ImageSize -> 400, 
   			ColorRules -> {0 -> None, 2 -> Targetinge, 4 -> Foundtinge, 6 -> Chosentinge}, 
   			PlotLabel -> Style["Blocking and Overlaps", "Subsection"]], 
  		MatrixPlot[ConsolBlocking, ColorRules -> {0. -> None}, ColorFunctionScaling -> False, ImageSize -> 400], 
  		MatrixPlot[Candidates[Chosen], ColorRules -> {0. -> None}, ColorFunction -> "SolarColors", ColorFunctionScaling -> False, ImageSize -> 400]},
    	Background -> White],TrackedSymbols->SelectedBlocks],
	(ClickBlock=Flatten@Position[Diagonal@BinCounts[{{#[[1]],Sources-#[[2]]}},{xbin},{ybin}],1];
	SelectedBlocks=BlockSelection[ClickBlock];ChoosO=(SelectedBlocks!= {}))&]];
CalculationTime+=TimeUsed[]-StartTime;
Carryon=DialogInput[TabView[{"Chooser"->Row[{Ontheleft,Ontheright},Spacer[10]],
							"Metabolites"->Dynamic[Column[{TextCell["Metabolite numbering for the blocking matrix ","Subtitle"],
								TextCell["Sinks","Subsection"],
								TextCell[Style[TableForm[Partition[CompoundStyles[SelectedBlocks][[1]],3,3,1,{}],TableSpacing->{1,5}],FontSize->10,Bold],"Output"],
								TextCell["Sources","Subsection"],
								TextCell[Style[TableForm[Partition[CompoundStyles[SelectedBlocks][[2]],3,3,1,{}],TableSpacing->{1,5}],FontSize->10,Bold],"Output"]},
								Background->LightYellow]],
							"Completed"->Dynamic[Column[{TextCell["Finalised Blocks","Subtitle"],TextCell["-------------------------"],
								TextCell[Style[TableForm[Map[Transpose,StyleBlock,{1}],TableSpacing->{4,5,0}],FontSize->10,Bold],"Output"]},
								Background->Foundtinge]],
							"Orphans"->Dynamic[Column[{TextCell["Isolated sinks","Subtitle"],TextCell["---------------------"],
								TextCell[Style[TableForm[Partition[Orphanstyles,3,3,1,{}],TableSpacing->{1,5}],FontSize->10,Bold],"Output"]},Background->LightPink]],
							"Refusals"->Dynamic[Column[{
								TextCell["Compounds previously refused status as externals","Subtitle"],
								TextCell["--------------------------------------------------------------------"],
								TextCell[Style[TableForm[Partition[Refusalstyles,3,3,1,{}],TableSpacing->{1,5}],FontSize->10,Bold],"Output"],
								TextCell["\nExternal specs that could not be matched","Subtitle"],
								TextCell["--------------------------------------------------------------------"],
								TextCell[Style[TableForm[Partition[NoMatch,3,3,1,{}]],FontSize->10,Bold],"Output"]},Background->LightYellow]],
							"Targets"->Dynamic[Column[{TextCell["Compounds of special interest","Subtitle"],
								TextCell["-----------------------------------------------"],
								TextCell[Style[TableForm[Partition[TargetCompounds,3,3,1,{}]],FontSize->10,Bold],"Output"]},Background->Targetinge]]}],
			WindowTitle->"Externals Selection: Round "~~ToString[Iterationcount],WindowElements->{"HorizontalScrollBar","VerticalScrollBar","MagnificationPopUp"},WindowSize->{All,Large},WindowFrameElements->{"ZoomBox"},WindowMargins->{{Automatic,Automatic},{10,0}}];
StartTime=TimeUsed[];If[SelectedBlocks== {},ChoiceOnly=False;];
{Chosen,GreyChoice,ChoiceOnly,Carryon}
])

Accept[i_,j_,Notrow_,Notcol_]:=Not[MemberQ[Notrow,i]&&MemberQ[Notcol,j]]

CompoundStyles[ChooseBlocks_]:=Module[{pickem,sinklist,sourcelist,Tinge},
	If[ChooseBlocks!= {},pickem=Transpose[MapAt[#+Sinks&,CrossBlocks,2]][[ChooseBlocks]];
		pickem=Map[Flatten,Transpose@Apply[Range,pickem,{2}]],pickem={{0},{0}}];
	sinklist=Table[{i,Nameseq[[i]]},{i,1,Sinks}];
	sourcelist=Table[{i,Nameseq[[i]]},{i,Sinks+1,Sinks+Sources}];
{Map[(Which[MemberQ[pickem[[1]],#[[1]]],Tinge=Chosentinge,MemberQ[TargetCompounds,#[[2]]],Tinge=Targetinge,True,Tinge=None];
	Style[ToString[#[[1]]]<>"  "<>ToString[#[[2]]],Background->Tinge])&,sinklist],
Map[(Which[MemberQ[pickem[[2]],#[[1]]],Tinge=Chosentinge,MemberQ[TargetCompounds,#[[2]]],Tinge=Targetinge,True,Tinge=None];
	Style[ToString[#[[1]]]<>"  "<>ToString[#[[2]]],Background->Tinge])&,sourcelist]}]

ExternalSearch[]:=Module[{Carryon,RowRefuse,ColRefuse,Choice,GreyChoice,LightGreys,
	ConstraintCount,ConstraintPairs,RightSides,Refused,Constraints,
	Selected,SelectOnly,NewRefused={},FinishedBlocks,IRrows,LeftoverB,LeftoverBlocks},
Off[LinearProgramming::"lpipp"];Off[LinearProgramming::"lpip"];Off[LinearProgramming::"lplop"];
Iterationcount=0;IsolationBlock={};FoundBlocks={};BlockSinks={};
IRB={};ChosenExternals={};
DecisionTrace={};MergeTrace={};
GreyChoice=GreyLimit=0.1;
StartEff = 0; PreMergeEff = 0; PostMergeEff = 0;UserChoose=True;

If[Length[Externalrows]==InternalCNodes,EmitSound[Crisis];DialogInput[{"Improper network - all metabolites are external!",DefaultButton[DialogReturn[]]}];Return[]];
If[RepeatMerge, Carryon = False; BlockMaker[Externalrows]; 
  StartEff = Eff[Effpower, {}];, Carryon = True];
While[Carryon, ++Iterationcount; Status = "Round "<>ToString[Iterationcount]<>": Calculating a new DAG";
 If[BlockMaker[Externalrows] == False, Carryon = False; Return[]];
 If[Iterationcount == 1, StartEff = Eff[Effpower, {}]];
 	Status="Round "<>ToString[Iterationcount]<>": Find overlaps and propose externals";
	RowRefuse=Flatten[Table[Position[Sourceseq,Refusals[[i]]],{i,1,Length[Refusals]}]];
	ColRefuse=Flatten[Table[Position[Sinkseq,Refusals[[i]]],{i,1,Length[Refusals]}]];
	Greys=Reap[Do[If[Accept[i,j,RowRefuse,ColRefuse]&&0<ConsolBlocking[[i,j]]<GreyMax,
		Sow[ConsolBlocking[[i,j]]]],{i,1,Sources},{j,1,Sinks}]];
	GreyLimit=0;
	If[Depth[Greys]>3,Greys=Sort@Greys[[2,1]];Howmany=Min[Length[Greys],Manygreys*Sinks];
		GreyChoice=Ceiling[Greys[[Howmany]],0.001],
		GreyChoice=0; Choice=Chosen={};Carryon=False;
		If[UserChoose,ChoiceDialog["Splitting finished - no more recognised overlaps.
					\nInspect the result by clicking Merger or Printed Output.",
					{" OK "->1},WindowFloating->True,ButtonBoxOptions->{ImageMargins->6}]]];
		While[GreyChoice!=GreyLimit,GreyLimit=GreyChoice;
			LightGreys=Reap[Do[If[Accept[i,j,RowRefuse,ColRefuse]&&0<ConsolBlocking[[i,j]]&&ConsolBlocking[[i,j]]<= GreyLimit,
				Sow[{i,j}]],{i,1,Sources},{j,1,Sinks}]];
			If[Depth[LightGreys]== 5,LightGreys=LightGreys[[2,1]];ConstraintCount=Length[LightGreys],LightGreys={};ConstraintCount=0];
			If[ConstraintCount>0,ConstraintPairs=Flatten[MapIndexed[{{First[#2],#1[[1]]+Sinks}->1,{First[#2],#1[[2]]}->1}&,LightGreys]];
				RightSides=Table[{1,1},{i,ConstraintCount}];
				Refused=Flatten[Table[Position[Nameseq,Refusals[[i]]],{i,1,Length[Refusals]}]];
				ConstraintPairs=Join[ConstraintPairs,Table[{ConstraintCount+i,Refused[[i]]}->1,{i,1,Length[Refused]}]];
				RightSides=Join[RightSides,Table[{0,0},{i,1,Length[Refused]}]];
				ConstraintCount+=Length[Refused];
				Constraints=SparseArray[ConstraintPairs,{ConstraintCount,Sources+Sinks}];
				Selected=LinearProgramming[Table[1,{i,Sources+Sinks}],Constraints,RightSides,0,Integers];
				Choice=Flatten[Position[Selected,1]];
				SelectOnly=False;SelectedBlocks={};

				If[UserChoose,EmitSound[Attention];{Chosen,GreyChoice,BlockChoice,Carryon}=SelectionDialog[Choice,GreyLimit],
						Chosen=Choice;BlockChoice=False];
				Switch[Carryon,2,UserChoose=False;Carryon=True,3,Carryon=False;Return[]];,

				Choice=Chosen={};Carryon=False;EmitSound[Attention];
				ChoiceDialog["Splitting terminated as no new externals can be found at the present overlap cutoff.
					 		\nInspect the result by clicking Merger or Printed Output.",
					 {" OK "->1},WindowFloating->True,ButtonBoxOptions->{ImageMargins->6}]
			];

(* If[!Carryon,Abort[]];  DEBUGGING *)

		];
		If[Carryon,NewExternals=Join[NewExternals,Nameseq[[Chosen]]];
			NewRefused=Complement[Choice,Chosen];
			AppendTo[DecisionTrace,{Iterationcount,If[BlockChoice,SelectedBlocks,"All"],If[NewRefused== {},"None",NewRefused]}];
			Refusals=Join[Refusals,Nameseq[[NewRefused]]];
			Externalrows=Flatten[Table[Position[InternalComplist,NewExternals[[i]]],{i,1,Length[NewExternals]}]];
			IRB=Select[Transpose[MapAt[#+Sinks&,CrossBlocks,2]],#[[1,1]]== #[[1,2]]||#[[2,1]]== #[[2,2]]&];
			If[BlockChoice,FinishedBlocks=Transpose[MapAt[#+Sinks&,CrossBlocks,2]][[Complement[Range[BlockCount],SelectedBlocks]]];IRB=Union[IRB,FinishedBlocks]];
			If[Length[IRB]== Dimensions[CrossBlocks][[2]],Carryon=False;EmitSound[Attention];
				ChoiceDialog["Finished, as there are no further blocks to process.\nInspect the result by clicking Merger or Printed Output.",{" OK "->1},WindowFloating->True,ButtonBoxOptions->{ImageMargins->6}]];
			IRBlocks=Apply[Range,IRB,{2}];IRBlocks=Map[Flatten,IRBlocks];
			IRBlocks=Table[Nameseq[[IRBlocks[[i]]]],{i,1,Length[IRBlocks]}];
			FoundBlocks=Join[FoundBlocks,IRBlocks];
			BlockSinks=Join[BlockSinks,Table[IRB[[i,1,2]]-IRB[[i,1,1]]+1,{i,1,Length[IRB]}]];
		(*If[IRBlocks!= {},Print[Row[{" A total of ",Length[IRBlocks]," irreducible blocks extracted in this round."}]],];*)
			IRrows=Flatten[Table[Position[InternalComplist,Flatten[IRBlocks][[i]]],{i,1,Length[Flatten[IRBlocks]]}]];
			Externalrows=Union[Externalrows,IRrows];
			Removedrows=Chosen;Externals=Nameseq[[Removedrows]];
			ChosenExternals=Join[ChosenExternals,Externals];
		];
	];

LeftoverB=Complement[Transpose[MapAt[#+Sinks&,CrossBlocks,2]],IRB];
LeftoverBlocks=Apply[Range,LeftoverB,{2}];LeftoverBlocks=Map[Flatten,LeftoverBlocks];
LeftoverBlocks=Table[Nameseq[[LeftoverBlocks[[i]]]],{i,1,Length[LeftoverBlocks]}];
BlockSinks=Join[BlockSinks,Table[LeftoverB[[i,1,2]]-LeftoverB[[i,1,1]]+1,{i,1,Length[LeftoverB]}]];
FinalBlocks=Join[FoundBlocks,LeftoverBlocks,{IsolationBlock}];FoundBlocks = {};
FinalBlockCount=Length[FinalBlocks] ;
AppendTo[BlockSinks,Length[IsolationBlock]];
Sinks=Total[BlockSinks];Sources=InternalCNodes-Sinks;
Sinkseq=Join[Sinkseq,IsolationBlock];Nameseq=Join[Sinkseq,Sourceseq];

FillOut[True];PreMergeExternals=Externals;PreMergeEff = Eff[Effpower, {}];
(* Perform any merge sequence loaded from file *)
If[RepeatMerge,Map[Combine,OldMerge]];
EmitSound[Attention];SplitDone=True;
]


End[] (* End Private Context *)

EndPackage[]