(* Mathematica Package *)

(* COPYRIGHT (C) 2010 Wynand Verwoerd 
under the terms of the GNU General Public License V3
This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
*)

BeginPackage["SubnetLayout`",{"NetSplitter`","FileReadWrite`"}]
(* Exported symbols added here with SymbolName::usage *)  
DrawFullnet::usage="Draws the full network before any splitting"
DrawSubnet::usage="Draws subnets with dynamic adjustment of fontsize"
SafeDrawSubnet::usage="Draws subnets, less dynamic adjustment avoids Math 6 crash"
GraphSpec::usage="  "

Begin["`Private`"] (* Begin Private Context *) 

DrawFullnet[]:=DynamicModule[{Inpath,Paper, dest,title,gspec,BlockRows,BlockStoi,nb,Asize,Orient},
	Status="Preparing layout for full network";
	Inpath=DirectoryName[Infile];Asize="A3";Orient="Landscape";
	title="Full network, excluding high connectivity metabolites";
	DontShow=Extract[Externalrows,Position[Externalrows,x_/;Connex[[x]]>ConnexMax,{1},Heads->False]];
	BlockRows=Complement[Range[Length[S]],DontShow];
	BlockStoi=S[[BlockRows,All]];BlockStoi=Normal@BlockStoi/.{x_/;x>0->1,x_/;x<0->2};
	gspec=Reap[Map[(If[#[[2]]>1,Sow[#[[1,1]]->#[[1,2]]]];
		If[OddQ[#[[2]]],Sow[#[[1,2]]->#[[1,1]]]];)&,Drop[ArrayRules[BlockStoi],-1]]][[2,1]];
	nb=CreateDocument[{TextCell[title],
			Framed[LayeredGraphPlot[gspec,MultiedgeStyle->False,SelfLoopStyle->None,
				PackingMethod->"LayeredLeft"(*"ClosestPacking"*),PlotStyle->Arrowheads[{{Small,0.5}}],
				ImageSize->Dynamic[Imsize[Asize,Orient],SynchronousUpdating->False],
				PlotRangePadding->None,
				AspectRatio->Dynamic[Aspect[Orient],SynchronousUpdating->False]],
			FrameMargins->2]},
		WindowMargins->{{40,40},{40,40}},WindowTitle->"Network Layout",
		WindowFrameElements->{"ZoomBox","ResizeArea"},ShowCellBracket->False];
	Paper=Grid[{{"Paper Size",RadioButtonBar[Dynamic[Asize],{"A4","A3","A2","A1","A0"}]},
				{"Orientation",RadioButtonBar[Dynamic[Orient],{"Tall","Portrait","Landscape","Wide"}]}},Background->LightYellow,Frame->All];
	dest=ChoiceDialog[Row[{"    ",Paper}],
				{"Printer"->"printer",
				"   PDF   \n   File"->ToFileName[Inpath,"NSPL_Fullnet.pdf"],
				"Cancel"->"exit"},
			WindowFloating->True,WindowMargins->{{0,Automatic},{Automatic,0}},
			WindowTitle->"Choose the layout size",ButtonBoxOptions->{ImageMargins->4}];
	NotebookClose[nb];
	If[dest=="exit",Return[]];
	nb=CreateDocument[{TextCell[title],
			Style[Framed[LayeredGraphPlot[gspec,MultiedgeStyle->False,SelfLoopStyle->None,
				PackingMethod->"LayeredLeft"(*"ClosestPacking"*),PlotStyle->Arrowheads[{{Small,0.5}}],
				ImageSize->Imsize[Asize,Orient],
				PlotRangePadding->None,AspectRatio->Aspect[Orient]],
				FrameMargins->2],Magnification->1]},
			ShowCellBracket->False];
	SetOptions[nb,PrintingOptions->{"PageSize"->PagePointsize[Asize,Orient],
			"PrintingMargins"->{{18,18},{0,18}},
			"PrintMultipleHorizontalPages"->True},
			PageFooters->{{None,None,None},{None,None,None}},
			PageHeaders->{{None,None,None},{None,None,None}}];
	If[dest== "printer",SetOptions[nb,PrintingOptions->{"PrintMultipleHorizontalPages"->PrintMultHorPages}]];
	If[dest!="exit",FrontEndExecute[FrontEndToken["SystemPrintOptionsDialog"]]];
	Which[dest== "exit",Null,dest== "printer",FrontEndExecute[FrontEndToken[nb,"PrintDialog"]],True,NotebookPrint[nb,dest]];
	NotebookClose[nb];
]

DrawSubnet[blockno_,blockname_,minsize_,showmax_]:=DynamicModule[
	{Inpath,Paper, dest,margins,paperorient,NodesnEdges,fontsize=7,Wider=1.0,Taller=1.0,graphtype=LayeredGraphPlot,
	CharSize,BoxWide=0.1,BoxHigh=0.1,BoxGrow,title,nb,filename,OutOptions,Disposal,whereto,Asize,Orient},
	Status="Preparing a network layout";
	Inpath=DirectoryName[Infile];
	filename=Which[blockno==0,"NSPL_Metanet",blockno==FinalBlockCount,"NSPL_Orphan net",True,"NSPL_Subnet "<>ToString[blockno]];
	Asize="A4";Orient="Portrait";
	title=blockname<>" (Max connect = "<>ToString[showmax]<>")";
	NodesnEdges=GraphSpec[blockno,minsize,fontsize,showmax];
	CharSize=Dynamic[fontsize/Imsize[Asize,Orient]]; 
	BoxGrow=Dynamic[If[MatchQ[graphtype,LayeredGraphPlot],3,1.8]];
	BoxWide=Dynamic[5*BoxGrow*Wider*CharSize]; BoxHigh=Dynamic[4*BoxGrow*Taller*CharSize];
	nb=CreateDocument[{TextCell[If[blockno== 0,"Metanetwork :"<>title,
		"Network Diagram for block "<>ToString[blockno]<>": "<>title]],
		Framed[Dynamic[Refresh[If[MatchQ[graphtype, GraphPlot3D],
 			GraphPlot3D[NodesnEdges, MultiedgeStyle -> False, VertexLabeling -> True,
 				 SelfLoopStyle->None, Boxed->False, PlotStyle->Arrowheads[{{Small, 0.5}}],
 				 ImageSize->Imsize[Asize,Orient],PlotRangePadding->None,AspectRatio->Aspect[Orient]],
 			graphtype[NodesnEdges,MultiedgeStyle->False,VertexLabeling->True,SelfLoopStyle->None,DirectedEdges->True,
			PackingMethod->Switch[Orient,"Tall","LayeredTop","Portrait","ClosestPacking","Landscape","LayeredLeft","Wide","LayeredLeft"],
			VertexRenderingFunction -> (If[StringMatchQ[#2[[1]], "SUBNET"~~__],
				{Hue[BlockHue[[GetBlock[#2]]], 0.15, 1], EdgeForm[Orange], Disk[#1, {4*BoxWide, BoxHigh*(1.475/Aspect[Orient])}], Text[#2, #1]} ,
				{LightYellow, EdgeForm[Orange],
					If[MatchQ[#2[[5]], Black|Purple],
						Disk[#1, {4*BoxWide, BoxHigh*(1.475/Aspect[Orient])}],
						Rectangle[#1-{4*BoxWide, BoxHigh*(1.475/Aspect[Orient])}, #1+{4*BoxWide, BoxHigh*(1.475/Aspect[Orient])}]],
				 Text[#2, #1]}
				]&),
			PlotStyle->Arrowheads[{{Small,0.5}}],ImageSize->Imsize[Asize,Orient],
			PlotRangePadding->None,AspectRatio->Aspect[Orient]]],
			TrackedSymbols->{Asize,Orient,fontsize,Wider,Taller,graphtype}]],
		FrameMargins->2]},
	WindowTitle->"Network Layout",WindowFrameElements->{"ZoomBox","ResizeArea","CloseBox"},
	ShowCellBracket->False];
	OutOptions = {
		{"exit" -> "None", "open"->"Window",ToFileName[Inpath, filename<>"3D.nb"]->"notebook" }, 
		{"exit" -> "None", "open"->"Window", ToFileName[Inpath, filename<>".pdf"]->"PDF File", 
    	 ToFileName[Inpath, filename<>".ps"]->"PS File", "printer" -> "Printer"}};
    Disposal=OutOptions[[2]];
	Paper=Grid[{{"Paper Size",RadioButtonBar[Dynamic[Asize],
		{"A4","A3","A2","A1","A0"}]},
		{"Orientation",RadioButtonBar[Dynamic[Orient],
		{"Tall","Portrait","Landscape","Wide"}]},
		{"Font Size",RadioButtonBar[Dynamic[fontsize,
			(NodesnEdges=Replace[NodesnEdges,Rule[FontSize,fontsize]->Rule[FontSize,#],Infinity];fontsize=#)&],
			{1,4,5,6,7,8,9,10}]},
		{"Node width", Slider[Dynamic[Wider],{0.1,5}]},
		{"Node height", Slider[Dynamic[Taller],{0.1,5}]},
		{"Layout type", RadioButtonBar[Dynamic[graphtype, 
			(Disposal=OutOptions[[If[MatchQ[#, GraphPlot3D], 1, 2]]]; graphtype=#)&],
			{LayeredGraphPlot->"Layered", GraphPlot->"Centred", GraphPlot3D->"3D Interactive"}]},
		{"Output to", Dynamic[SetterBar[Dynamic[whereto], Disposal, Appearance -> "Row", 
			Background -> LightCyan, ImageMargins -> 3]]}
			},
		Background->LightYellow,Frame->All];
	dest = ChoiceDialog[Row[{"                 ", Paper}],
			WindowFloating->True, WindowMargins->{{0,Automatic},{Automatic,0}}, 
			(*WindowSize -> {500, 250},*) WindowTitle -> "Layout Configuration"];
	If[whereto!= "open",NotebookClose[nb]];
	If[dest==False||whereto=="open"||whereto=="exit",Return[]];
	CharSize=fontsize/Imsize[Asize,Orient]; 
	BoxGrow=If[MatchQ[graphtype,LayeredGraphPlot],3,1.8];
	BoxWide=5*BoxGrow*Wider*CharSize; BoxHigh=4*BoxGrow*Taller*CharSize;	
	nb=CreateDocument[{TextCell[If[blockno== 0,"Metanetwork :"<>title,"Network Diagram for block "<>ToString[blockno]<>": "<>title]],
		Style[Framed[If[MatchQ[graphtype, GraphPlot3D],
 			GraphPlot3D[NodesnEdges, MultiedgeStyle -> False, VertexLabeling -> True,
 				 SelfLoopStyle->None, Boxed->False, PlotStyle->Arrowheads[{{Small, 0.5}}],
 				 ImageSize->Imsize[Asize,Orient],PlotRangePadding->None,AspectRatio->Aspect[Orient]],
 			graphtype[NodesnEdges,MultiedgeStyle->False,VertexLabeling->True,SelfLoopStyle->None,DirectedEdges->True,
			PackingMethod->Switch[Orient,"Tall","LayeredTop","Portrait","ClosestPacking","Landscape","LayeredLeft","Wide","LayeredLeft"],
			VertexRenderingFunction -> (If[StringMatchQ[#2[[1]], "SUBNET"~~__],
				{Hue[BlockHue[[GetBlock[#2]]], 0.15, 1], EdgeForm[Orange], Disk[#1, {4*BoxWide, BoxHigh*(1.475/Aspect[Orient])}], Text[#2, #1]} ,
				{LightYellow, EdgeForm[Orange],
					If[MatchQ[#2[[5]], Black|Purple],
						Disk[#1, {4*BoxWide, BoxHigh*(1.475/Aspect[Orient])}],
						Rectangle[#1-{4*BoxWide, BoxHigh*(1.475/Aspect[Orient])}, #1+{4*BoxWide, BoxHigh*(1.475/Aspect[Orient])}]],
				 Text[#2, #1]}
				]&),
			PlotStyle->Arrowheads[{{Small,0.5}}],ImageSize->Imsize[Asize,Orient],
			PlotRangePadding->None, AspectRatio->Aspect[Orient]]], FrameMargins->2],
		Magnification->1]},	ShowCellBracket->False];
	margins=If[Orient=="Tall"||Orient=="Portrait",{{18,18},{0,18}},{{18,18},{18,18}}];
	paperorient=If[Orient=="Tall"||Orient=="Portrait","Portrait","Landscape"];
	SetOptions[nb,PrintingOptions->{
		"PageSize"->PagePointsize[Asize,Orient],
		"PaperOrientation"->paperorient,
		"PrintingMargins"->margins,
		"PrintMultipleHorizontalPages"->True},
		PageFooters->{{None,None,None},{None,None,None}},
		PageHeaders->{{None,None,None},{None,None,None}}];
	If[whereto== "printer",SetOptions[nb,
		PrintingOptions->{"PrintMultipleHorizontalPages"->PrintMultHorPages}]];
	If[StringFreeQ[whereto,__ ~~ "3D.nb"],FrontEndExecute[FrontEndToken["SystemPrintOptionsDialog"]]];
	Which[whereto== "printer",FrontEndExecute[FrontEndToken[nb,"PrintDialog"]],StringMatchQ[whereto,__ ~~ "3D.nb"],NotebookSave[nb,whereto],True,NotebookPrint[nb,whereto]];
	NotebookClose[nb];
]

SafeDrawSubnet[blockno_,blockname_,minsize_,showmax_]:=DynamicModule[{Inpath,Paper, dest,NodesnEdges,fontsize=6,title,nb,filename,Asize,Orient},
	Status="Preparing a network layout";
	Inpath=DirectoryName[Infile];Asize="A4";Orient="Portrait";
	filename=Which[blockno==0,"NSPL_Metanet",blockno==FinalBlockCount,"NSPL_Orphan net",True,"NSPL_Subnet "<>ToString[blockno]];
	title=blockname<>" (Max connect = "<>ToString[showmax]<>")";
	NodesnEdges=GraphSpec[blockno,minsize,6,showmax];
	nb=CreateDocument[{TextCell[If[blockno== 0,"Metanetwork :"<>title,"Network Diagram for block "<>ToString[blockno]<>": "<>title]],
		Style[Framed[LayeredGraphPlot[NodesnEdges,MultiedgeStyle->False,VertexLabeling->True,SelfLoopStyle->None,
			PackingMethod->"LayeredLeft"(*"ClosestPacking"*),
			(*VertexRenderingFunction->({Text[#2,#1]}&),*) 
			PlotStyle->Arrowheads[{{Small,0.5}}],
			ImageSize->Dynamic[Imsize[Asize,Orient],SynchronousUpdating->False],
			PlotRangePadding->None,
			AspectRatio->Dynamic[Aspect[Orient],SynchronousUpdating->False]],
			FrameMargins->2],Magnification->1]},
		WindowTitle->"Network Layout",
		WindowFrameElements->{"ZoomBox","ResizeArea","CloseBox"},ShowCellBracket->False];
	SetOptions[nb,PageFooters->{{None,None,None},{None,None,None}},
		PageHeaders->{{None,None,None},{None,None,None}}];
	Paper=Grid[{{"Paper Size",RadioButtonBar[Dynamic[Asize],
		{"A4","A3","A2","A1","A0"}]},
		{"Orientation",RadioButtonBar[Dynamic[Orient],{"Tall","Portrait","Landscape","Wide"}]},
		{"Font Size",RadioButtonBar[Dynamic[fontsize],{4,5,6,7,8}]}},
		Background->LightYellow,Frame->All];
	dest=ChoiceDialog[Row[{"    ",Paper}],
		{"   PDF   \n   File"->ToFileName[Inpath,filename<>".pdf"],
		"Postscript\nFile"->ToFileName[Inpath,filename<>".ps"],
		"Printer"->"printer","Keep on Screen"->"open",
		"Cancel"->"exit"},
		WindowFloating->True,WindowMargins->{{0,Automatic},{Automatic,0}},
		WindowTitle->"Choose the layout size",ButtonBoxOptions->{ImageMargins->4}];
	SetOptions[nb,PrintingOptions->{"PageSize"->PagePointsize[Asize,Orient],
		"PrintingMargins"->{{18,18},{0,18}},
		"PrintMultipleHorizontalPages"->True}];
	If[dest== "printer",SetOptions[nb,PrintingOptions->{"PrintMultipleHorizontalPages"->PrintMultHorPages}]];
	If[dest!="exit"&&dest!="open",FrontEndExecute[FrontEndToken["SystemPrintOptionsDialog"]]];
	Which[dest== "exit",Null,dest=="open",Null,
		dest== "printer",FrontEndExecute[FrontEndToken[nb,"PrintDialog"]],
		True,NotebookPrint[nb,dest]];
	If[dest!= "open",NotebookClose[nb]];
]

(*TextSplit[text_,truncate_,line_]:=Module[{long,short}, 
	long=Min[StringLength[text],truncate];
	short=StringTake[text,long];
	StringInsert[short,"\n",Range[line+1,long,line]]]*)
	
TextSplit[text_, maxlength_, maxlines_] := 
 Module[{lines, newlines}, lines = StringSplit[text, "\n" -> "\n"];
 	newlines = Flatten@Map[
  		StringSplit[StringInsert[#,"\n", Range[maxlength + 1, StringLength[#], maxlength] ],"\n"->"\n"] &,
  	 		lines];
	StringJoin[newlines[[1 ;; Min[Length[newlines], 2*maxlines - 1]]]]]
  
CompLabel[Comp_,Peri_,LabelFontSize_]:=
	(Which[MemberQ[Peri[[4]],Comp],Tinge=Cyan,
		MemberQ[Peri[[2]],Comp],Tinge=Red,
		MemberQ[TargetIDs,Comp],Tinge=Magenta,
		MemberQ[Peri[[1]],Comp],Tinge=Green,
		MemberQ[Peri[[3]],Comp],Tinge=Blue,True,Tinge=Cyan];
	Style[TextSplit[Comp,10,3],FontSize->LabelFontSize,Bold,LineSpacing->{1.0,0},Tinge,Background->White])
	
ReactLabel[Rxn_,meta_,LabelFontSize_]:=(Tinge=
	If[(meta && StringCount[Rxn,"SUBNET"]==0),Purple,Black];
	Style[TextSplit[Rxn,10,3],FontSize->LabelFontSize,Bold,LineSpacing->{1.0,0},Tinge,Background->White])

GetBlock[LabelStyle_] := 
 ToExpression[
  StringCases[LabelStyle[[1]], "SUBNET " ~~ x__ ~~ "\n" -> x][[1]]]
  
metas[row_,cols_]:=Module[{npos,nneg},
	npos=Count[Normal@S[[row,cols]],x_/;x>0];
	nneg=Count[Normal@S[[row,cols]],x_/;x<0];
	Sign[npos]+2Sign[nneg]]

GraphSpec[Blockno_,minsize_,fontsize_,ShowThreshold_]:=
	Module[{DontShow,Blocklisting,BlockIntRows,BlockIntCols,Revno,Reversals,
		BlockRows,BlockExtRows,BlockStoi,Inflows,Outflows,Crossflows,BlockComps,
		BlockRxns,Harvest,thisrow,pos,otherblock,newblock,DropComps,
		MetaComps,MetaRows,MetaRxns,MetaRules,MetaStoi,Onelinks,MetaRxnVertices,
		DetachedOrphans,Selfloops,Periphery,MetaCount=BlockCount,LinkingCols,
		OrphanRows,OrphanCols,Orphanlinks={},OrphanRules,OrphanBlocks={},pp},
	If[Length[FinalBlocks]==0,Return[{}]];
	DontShow=Extract[Externalrows,Position[Externalrows,x_/;Connex[[x]]>ShowThreshold,{1},Heads->False]];
	If[Blockno>0,BlockIntRows=Flatten[Table[Position[Complist,FinalBlocks[[Blockno,i]]],
			{i,1,Length[FinalBlocks[[Blockno]]]}]];
		BlockIntCols=Nonzerows[Normal[Transpose[S[[BlockIntRows]]]]];
		BlockRows=Nonzerows[Normal[S[[All,BlockIntCols]]]];
		BlockExtRows=Complement[BlockRows,BlockIntRows];
		Inflows=BlockExtRows\[Intersection]InputExternalrows;
		Outflows=BlockExtRows\[Intersection]OutputExternalrows;
		Crossflows=Complement[BlockExtRows,Inflows\[Union]Outflows];
	(*Exclude high connectivity externals to avoid cluttering the graph*)
		Inflows=Complement[Inflows,DontShow];
		Outflows=Complement[Outflows,DontShow];
		Crossflows=Complement[Crossflows,DontShow];
		BlockExtRows=Join[Inflows,Crossflows,Outflows];
		BlockRows=Join[BlockIntRows,BlockExtRows];
		BlockStoi=S[[BlockRows,BlockIntCols]];
		BlockStoi=Normal@BlockStoi/.{x_/;x>0->1,x_/;x<0->2};
		Revno=Length[Reversibles];
		If[Revno>0,Reversals=Flatten@Table[Position[BlockIntCols,Reversibles[[i]]],{i,1,Revno}];
			BlockStoi[[All,Reversals]]=BlockStoi[[All,Reversals]]/.x_/;x>0->3];
		BlockComps=CompIDs[[BlockRows]];
		BlockRxns=Reactlist[[BlockIntCols]];
		Periphery={CompIDs[[Inflows]],CompIDs[[Crossflows]],CompIDs[[Outflows]],{}};
		Reap[Map[
		   (If[#[[2]]>1,Sow[CompLabel[BlockComps[[#[[1,1]]]],Periphery,fontsize]->ReactLabel[BlockRxns[[#[[1,2]]]],False,fontsize]]];
			If[OddQ[#[[2]]],Sow[ReactLabel[BlockRxns[[#[[1,2]]]],False,fontsize]->CompLabel[BlockComps[[#[[1,1]]]],Periphery,fontsize]]];)&,
		Drop[ArrayRules[BlockStoi],-1]]][[2,1]],

	(* Now the Blockno=0 case, i.e. creating the metanet. First make it all-inclusive *)
	Harvest=Reap[Do[Blocklisting=FinalBlocks[[blockno]];
		If[Length[Blocklisting]>= minsize,
			BlockIntRows=Flatten[Table[Position[Complist,Blocklisting[[i]]],{i,1,Length[Blocklisting]}]];
			BlockIntCols=Nonzerows[Normal[Transpose[S[[BlockIntRows]]]]];
			BlockRows=Nonzerows[Normal[S[[All,BlockIntCols]]]];
			BlockExtRows=Complement[BlockRows,BlockIntRows];
		(*For each external, create its connection to this block,  *)
			Map[(ToFrom=metas[#,BlockIntCols];
				If[ToFrom>0,Sow[{#,blockno}->ToFrom,link]];)&,BlockExtRows];
		(*Next, if it overlaps with internals in another block, create connections to both blocks
   			But exclude the orphan block, orphans are individually added below *)
			Do[thisrow = BlockExtRows[[i]]; pos = Position[IntrowsperBlock, thisrow]; 
 				otherblock = If[pos != {} && pos[[1, 1]] != FinalBlockCount, pos[[1, 1]], 0]; 
 				If[otherblock>0&&Length[FinalBlocks[[otherblock]]]>= minsize,
					Sow[thisrow,overlap];Sow[{thisrow,blockno}->2,link];
					Sow[{thisrow,otherblock}->2,link]],{i,1,Length[BlockExtRows]}];
		],{blockno,1,BlockCount}]];
	{MetaRules,OverlapRows}=Switch[Length[Harvest[[2]]],2,Harvest[[2]],1,{Harvest[[2,1]],{}},0,{{},{}}];
	MetaStoi=SparseArray[MetaRules,{CNodes,BlockCount}];

	(* But now compile a list of compounds to be eliminated from the metanet *)
	DropComps=Flatten[Position[Map[Nonzeroes,Normal@MetaStoi],x_/;x== 0]];
	Onelinks=Flatten[Position[Map[Nonzeroes,Normal@MetaStoi],x_/;x== 1]];
	If[minsize < 2 && IsolationBlock != {},
		OrphanRows=Flatten[Table[Position[Complist,IsolationBlock[[i]]],{i,1,Length[IsolationBlock]}]];
		OrphanCols=Nonzerows[Normal[Transpose[S[[OrphanRows]]]]];
		Orphanlinks=Nonzerows[Normal[S[[All,OrphanCols]]]],
		OrphanRows={}];
	DropComps=Union[DropComps,Complement[Onelinks,Orphanlinks],DontShow];

	(* Do the elimination and allocate remaining externals to their categories *)
	MetaStoi=Delete[MetaStoi,Split[DropComps]];
	MetaComps=Delete[CompIDs,Split[DropComps]];
	MetaRows=Complement[Range[CNodes],DropComps];
	If[MetaRows== {},
		DialogInput[DialogNotebook[{TextCell["Metanet cannot be drawn - subnets are not connected\nexcept trivially by high connectivity metabolites"],Button["OK",DialogReturn[]]}]];
		Return[{}];,
		Inflows=MetaRows\[Intersection]InputExternalrows;
		Outflows=MetaRows\[Intersection]OutputExternalrows;
		Crossflows=Complement[MetaRows,Inflows\[Union]Outflows];
		Periphery={CompIDs[[Inflows]],CompIDs[[Crossflows]],CompIDs[[Outflows]],CompIDs[[OverlapRows]]};];

	(*This section adds orphan nodes where they connect to externals that are displayed *)
	If[OrphanRows!= {},
		pp=Split[Position[Normal@S[[MetaRows,All]],x_/;x!=0],#1[[1]]== #2[[1]]&];
		LinkingCols=Table[pp[[i,j,2]],{i,1,Length[pp]},{j,1,Length[pp[[i]]]}];
		Selectlinks[orphan_,metarow_]:=Extract[LinkingCols[[metarow]],Position[Normal[S[[OrphanRows[[orphan]],LinkingCols[[metarow]]]]],x_/;x!=0]];
		Harvest=Reap[Do[newblock=True;
			Do[ToFrom=metas[MetaRows[[metarow]],Selectlinks[orphan,metarow]];
				If[ToFrom>0,newblock=If[newblock,MetaCount++;
						AppendTo[OrphanBlocks,CompIDs[[OrphanRows[[orphan]]]]];False];
					Sow[{metarow,MetaCount}->ToFrom]];,
					{metarow,Length[MetaRows]}],{orphan,Length[OrphanRows]}]];
		OrphanRules=If[Depth[Harvest]== 6,Harvest[[2,1]],{}];
		MetaRules=Join[Drop[ArrayRules[MetaStoi],-1],OrphanRules];
		MetaStoi=SparseArray[MetaRules];];

	(*The next section inserts selfloops to make Mathematica display isolated subnets*)
	MetaRxnVertices=Map[If[Length[FinalBlocks[[#]]]>= minsize,
		"SUBNET "<>ToString[#]<>"\nN="<>ToString[Length[FinalBlocks[[#]]]],"skip"]&,
		Range[BlockCount]];
	MetaRxnVertices=Select[MetaRxnVertices,#!= "skip"&];
	DetachedOrphans=Length[IsolationBlock]-Length[OrphanBlocks];
	If[DetachedOrphans>= minsize,
		AppendTo[MetaRxnVertices,"DETACHED\nORPHANS\n"<>"N="<>ToString[DetachedOrphans]]];
	Selfloops=Map[Sow[ReactLabel[#,True,fontsize]->ReactLabel[#,True,fontsize]]&,MetaRxnVertices];
	(*Finally generate the graph specification from MetaStoi*)
	MetaRxns=Join[Map["SUBNET "<>ToString[#]<>"\nN="<>ToString[Length[FinalBlocks[[#]]]]&,Range[BlockCount]],
		OrphanBlocks];
	Join[Selfloops,Reap[Map[(If[#[[2]]>1,
		Sow[CompLabel[MetaComps[[#[[1,1]]]],Periphery,fontsize]->ReactLabel[MetaRxns[[#[[1,2]]]],True,fontsize]]];If[OddQ[#[[2]]],
		Sow[ReactLabel[MetaRxns[[#[[1,2]]]],True,fontsize]->CompLabel[MetaComps[[#[[1,1]]]],Periphery,fontsize]]];)&,
		Drop[ArrayRules[MetaStoi],-1]];][[2,1]]]
	]
]

Aspect[Or_]:=Switch[Or,"Tall",2.25,"Portrait",1.475,"Landscape",0.67,"Wide",0.35]

Imsize[As_,Or_]:=Switch[Or,
	"Tall",Switch[As,"A4",660,"A3",950,"A2",1370,"A1",1950,"A0",2800],
	"Portrait",Switch[As,"A4",660,"A3",950,"A2",1370,"A1",1950,"A0",2800],
	"Landscape",Switch[As,"A4",980,"A3",1400,"A2",2000,"A1",2900,"A0",4100],
	"Wide",Switch[As,"A4",1860,"A3",2700,"A2",3940,"A1",5700,"A0",8100]]*0.8

PagePointsize[As_,Or_]:=Switch[Or,
	"Tall",Switch[As,"A4",{595,841},"A3",{841,1190},"A2",{1190,1683},"A1",{1683,2380},"A0",{2380,3366}],
	"Portrait",Switch[As,"A4",{595,841},"A3",{841,1190},"A2",{1190,1683},"A1",{1683,2380},"A0",{2380,3366}],
	"Landscape",Switch[As,"A4",{841,595},"A3",{1190,841},"A2",{1683,1190},"A1",{2380,1683},"A0",{3366,2380}],
	"Wide",Switch[As,"A4",{841,595},"A3",{1190,841},"A2",{1683,1190},"A1",{2380,1683},"A0",{3366,2380}]]

End[] (* End Private Context *)

EndPackage[]