(*
Created on 2011/01/20

Imported from
H:\My Documents\Research\Projects\SubNetExtraction\netsplitter\NetSplitter.nb
*)

(* Mathematica Package *)
BeginPackage["MainDialog`",{"NetSplitter`","FileReadWrite`","SelectExternals`","Printout`","SubnetLayout`","RandomWalkBlocker`"}]
(* Exported symbols added here with SymbolName::usage *)  
StatusWindow::usage="Displays a status window"
Splash::usage="Splash screen with copyright and contact details"
ControlPanel::usage="Displays the main user interface"

Begin["`Private`"] 
(* Begin Private Context *)


StatusWindow[]:=
	keepup=CreateDialog[{TextCell["What is happening? "],TextCell[Dynamic[Status]],TextCell["  "]},
		WindowTitle->"Netsplitter Status",WindowFloating->True,WindowMargins->{{Automatic,0},{Automatic,0}},
		WindowClickSelect->False];

Splash[logo1_,logo2_]:=DialogInput[Panel[Column[
			{TextCell[Style["NetSplitter ","Title"]],Null,TextCell[Style[NSPLVersion,"Subtitle"]],TextCell[" "],
			TextCell[Style["Copyright under GNU General Public License V3:\n\tWynand S. Verwoerd\n\tCenter for Advanced Computational Solutions \n\tLincoln University, Christchurch\n\tNew Zealand ",
				FontFamily->"Times",FontSlant->"Italic",FontWeight->"Bold",FontSize->14]],
			TextCell[Style[" "]],TextCell[Style[" Contact: wynand.verwoerd@lincoln.ac.nz  ","Subsubtitle"]],
			GraphicsRow[{logo1,DefaultButton[],logo2},ImageSize->{580,72},AspectRatio->0.2]},Alignment->Center],
		BaseStyle->Darker[Purple],ImageSize->600,Alignment->Center],
	WindowFloating->True,Background->Magenta,(*WindowSize->{615,380},*)
	WindowFrameElements->{"CloseBox"},WindowElements->{},ShowCellBracket->False];

ControlPanel[]:=
	(Infile = ToFileName[StartupDirectory,
		If[StartupDirectory == NotebookDirectory[], "Demo.tsv", "Browse.To.File"]];
	 Onthetop=Panel[Grid[{{
		"Stoichiometry matrix input file: ",InputField[Dynamic[Infile],String,FieldSize->40],
		Button["Browse..",DialogReturn[buttonchoice=7]]},
		{Tooltip["External metabolites input file: ","Optional; choose from Subnets directory to resume previous run"],
			InputField[Dynamic[Exfile],String,FieldSize->40],
		Button["Browse..",DialogReturn[buttonchoice=8]]},
		{Tooltip["Flux values input file: ","Optional"],
			InputField[Dynamic[Fluxfile],String,FieldSize->40],
		Button["Browse..",DialogReturn[buttonchoice=10]]},
		{Tooltip["Target metabolites input file: ","Optional"],
			InputField[Dynamic[Targetfile],String,FieldSize->40],
		Button["Browse..",DialogReturn[buttonchoice=9]]},
		{Tooltip[Button["Refresh Targets",DialogReturn[buttonchoice=5],Enabled->Dynamic[SplitDone]],
			"Display a new set of targets in Merge/Printout"],
		 Tooltip[Row[{"             Expand reversible reactions  ",Checkbox[Dynamic[ReverseThem]]}],
		 	"Reversible reactions will be duplicated in direction opposite to input specification"],
		 Tooltip[Button["Convert input",DialogReturn[buttonchoice=11]],
		 	"Only converts input file to alternative format (.tsv / .sbml )"]}},
		Spacings->{0,2}],Style["Input Selection","Subtitle"]];
	Inthemid=Panel[Grid[{{
		Tooltip["Omit reactions with ID containing:","E.g., to leave out biomass reactions. Alternatives can optionally be given in a format like BIO|BC|Biom (no spaces) "],
			InputField[Dynamic[Dropstring],String,FieldSize->30]},
		{"Vector similarity measure: ",PopupMenu[Dynamic[DistChoice],{DiceDissimilarity->"Dice",JaccardDissimilarity->"Jaccard",SokalSneathDissimilarity->"Sokal-Sneath"},SokalSneathDissimilarity]},
		{"Intergroup distance: ",PopupMenu[Dynamic[LinkChoice],{"Single"->"Smallest","Average","Complete"->"Largest","Weighted","Centroid","Median","Ward"},"Single"]},
		{"Max connectivity for internal metabolites:",InputField[Dynamic[ConnexMax],Number,FieldSize->3]}}],
		Style["Algorithmic Options","Subtitle"]];
	suborfull=Dynamic[If[SplitDone,"Subnet Diagram","Full Network"]];
	Onthebot=Panel[Grid[{{
		DefaultButton["Split Network",DialogReturn[buttonchoice=1]],
		Button["Merge Subnets",DialogReturn[buttonchoice=2],Enabled->Dynamic[SplitDone]],
		Button["Save  Subnets",DialogReturn[buttonchoice=3],Enabled->Dynamic[SplitDone]],
		Button["Printed Output",DialogReturn[buttonchoice=4],Enabled->Dynamic[SplitDone]],
		Button[Dynamic[suborfull],DialogReturn[buttonchoice=6]],
		CancelButton["Exit",DialogReturn[buttonchoice=12]]}},Spacings->2],
		Style["Actions","Subtitle"]];
	buttonchoice = 1;
	While[buttonchoice<12,CalculationTime+=TimeUsed[]-StartTime;
		Status="Waiting for action";(*EmitSound[Attention];*)
		DialogInput[Column[{Onthetop,Inthemid,Onthebot},Spacer[20]],WindowTitle->"NetSplitter Control Panel ",WindowSize->All];
		StartTime=TimeUsed[];
		Switch[buttonchoice,
			1,CalculationTime=0;If[GetInputS[],ExternalSearch[];],
			2,Combine[{}],
			3,SubExport[],
			4,Printout[],
			5,TargetCompounds={};Targetfile=If[StringLength[Targetfile]== 0," ",Targetfile];
				Targetmess=If[ToString[FileType[Targetfile]]== "File",GetTargets[Targetfile],
					"\nNo target metabolites input file was read."];,
			6,If[SplitDone,CheckedInput=False;
				(*AscendList=Sort@Map[Length,Most@FinalBlocks];metanodes=30;
				minblocksize=If[Length[AscendList]>metanodes,AscendList[[-metanodes]],1];*)
				minblocksize=1;showmax=ConnexMax;
				Blockchoice=Pick[Range[FinalBlockCount],Map[Length[#]>0&,FinalBlocks]];PrependTo[Blockchoice,0];
				While[!CheckedInput,Descript=ConstantArray["",BlockCount];
					PrependTo[Descript,"Subnet linking by shared externals"];
					AppendTo[Descript,"Subnets of size 1 (orphans) "];
					Choicetext=#->"Subnet "<>ToString[#]&/@Blockchoice;
					Choicetext=ReplacePart[Choicetext,1->0->"Metanet"];
					If[IsolationBlock!={},Choicetext=ReplacePart[Choicetext,-1->FinalBlockCount->"Orphans"]];
					res=DialogInput[{name=Descript[[1]],blockno=0},Grid[{{
						TextCell["Block selection"],
						TextCell["Subnetwork Title (optional)"],
						TextCell["Smallest Metanet block"],
						InputField[Dynamic[minblocksize],Number,FieldSize->2]},
						{PopupMenu[Dynamic[blockno,(blockno=#;name=If[#<= FinalBlockCount,Descript[[#+1]],Descript[[2]]])&],Choicetext,FieldSize->8],
						InputField[Dynamic[name],String,FieldSize->40],
						TextCell["Max Connectivity"],
							InputField[Dynamic[showmax],Number,FieldSize->2]},{TextCell[""],TextCell[""],
								DefaultButton["Draw subnetwork",DialogReturn[{blockno,name,minblocksize,showmax}]],
								CancelButton[]}},
						Spacings->2],
					WindowTitle->"Choose the block to be drawn"];
					If[Length[res]== 4,Which[
						Not[0<= res[[1]]<= FinalBlockCount],
							DialogInput[DialogNotebook[{TextCell["Invalid block number: "<>ToString[res[[1]]]<>"\nMust be between 0 and  "<>ToString[FinalBlockCount]],
								Button["OK",DialogReturn[1]]}]];,
						res[[1]]== 0&&Count[Map[Length,Most@FinalBlocks],x_/;x>res[[3]]]<2,
							DialogInput[DialogNotebook[{TextCell["Cannot draw meta-network for fewer than two blocks. \nTry reducing minimum block size"],
								Button["OK",DialogReturn[1]]}]];,
						True,
							CheckedInput=True;
							If[res[[1]]== 0,res[[3]]=minblocksize=Max[minblocksize,1];res[[2]]=res[[2]]<>", showing blocks larger than "<>ToString[minblocksize-1]];
							If[SafeDraw,SafeDrawSubnet@@res,DrawSubnet@@res];],
					CheckedInput=True;]],
				buttonchoice=1;GetInputS[];DrawFullnet[]],
			7,If[ToUpperCase[Infile]=="DEMO",Infile=ToFileName[NotebookDirectory[],"Demo.tsv"]];
				ret=SystemDialogInput["FileOpen",{Infile,{"Network "->{"*.sbml","*.xml","*.bpt","*.tsv"},"All"->{"*.*"}}}];
				If[ret!= "$Canceled",Infile=ret;Exfile=DirectoryName[Infile];Targetfile=DirectoryName[Infile];Fluxfile=DirectoryName[Infile];];SplitDone=False;,
			8,ret=SystemDialogInput["FileOpen",{Exfile,{"Text "->{"*.txt","*.tsv","*."},"All"->{"*.*"}}}];
				If[ret!= "$Canceled",Exfile=ret];,
			9,ret=SystemDialogInput["FileOpen",{Targetfile,{"Text "->{"*.txt","*.tsv","*."},"All"->{"*.*"}}}];
				If[ret!= "$Canceled",Targetfile=ret];,
			10,ret=SystemDialogInput["FileOpen",{Fluxfile,{"Spreadsheet "->{"*.xls","*.dif","*.csv"},"Text "->{"*.txt","*.tsv"},"All"->{"*.*"}}}];
				If[ret!= "$Canceled",Fluxfile=ret];,
			11,Convert=True;GetInputS[];Convert=False;]];
	NotebookClose[keepup];)



End[] 
(* End Private Context *)

EndPackage[]
