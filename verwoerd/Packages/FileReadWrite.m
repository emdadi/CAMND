
(* Mathematica Package *)

(* COPYRIGHT (C) 2010 Wynand Verwoerd 
under the terms of the GNU General Public License V3
This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
*)

BeginPackage["FileReadWrite`",{"NetSplitter`"}]
(* Exported symbols added here with SymbolName::usage *)  
Nonzeroes::usage="Nonzeroes[y] gives the number of non-zero entries in vector y"
Zerows::usage="Zerows[A] gives a list of zero rows in matrix A"
Nonzerows::usage="Nonzerows[A] gives a list of non-zero rows in matrix A"
LastWord::usage="Utility to extract the last whitespace separated word from a string"
String2Number::usage="Utility to extract number, including e.g. 3.4E-03"
GetInputS::usage="GetInputS calls readers for input files and digests the results"
GetExternals::usage="GetExternals[Exfile] reads external compounds"
GetTargets::usage="GetTargets[Targetfile] reads target compounds"
SubExport::usage="SubEport creates individual subnets and controls their storage in separate files"
TSVreader::usage="Reads text format stoichiometry input files" 
BIOPTreader::usage="Reads BIOPT format stoichiometry input files"
SBMLreader::usage="Reads stoichiometry input files in Sys Bio Markup Lang format"
TSVExport::usage="Writes text file specification for a single subnet"
SBMLExport::usage="Writes SBML file specification for a single subnet"

(* Variables only shared with the Printout package, defined here to avoid cluttering
   the list of generally shared "semi-global" variables defined in the Netsplitter` context *)
Fluxmess="";Externalmess;Targetmess;ModelNote;
RemovedRXN;RemovedComp;DropRXN;
Reversecount;Extras;
(*fulltext;FluxUnits;FluxBounds;ObjWeights;FluxValues;FluxSigns;*)

Begin["`Private`"] (* Begin Private Context *) 

Nonzeroes[y_]:=Count[y,x_/;x!= 0]

Zerows[A_]:=Flatten[Position[Map[Nonzeroes,A],x_/;x==0]]

Nonzerows[A_]:=Flatten[Position[Map[Nonzeroes,A],x_/;x!= 0]]

String2Number[string_]:=If[StringMatchQ[StringTrim[string],NumberString ~~ ("E" | "e" ~~ NumberString) ...], 
 ToExpression[StringReplace[string, "E" | "e" -> "*10^"]], False]
      
LastWord[stringarr_]:=Map[Reverse,StringReverse[StringSplit[StringReverse[stringarr],Whitespace,2]],{Depth[stringarr]-1}]

GetInputS[]:=Module[{Extens,duplis,duplicates,duplisat={},NewIDs,Extrarows,
	MetabCompartlist,Expandedlist,ncomp,nreact,nrows,ncols,ReadFluxes,flip,noflip,rno,
	shortlist,suflist,compartmax,cpartid,falsereverse,nabsent,ReadOK},Off[Remove::"remal"];
	SplitDone=False;Status="Processing the input files";StartTime=TimeUsed[];Fluxes={};Fluxmess="";
	Extens=If[ToUpperCase[Infile]=="DEMO",Infile=ToFileName[NotebookDirectory[],"Demo.tsv"];"TSV",ToUpperCase@StringDrop[Infile,{1,Last[StringPosition[Infile,"."]][[1]]}]];
	ReadOK=Switch[Extens,"TSV",TSVreader[Infile],"BPT",BIOPTreader[Infile],"SBML",SBMLreader[Infile],"XML",SBMLreader[Infile],
		_,EmitSound[Crisis];
			DialogInput[DialogNotebook[{TextCell["Unable to read input file \nOnly .SBML, .XML or .TSV files can be processed \nUse Browse button to find valid input."],
			Button["OK",DialogReturn[1]]}]];False];
	If[ReadOK,DialogInput[DialogNotebook[{TextCell["Successfully read ID: "<>ModelID<>
		"\nModel name: "<>ModelName<>"\n\nNotes on the model: \n"<>ModelNote,PageWidth -> 400],
			Button["OK",DialogReturn[1]]},WindowTitle->"Input Model"]];,Return[False]];

(* Compound/Reaction numbers converted to strings, embedded quotes removed *)
	Complist=ToString/@Complist;Complist=Map[StringReplace[#,"\""->""]&,Complist];
	StoichioExternals=ToString/@StoichioExternals;StoichioExternals=Map[StringReplace[#,"\""->""]&,StoichioExternals];
	CompNames=Map[StringReplace[#,"\""->""]&,CompNames];
	Reactlist=ToString/@Reactlist;Reactlist=Map[StringReplace[#,"\""->""]&,Reactlist];
	Reverselist=ToString/@Reverselist;Reverselist=Map[StringReplace[#,"\""->""]&,Reverselist];
(* Truncate long names on the compound list *)
	Complist=Replace[Complist,x_/;StringLength[x]>50:> StringTake[x,49]<>" #"(*<>ToString[n++]*),1];

(* Do a few data integrity checks *)
	ncomp=Length[Complist]; nreact=Length[Reactlist]; nrows=Dimensions[S][[1]];ncols=Dimensions[S][[2]];
	duplis=Length[Complist]-Length@DeleteDuplicates[Complist];
	If[ncomp==0||nreact==0||nrows==0||ncols==0,EmitSound[Crisis];DialogInput[DialogNotebook[{TextCell[
		"Invalid data : Input file\n"<>Infile<>"\ndid not yield a metabolite list, reaction list and stoichiometry matrix."],
		Button["OK",DialogReturn[1]]}]];Return[False]];
	If[duplis > 0, duplicates = Map[#[[1]] &, Select[Tally[Complist], #[[2]] > 1 &]];
		If[Total[StringCount[duplicates, " #"]] > 1,
		(*	Revert duplicates created by truncation to the full ID *)
			duplisat = Flatten@Map[Position[Complist, #] &, duplicates];
			Complist[[duplisat]] = CompIDs[[duplisat]];
			StoichioExternals = DeleteCases[StoichioExternals, x_ /; StringMatchQ[x, CompIDs[[duplisat]]~~" "~~___]];
		(* Remove troublemakers before truncating names in the s. externals, 
			then below just CompIDS are added back for them  *)
			duplis = Length[Complist] - Length@DeleteDuplicates[Complist];]
	];
(* Truncate long names for externals as well, but avoid any duplicate found above *)
	StoichioExternals=Replace[StoichioExternals,x_/;StringLength[x]>50:>StringTake[x,49]<>" #"(*<>ToString[n++]*),1];
	StoichioExternals=Join[StoichioExternals,CompIDs[[duplisat]]];
	(* Now repeat duplicates test to check if reversion solved the problem *)
	If[duplis>0,EmitSound[Crisis];
		DialogInput[DialogNotebook[{TextCell["Invalid data - Duplicate metabolite IDs: 
			\n"<>ToString@Take[Commonest[Complist],Min[3,duplis]]<>" 
			\netc. were found. Check that metabolites listed as external 
			\nin the stoichiometry file are not listed under internals as well!"],
			Button["OK",DialogReturn[1]]}]];Return[False]];
	duplis=Length[Reactlist]-Length@DeleteDuplicates[Reactlist];
	If[duplis>0,EmitSound[Crisis];DialogInput[DialogNotebook[{TextCell["Invalid data - Duplicate reaction IDs: 
		\n"<>ToString@Take[Commonest[Reactlist],Min[3,duplis]]<>" 
		\netc. were found."],Button["OK",DialogReturn[1]]}]];Return[False]];
	If[ncomp!= nrows,EmitSound[Crisis];DialogInput[DialogNotebook[{TextCell[
		"Invalid data : "<>ToString[ncomp]<>" metabolite names listed but "<>ToString[nrows]<>" rows in S matrix "],
		Button["OK",DialogReturn[1]]}]];Return[False]];
	If[nreact!= ncols,EmitSound[Crisis];DialogInput[DialogNotebook[{TextCell[
		"Invalid data : "<>ToString[nreact]<>" reaction names listed but "<>ToString[ncols]<>" columns in S matrix "],
		Button["OK",DialogReturn[1]]}]];Return[False]];
	falsereverse=Complement[Reverselist,Reactlist];nabsent=Length[falsereverse];
	If[nabsent>0,EmitSound[Crisis];DialogInput[DialogNotebook[{TextCell[
		"Invalid data - listed reversibles "<>ToString@Take[falsereverse,Min[3,nabsent]]<>" etc. not found on the full reaction list "],
		Button["OK",DialogReturn[1]]}]];Return[False]];

(* CreateDocument[MatrixPlot[S]]; *)

	Newlist={};Refusals={};Stepcount=0;RepeatMerge=False;OldMerge={};
	Exfile=If[StringLength[Exfile]== 0," ",Exfile];
	Externalmess=If[ToString[FileType[Exfile]]== "File",
		GetExternals[Exfile];
		"\n"<>ToString[Length[Newlist]]<>" external metabolites read from input file: "<>Exfile,"\nNo external metabolites input file was read."];

	RemovedRXN={};RemovedComp={};DropComps={};
	DropRXN=Position[Reactlist,x_/;StringMatchQ[ToString[x],Map[___~~#~~___&,StringSplit[Dropstring,"|"]]],{1}];
	If[DropRXN!= {},
		RemovedRXN=Reactlist[[Flatten[DropRXN]]];
		Reactlist=Delete[Reactlist,DropRXN];
		If[ReactionXML!={},ReactionXML=Delete[ReactionXML,DropRXN]];
		S=Transpose@Delete[Transpose[S],DropRXN];
		DropComps=Split@Zerows[S];
		S=Delete[S,DropComps];
		RemovedComp=Complist[[Flatten[DropComps]]];
		Complist=Delete[Complist,DropComps];
		CompIDs=Delete[CompIDs,DropComps];
		CompNames=Delete[CompNames,DropComps];
		Compartlist=Delete[Compartlist,DropComps];
		StoichioExternals=StoichioExternals\[Intersection]Complist;
		Reverselist=Reverselist\[Intersection]Reactlist];
	Reversecount=Length[Reverselist];
	Reversibles=Flatten@Map[Position[Reactlist,#]&,Reverselist];
	(*Keep a copy to use when saving subnets
	OriginalReversibles=Reversibles;OriginalReverselist=Reverselist;*)
	
	If[CompartIDs == {DefaultCompartment},
		suflist = StringCases[CompIDs, __~~"_"~~suffix:WordCharacter.. -> suffix]; 
  		shortlist = SortBy[Cases[Tally@Flatten[suflist], {_, _?(# > 2 &)}], Last];
		compartmax = Min[12, Length[shortlist]]; 
		If[compartmax>0,
		shortlist = Take[shortlist, {-1, -compartmax, -1}][[All, 1]]; 
		cpartid = DialogInput[{cpartid = {}}, 
			Column[{"No explicit compartment specs found.\n" <> 
       		"Multiple occurrences of the following suffixes on\n"<>
       		"metabolite IDs might indicate compartment allocations:", 
			CheckboxBar[Dynamic[cpartid], shortlist, Appearance -> "Vertical"], 
			"\nChoose any by ticking, or click Cancel to exit.",
			Row[{DefaultButton["Proceed", DialogReturn[cpartid]], 
         		Button["Cancel", DialogReturn[{}]]}]}]];
  		If[cpartid != {}, 
			suflist = Flatten[suflist /. {} -> DefaultCompartment];
			Compartlist = Replace[Flatten[suflist], {x_/;!MemberQ[cpartid, x] -> DefaultCompartment}, 1]; 
   			Complist = MapThread[StringJoin, {CompIDs, ConstantArray[" ", Length[CompIDs]], Compartlist}]];
   			CompartIDs=DeleteDuplicates[Compartlist];
			CompartOuters=ConstantArray["Missing",Length[CompartIDs]];
			CompartSizes=ConstantArray["1.0",Length[CompartIDs]];
			CompartUnits=ConstantArray["volume",Length[CompartIDs]];
			CompartConstant=ConstantArray["Missing",Length[CompartIDs]];
			CompartDims=ConstantArray["3",Length[CompartIDs]];
	]];

	If[Convert,InputExternalrows={};OutputExternalrows={};
	Inpath=DirectoryName[Infile];CurDir=Directory[];SetDirectory[Inpath];
	FileTitle={ModelID<>"_Converted",ModelName<>" Converted from "<>Extens<>" by Netsplitter."};
	Status="Saving export file";
	(* SBMLExport is normally given a list of internal metabolites only. To do that use
	Switch[Extens,
		"TSV",SBMLExport[FileTitle,Complement[Complist,StoichioExternals]],
		"BPT",SBMLExport[FileTitle,Complement[Complist,StoichioExternals]],
		"SBML",TSVExport[FileTitle,Complement[Complist,StoichioExternals]]];
	However, one-sided reactions get left out this way if their single reactant
	is also declared as external, because there is then nothing to connect it 
	to the rest of the network. 
	To avoid that when converting an input file and preserve fidelity in SBML,
	call the function giving the full list. That seems to work OK, but watch this in
	case of problems.	
	*)
	Switch[Extens,
		"TSV",SBMLExport[FileTitle,Complist],
		"BPT",SBMLExport[FileTitle,Complist],
		"SBML",TSVExport[FileTitle,Complist]];
	SetDirectory[CurDir];Return[True]];

	Fluxfile=If[StringLength[Fluxfile]== 0," ",Fluxfile];
	DoFlux=If[ToString[FileType[Fluxfile]]== "File",If[StringMatchQ[Infile,___~~"Block"~~NumberString~~"."~~__],
		EmitSound[Crisis];
		ChoiceDialog[Column[{TextCell["ARE YOU SURE THE FLUX FILE SHOULD BE APPLIED?"],
			TextCell["\nThe input file name "<> Infile<>" corresponds to a subnet file created by Netsplitter in a previous run."],
			TextCell["\nIf so, flux information used in that run is already incorporated and repeating it will cause chaos."]}],
		{"Ignore flux file"->False,"Load flux file"->True},
		WindowFloating->True,ButtonBoxOptions->{ImageMargins->20}],True],False];
	Fluxmess=Fluxmess<>If[DoFlux,
	  ReadFluxes=Which[
		StringCount[ToUpperCase[Fluxfile],".XLS"]== 1,Import[Fluxfile][[1]],
		StringCount[ToUpperCase[Fluxfile],".DIF"]== 1,Import[Fluxfile],
		StringCount[ToUpperCase[Fluxfile],".CSV"]== 1,Import[Fluxfile],
		StringCount[ToUpperCase[Fluxfile],".TXT"]== 1,Import[Fluxfile,"Table"],
		StringCount[ToUpperCase[Fluxfile],".TSV"]== 1,Import[Fluxfile,"Table"],
		True,{{"Invalid ","format"}}][[All,1;;2]];
	  If[ReadFluxes== {{"Invalid ","format"}},"\nWARNING - unable to process flux file format !",
	  	FluxValues=ReadFluxes;
	  	ReadFluxes = Map[MapAt[Sign, #, 2] &, FluxValues];
	  	FluxSigns=DeleteDuplicates[Join[ReadFluxes,FluxSigns], SameQ[#1[[1]], #2[[1]]] &];
	  	"\n"<>ToString[Length[ReadFluxes]]<>" Flux values were read from file "<>Fluxfile],
			"\nNo flux values input file was read."];
	flip=0;noflip=0;	
	If[FluxSigns != {},Reversibles = Flatten@Map[Position[Reactlist, #] &, Reverselist]; 
 		Map[(If[MemberQ[Reactlist, #[[1]]], rno = Position[Reactlist, #[[1]]][[1, 1]];
		Which[#[[2]] < 0, S[[All, rno]] = -S[[All, rno]]; Reversibles = DeleteCases[Reversibles, rno]; flip++;,
			  #[[2]] > 0, Reversibles = DeleteCases[Reversibles, rno]; noflip++;]])&, FluxSigns];
		Reverselist = Reactlist[[Reversibles]];
	Fluxmess = 	Fluxmess <> 
				"\nAs a result, " <> ToString[noflip] <> " reaction directions were confirmed and " <>
				 ToString[flip] <> " were reversed, leaving " <> ToString[Length[Reversibles]] <> 
   				" as undetermined reversibles."];

CNodes=Dimensions[S][[1]];RNodes=Dimensions[S][[2]];
Externals=StoichioExternals;
Externalrows=Flatten[Table[Position[Complist,Externals[[i]]],{i,1,Length[Externals]}]];StExrows=Externalrows;



CR=Abs[Sign[S]]-Sign[S];
RC=Transpose[Abs[Sign[S]]+Sign[S]];
(*Append duplicate rows/cols for reversible reactions *) If[ReverseThem,
RC=Join[RC,Transpose@CR[[All,Reversibles]]];CR=Transpose@Join[Transpose[CR],RC[[Reversibles]]]];
InternalComplist=Complist;InternalCNodes=CNodes;
InputExternalrows=Flatten[Position[Total[RC],0]]; OutputExternalrows=Flatten[Position[Total[Transpose[CR]],0]];
StructuralExternalrows=Union[InputExternalrows,OutputExternalrows];Externalrows=Union[Externalrows,StructuralExternalrows];

Connex=Total[Transpose[Abs[Sign[S]]]];ConnectionExternalrows=Ordering[Connex,Count[Connex,x_/;x>ConnexMax],Greater];
If[Refusals!= {},Refusedrows=Flatten[Table[Position[Complist,Refusals[[i]]],{i,1,Length[Refusals]}]];
	ConnectionExternalrows=DeleteCases[ConnectionExternalrows,x_/;MemberQ[Refusedrows,x]]];
BasicExternals=Complist[[Union[StructuralExternalrows,ConnectionExternalrows]]];
Externalrows=Union[Externalrows,ConnectionExternalrows];
Removedrows=Externalrows;Externals=Complist[[Externalrows]];

(* External compound spec (ID Compartment) from file is expanded to a list 
	by regular expression matching to the full list of Id's and compartments *)
MetabCompartlist = MapThread[StringJoin,{CompIDs,ConstantArray[" ",Length[CompIDs]],Compartlist}];
Expandedlist=Cases[MetabCompartlist, x_ /; StringMatchQ[x, RegularExpression/@Newlist]];
NoMatch=Complement[Newlist,Expandedlist,  SameTest->(StringMatchQ[#2, RegularExpression[#1]] &)];
Newlist=Expandedlist;
Expandedlist=Cases[MetabCompartlist, x_ /; StringMatchQ[x, RegularExpression/@Refusals]];
NoMatch=Join[NoMatch,Complement[Refusals,Expandedlist,  SameTest->(StringMatchQ[#2, RegularExpression[#1]] &)]];
Refusals=Expandedlist;
Newlist = Complement[Newlist, Refusals];
Extras = {}; Extrarows = {};
If[Newlist != {}, NewIDs = StringSplit[Newlist, Whitespace, 2][[All, 1]];
  Extrarows = Flatten[Table[Position[CompIDs, NewIDs[[i]]], {i, 1, Length[NewIDs]}]];
  Extrarows = Complement[Extrarows, Externalrows];
  Externalrows = Join[Externalrows, Extrarows];
  Removedrows = Externalrows; Externals = Complist[[Externalrows]]];
Extras=Complist[[Extrarows]];

(* The following code can be uncommented to force metabolites read from the externals file to
be taken as external and stop them from being reincorporated (for debugging purposes.)
This is not easily done in an existing Smatrix.TSV input file because it would require reordering the 
matrix so the newly classified externals come last. For SBML input there is no problem, simply add
"BoundaryCondition=True" to any metabolites that are forced to be external  *)

(*StoichioExternals=Union[StoichioExternals,Extras]; Extras={};*)


TargetCompounds={};Targetfile=If[StringLength[Targetfile]== 0," ",Targetfile];
Targetmess=If[ToString[FileType[Targetfile]]== "File",GetTargets[Targetfile],"\nNo target metabolites input file was read."];

Return[True]
]

GetExternals[Exfile_]:=Module[{fulltext,Nlines,ReadConfig,Configline,configno,config,
	ReadSteps,RSteps,ReadMerge,Exclusion,takeit,i},
	fulltext=Import[Exfile,"Lines"];Nlines=Length[fulltext];
(*Read the configuration from Externals file*)
	ReadConfig=False;
	Do[If[StringMatchQ[fulltext[[i]],"%%%%%%%%\tConfiguration\t%%%%%%%%"],ReadConfig=True];
		If[ReadConfig&&StringMatchQ[fulltext[[i]],"% "~~NumberString~~__],
			Configline=StringReplace[fulltext[[i]],"\t%"~~__->""];
			Configline=StringSplit[Configline,Whitespace,3];
			configno=ToExpression[Configline[[2]]];
			config=If[Length[Configline]== 3,Configline[[3]],""];
			Switch[configno,
				1,ReverseThem=ToExpression[config];,
				2,Dropstring=config;,
				3,DistChoice=ToExpression[config],
				4,LinkChoice=ToExpression[config],
				5,ConnexMax=ToExpression[config],
				6,Fluxfile=config]];
		If[StringMatchQ[fulltext[[i]],"%%%%%%%%\tChoices\t%%%%%%%%"],
	Break[]],{i,1,Nlines}];

(*Read the merge history from Externals file*)
	ReadMerge=False;
	Do[If[StringMatchQ[fulltext[[i]],"%%%%%%%%\tMerges\t%%%%%%%%"],ReadMerge=True];
		If[ReadMerge&&StringMatchQ[fulltext[[i]],"% "~~NumberString~~__],
			AppendTo[OldMerge,ToExpression[StringSplit[fulltext[[i]],"\t"][[2;;3]]]];
			If[StringMatchQ[fulltext[[i]],StartOfString~~Except["%"]~~___],Break[]]
	],{i,1,Nlines}];
	OldMerge=Map[{If[Length[#[[1]]]>0,#[[1]],{}],If[Length[#[[2]]]>0,#[[2]],{}]}&,OldMerge];
	If[OldMerge!= {},ReadSteps=Length[OldMerge];RSteps = Null;
		EmitSound[Attention];
		RepeatMerge=ChoiceDialog[Column[{
			TextCell["A merge history of " <> ToString[ReadSteps] <> " steps was read from the externals file."], 
   			TextCell["\nSteps up to the following number are reapplied to this run."], 
			InputField[Dynamic[RSteps], Number, FieldSize -> 12, 
    			FieldHint -> ToString[ReadSteps] <> "       (Click to change)", 
    			FieldHintStyle -> {Red, Opacity[1], Plain}], 
			TextCell["\nChoosing None allows externals selection to be continued.
					\nChoosing OK but with merge steps = 0 skips selection to perform a new merging cycle."]
    	  }],{"\tOK\t" -> True, "\tNone\t" -> False}, 
    	WindowFloating -> True, ButtonBoxOptions -> {ImageMargins -> 20}];
	If[ToString[RSteps] != "Null", ReadSteps = Max[0, Min[RSteps, ReadSteps]]];
	If[RepeatMerge,OldMerge=OldMerge[[1;;ReadSteps]]]
	];

(*Read the externals and refusals from Externals file.
  Keep only first 2 fields, and put .* to match any compartment if the 2nd is missing.*)
	Exclusion=False;
	For[i=1,i<= Nlines,i++,
		takeit=StringDrop[fulltext[[i]],Flatten[StringPosition[fulltext[[i]],"%"~~___,Overlaps->False],1]];
		takeit=StringTrim[takeit];
		If[takeit!="",takeit=StringJoin[Riffle[PadRight[StringSplit[takeit, Whitespace], 2, ".*"][[{1, 2}]], " "]]];
		If[ToUpperCase[takeit]== "<REFUSALS> .*",Exclusion=True;takeit=""];
		If[StringLength[takeit]>0,If[Exclusion,AppendTo[Refusals,takeit],AppendTo[Newlist,takeit]]]];

]

GetTargets[Targetfile_]:=Module[{fulltext,takeit,Newlist={},MetabCompartlist,Targetrows={},i},
	fulltext=Import[Targetfile,"Lines"];
For[i=1,i<= Dimensions[fulltext][[1]],i++,takeit=StringDrop[fulltext[[i]],Flatten[StringPosition[fulltext[[i]],"%"~~___,Overlaps->False],1]];
	If[takeit!="",takeit=StringJoin[Riffle[PadRight[StringSplit[takeit, Whitespace], 2, ".*"][[{1, 2}]], " "]]];
	takeit=StringTrim[takeit];
	If[StringLength[takeit]>0,AppendTo[Newlist,takeit]]];
	TargetIDs={};TargetCompounds={};
	If[Newlist!= {},
		MetabCompartlist = MapThread[StringJoin,{CompIDs,ConstantArray[" ",Length[CompIDs]],Compartlist}];
		Newlist=Cases[MetabCompartlist, x_ /; StringMatchQ[x, RegularExpression/@Newlist]];
		Targetrows=Flatten[Table[Position[MetabCompartlist,Newlist[[i]]],{i,1,Length[Newlist]}]]];
	TargetCompounds=Complist[[Targetrows]];TargetIDs=CompIDs[[Targetrows]];
Return["\n"<>ToString[Length[TargetCompounds]]<>" target metabolites read from input file: "<>Targetfile]
]

SubExport[]:=Module[{CurDir,IDlist,rowlist,Extlisting,Refuselisting,Configuration,DecisionHistory,MergeHistory},Status="Exporting subnetwork files";Off[DeleteDirectory::"nodir"];Off[DeleteDirectory::"dirne"];
Inpath=DirectoryName[Infile];CurDir=Directory[];SetDirectory[Inpath];
DeleteDirectory["Subnetworks",DeleteContents->True];
CreateDirectory["Subnetworks"];SetDirectory["Subnetworks"];
PreMergeExternals=Complement[PreMergeExternals,BasicExternals];
IDlist = StringSplit[PreMergeExternals, Whitespace, 2][[All, 1]];
rowlist=Flatten[Table[Position[CompIDs, IDlist[[i]]],{i,1,Length[IDlist]}]];
Extlisting =  MapThread[StringJoin, 
	{IDlist, ConstantArray[" ",Length[IDlist]], Compartlist[[rowlist]],
	 ConstantArray["\t\t% ", Length[IDlist]], CompNames[[rowlist]]}];
IDlist = StringSplit[Refusals, Whitespace, 2][[All, 1]];
rowlist=Flatten[Table[Position[CompIDs, IDlist[[i]]],{i,1,Length[IDlist]}]];
Refuselisting =  MapThread[StringJoin, 
	{Refusals, ConstantArray["\t\t% ", Length[rowlist]], CompNames[[rowlist]]}];	 
Configuration={StringJoin[{"%%%%%%%%\tConfiguration\t%%%%%%%%\n",
							"% 1 "<> ToString[ReverseThem] <>"\t%(Reversible reactions duplicated)\n",
							"% 2 "<> Dropstring<>"\t%(Reaction omission)\n",
							"% 3 "<>ToString[DistChoice]<>"\t%(Distance function)\n",
							"% 4 "<>ToString[LinkChoice]<>"\t%(Intercluster distance)\n",
							"% 5 "<>ToString[ConnexMax]<>"\t%(Maximal internal connectivity)\n",
							"% 6 "<>Fluxfile<>"\t%(Flux values path)"}]};
DecisionHistory={StringJoin[Prepend[Riffle[Map[Riffle[#,"\t"]&,
	Map[ToString,DecisionTrace,{2}]],"\n% "],"%%%%%%%%\tChoices\t%%%%%%%%\n% "]]};
MergeHistory={StringJoin[Prepend[Riffle[Map[Riffle[#,"\t"]&,
	Map[ToString,MergeTrace,{2}]],"\n% "],"%%%%%%%%\tMerges\t%%%%%%%%\n% "]]};
Export[Outfile,Join[{"%%% List of external compounds to recreate the subnetworks found from",
	"%%% Input file "<>Infile<>" calculated at "<>DateString[]},
	Configuration,DecisionHistory,MergeHistory,Extlisting,
	{"<Refusals>"},Refuselisting],{"Text","Lines"}];
If[FinalBlockCount>1,
	Blocktitle=Table[{ModelID<>"_Block_"<>ToString[i],"Subnet of "<>ModelName},{i,1,FinalBlockCount}];
	Blocktitle[[-1]]={ModelID<>"_Orphan_Block","Minor subnets of "<>ModelName};
	dest=ChoiceDialog["Which output file format?",{"SBML"->"sbml","Text file"->"tsv","Cancel"->"exit"},WindowFloating->True,ButtonBoxOptions->{ImageMargins->6}];
	Which[dest== "exit",Null,
		  dest== "tsv",MapThread[TSVExport,{Blocktitle,FinalBlocks}],
		  dest== "sbml",MapThread[SBMLExport,{Blocktitle,FinalBlocks}]];];
SetDirectory[CurDir];
]
      
LastWord[stringarr_]:=Map[Reverse,StringReverse[StringSplit[StringReverse[stringarr],Whitespace,2]],{Depth[stringarr]-1}]

TSVreader[filepath_]:=Module[{fulltext,Splitlist,InputName,fullname},
	fulltext=StringSplit[Import[filepath,"String"],"<"~~LetterCharacter..~~">"~~Whitespace];
	InputName=ImportString[fulltext[[1]],"TSV"][[1,1]];
	InputName=StringTrim[InputName]; If[InputName=="",InputName="NoName"];
	fullname = StringSplit[InputName, Whitespace, 2]; 
	{ModelID, ModelName} = If[Length[fullname] == 1, {InputName, InputName}, fullname];
	ModelNote="TSV input file - no notes.";
	FuncDefs = {};UnitDefs = {};CompartTypes = {};SpecieTypes = {};
	Parameters = {};InitAss = {};Rules = {};Constraints = {};Events = {};
	ReactionXML = {}; SpeciesXML={}; 
	FluxValues={}; FluxSigns={}; FluxBounds={}; ObjWeights = {};
	Reactlist=Flatten[ImportString[fulltext[[2]],"TSV"]];
	Reverselist=Flatten[ImportString[fulltext[[3]],"TSV"]];
	ReactNames=ConstantArray["NoName",Length[Reactlist]];
	Complist=Flatten[ImportString[fulltext[[4]],"TSV"]];
	Complist=Map[StringReplace[#, {(StartOfString ~~Whitespace) | (Whitespace ~~ EndOfString) :>  "",Whitespace->" "}]&,Complist];
	StoichioExternals=Flatten[ImportString[fulltext[[5]],"TSV"]];
	StoichioExternals=Map[StringReplace[#, {(StartOfString ~~Whitespace) | (Whitespace ~~ EndOfString) :>  "",Whitespace->" "}]&,StoichioExternals];
	Complist=Join[Complist,StoichioExternals];

(*Splitlist=Transpose@LastWord[Complist];*)
	Splitlist=Transpose@StringSplit[Complist,Whitespace,2];
	{CompIDs,Compartlist}=If[Length[Splitlist]== 1,{Splitlist[[1]],ConstantArray[DefaultCompartment,Length[Splitlist[[1]]]]},Splitlist];
	CompartIDs=DeleteDuplicates[Compartlist];
	CompartOuters=ConstantArray["Missing",Length[CompartIDs]];
	CompartSizes=ConstantArray["1.0",Length[CompartIDs]];
	CompartUnits=ConstantArray["volume",Length[CompartIDs]];
	CompartConstant=ConstantArray["Missing",Length[CompartIDs]];
	CompartDims=ConstantArray["3",Length[CompartIDs]];
	If[Length[CompartIDs]<Length[Compartlist]/2,CompNames=ConstantArray["NoName",Length[CompIDs]],
		CompNames=Compartlist;
		Compartlist=ConstantArray[DefaultCompartment,Length[CompIDs]];
		CompartIDs={DefaultCompartment};CompartOuters={"Missing"};];

	S=ImportString[fulltext[[6]],"MTX"];
	Remove[fulltext];True
	]

ExtractCoefficients[equationlist_, leftside_, reactlist_, 
  complist_] :=
 (* equationlist is a list of all left(right)sides of the reaction
 	equations, each as a string in the form A + 2 B + 5 C, where
 	omitted stoichiometry coefficients are assumed as 1
 It returns a set of rules for entries in a sparse matrix  *)
 Module[{pairs, extractcoef, multiply, rules = {}, eqsep}, 
  eqsep = Whitespace ~~ "+" ~~ Whitespace; 
  multiply = If[leftside, -1, 1]; 
  Do[pairs = StringTrim@Flatten@StringSplit[equationlist[[reaction]], eqsep]; 
   extractcoef = Flatten[StringCases[pairs, 
   	{coef : NumberString ~~ Whitespace ~~ comp__ :> {String2Number[coef], comp}, 
   	 Whitespace ... ~~ comp__ -> {1, comp}}]
   	 					, 1];
   extractcoef[[All, 2]] = StringReplace[extractcoef[[All, 2]], Whitespace -> "_"];
  (* Map[If[Position[complist,#]=={},Print[#]]&,extractcoef[[All,2]]];*)
   rules = Join[rules, Cases[extractcoef, 
   	{coef_, comp_} :> {Position[complist, comp] /. List -> Sequence, reaction} -> multiply*coef]];,
  {reaction, Length[reactlist]}];
  rules]

BIOPTreader[filepath_] := 
 Module[{fulltext,percat,blockpairs,deletelines, DividersAt, Headers, Known, Recognized, 
   Unrecognized, removesections, removals, eqsep, eqtoken, eqrev, 
   rxnstart, constart, metabstart, exstart, maxstart, minstart, flunstart,
   obstart, dobstart, end, reactionpattern, Lefts, Rights, fluxval,fluxread, BoundFluxes,rightpairs, leftpairs, 
   ExIDs, exat, oblist, objtoo, prefixterms, signs, terms, splitterms,
   tofix, tofixat, unique,duplicates,excess,absent}, 
   (* Function to calculate the end of a file section *)
  finish[part_] := If[part == 0, 0, Min@Cases[DividersAt, _?(# > part &)] - 1];
  (* Define recognised equation tokens *)
  eqsep = Whitespace ~~ "+" | "-" ~~ Whitespace;
  eqtoken = "->" | "<->" | "-->" | "=>" | "<=>" | "<==>";
  eqrev = "<->" | "<==>";
  FuncDefs = {}; UnitDefs = {}; CompartTypes = {}; SpecieTypes = {};
  Parameters = {}; InitAss = {}; Rules = {}; Constraints = {}; 
  Events = {};
  ReactionXML = {}; SpeciesXML = {};
  
  fulltext = Import[filepath, {"Text", "Lines"}];
  ModelNote = StringJoin[TakeWhile[fulltext, StringMatchQ[#, StartOfLine ~~ "#" ~~ __] &]];
  ModelID = "BioOpt"; ModelName = "No Name"; 
  (* Remove % blocks and lines that start with # from further processing *)
  percat = Quiet@Flatten@Position[fulltext,x_ /; StringMatchQ[x,___~~"%"~~___]];
  If[OddQ[Length[percat]], 
 	DialogInput[{"Data error - An odd number of % signs found, so blocks cannot be defined for deletion!",
 		 DefaultButton[DialogReturn[False]]}], 
 	blockpairs = Partition[percat, 2]; 
 	deletelines = Flatten@Map[Range[#/. List->Sequence]&, blockpairs]; 
 	fulltext = Delete[fulltext, Partition[deletelines, 1]]];
  fulltext = DeleteCases[fulltext, _?(StringMatchQ[#, StartOfLine ~~ "#" ~~ __] &)];
  fulltext=StringReplace[fulltext, "\",\"" -> ","]; (* Fix for EXCEL that puts quotes around commas *)
  (* Identify file sections *)
  Headers = Cases[fulltext, _String?(StringMatchQ[#, "-" ~~ (WordCharacter | Whitespace) ..] &)]; 
  DividersAt = Flatten@Position[fulltext, _String?(StringMatchQ[#,"-" ~~ (WordCharacter | Whitespace) ..] &)]; 
  AppendTo[DividersAt, Length[fulltext] + 1];

  (*CAUTION: If recognition patterns are added under "Known", 
  	corresponding start variables have to be added as indicated below *)
  Known = {
    _String?(StringMatchQ[#, "-REACTION" ~~ ___, IgnoreCase -> True] &),
    _String?(StringMatchQ[#, "-CONSTRAIN" ~~ ___, IgnoreCase -> True] &),
    _String?(StringMatchQ[#, "-METABOLITE" ~~ ___, IgnoreCase -> True] &),
    _String?(StringMatchQ[#, "-EXTERNAL" | "-EXTRACELLULAR" ~~ ___, IgnoreCase -> True] &),
    _String?(StringMatchQ[#, "-MAX" ~~ ___, IgnoreCase -> True] &),
    _String?(StringMatchQ[#, "-MIN" ~~ ___, IgnoreCase -> True] &),
    _String?(StringMatchQ[#, "-FLUXUNIT" ~~ ___, IgnoreCase -> True] &),
    _String?(StringMatchQ[#, "-OBJ" ~~ ___, IgnoreCase -> True] &),
    _String?(StringMatchQ[#, "-DESIGNOBJ" ~~ ___, IgnoreCase -> True] &)}; 

  Recognized = Flatten@Map[Cases[fulltext, #] &, Known]; 
  If[Length[Recognized]==0,EmitSound[Crisis];
  	DialogInput[{"No sections recognised - \nThis is probably not a file in BioOpt format !",
  		DefaultButton[DialogReturn[]]}];Return[False]];
  duplicates = Cases[Tally[Recognized], {dup_, count_ /; count > 1} -> dup];
  If[Length[duplicates] > 0, EmitSound[Crisis]; 
  	DialogInput[{"Data error - The following section headings were duplicated:\n" <>
  		ToString[duplicates] <>  " \nLoading of data file is aborted.", 
    	DefaultButton[DialogReturn[]]}]; Return[False]];
  Unrecognized = Complement[Headers, Recognized]; 
  If[Unrecognized != {}, DialogInput[{
    "The following sections in the data file\n are not recognized and will be ignored:\n" <>
    	 ToString[Unrecognized],
   		DefaultButton[DialogReturn[]]}]]; 
  removesections = Flatten[Map[Position[Headers, #] &, Unrecognized]];
  removals = Flatten/@Map[{{DividersAt[[#]]}, {finish[DividersAt[[#]]]}} &, removesections];
  fulltext = Fold[Drop, fulltext, removals]; 
  DividersAt = Map[Position[fulltext, #] &, Known] /. {} -> {{-1}}; 
  DividersAt = Append[DividersAt[[All, 1, 1]], Length[fulltext] + 1];

  (* Start positions for sections - one variable for each recognition
  	section above *)
  {rxnstart, constart, metabstart, exstart, maxstart, minstart,flunstart, 
    obstart, dobstart, end} = DividersAt + 1;
  If[rxnstart==0,EmitSound[Crisis];
  	DialogInput[{"The required section -REACTIONS not recognised - \nThis is not a valid BioOpt file!",
  		DefaultButton[DialogReturn[]]}];Return[False]];

  (* Extract reaction IDs, names and fluxes from the reaction section, as well as LHS and RHS *)
 (* reactionpattern=Shortest[id__]~~":"~~(name__~~":"~~Whitespace)...~~lhs___~~Whitespace~~eqtoken~~rhs___; *)
  reactionpattern = Shortest[id__] ~~ ":" ~~ (Shortest[name__] ~~ ":") ... ~~ 
  					(lhs___ /;StringFreeQ[lhs, " :"]) ~~ Whitespace ~~ eqtoken ~~ rhs___;
  {Reactlist, ReactNames, Lefts} = Transpose@Flatten[
  	StringCases[fulltext[[rxnstart ;; finish[rxnstart]]], 
    	reactionpattern -> {StringTrim[id], StringTrim[name], lhs}]/. "" -> " ", 1]; 
  ReactNames = Replace[ReactNames, _?(StringFreeQ[#, __] &) -> "NoName", 1];
  Reverselist = Flatten@StringCases[fulltext[[rxnstart ;; finish[rxnstart]]], 
  		Shortest[react__]~~":"~~___~~Whitespace~~eqrev~~___	:> StringTrim[react]];
  RHS = StringSplit[Flatten[StringCases[fulltext[[rxnstart ;; finish[rxnstart]]], 
      reactionpattern -> rhs] /. "" -> " ", 1], Whitespace ~~ ":"];
  {Rights, fluxval} = Transpose[RHS /.
  	 {{x_?(NumberQ@String2Number[#] &)} -> {" ",x}, {x_?(! NumberQ@String2Number[#] &)} -> {x," "}}];
  FluxValues = MapThread[{#1, String2Number[#2]} &, {Reactlist, fluxval}];
  FluxValues = DeleteCases[FluxValues, {__, __?(!NumericQ[#]&)}];
  (* Split LHS and RHS into terms, each consisting of coefficient and name *)
  leftpairs = StringTrim@Flatten@StringSplit[Lefts, eqsep];
  rightpairs = StringTrim@Flatten@StringSplit[Rights, eqsep];
  (* First make a list of all IDs found in reaction equations *)
  CompIDs = Union[
  	StringTrim[leftpairs, NumberString ~~ Whitespace], 
    StringTrim[rightpairs, NumberString ~~ Whitespace]];
  CompIDs=DeleteCases[CompIDs, ""];
  CompNames = ConstantArray["NoName", Length[CompIDs]];
  Compartlist =ConstantArray[DefaultCompartment, Length[CompIDs]];
  
  (* If there is a metabolites section, get an extended metabolite list from there *)
  If[metabstart > 0, {Complist, CompNames, Compartlist} = Transpose@Flatten[
  	StringCases[fulltext[[metabstart ;; finish[metabstart]]], 
      Shortest[id__]~~":"~~Whitespace~~(name__~~":"~~Whitespace)...~~(compart__)...
      	-> {StringTrim[id], StringTrim[name], StringTrim[compart]}], 1];
  absent = Complement[CompIDs, Complist]; 
  Complist = StringReplace[Complist, Whitespace -> "_"];
  CompNames = Replace[CompNames, _?(StringFreeQ[#,__]&)->"NoName", 1];
  Compartlist = Replace[Compartlist, _?(StringFreeQ[#,__]&)->DefaultCompartment, 1];
  If[absent != {}, 
   DialogInput[{"Warning - metabolites\n" <> ToString[absent] <> 
      "\nappear in reaction equations but are not listed in the METABOLITES section.
       \nDefault names and compartments are used for them.", DefaultButton[DialogReturn[]]}]
   ];
  (* Remove blanks from metabolite IDs to avoid trouble with composite naming used in Netspitter *)
  CompIDs = Join[Complist, StringReplace[absent, Whitespace -> "_"]]; 
  CompNames = Join[CompNames, ConstantArray["NoName", Length[absent]]];
  Compartlist = Join[Compartlist, ConstantArray[DefaultCompartment, Length[absent]]];
  ,CompIDs = StringReplace[CompIDs, Whitespace -> "_"];
  ];
 
 (* Now extract the coefficients and slot them into the appropriate row and column *)
  S = SparseArray[Join[
  	ExtractCoefficients[Lefts, True, Reactlist, CompIDs], 
    ExtractCoefficients[Rights, False, Reactlist, CompIDs]]];

  (* Duplication between metabolite IDs and reaction IDs not allowed by SBML, 
  	so fix this by adding a "-RXN" suffix to relevant reaction IDs*)
  tofix = Intersection[Reactlist, CompIDs];
  If[tofix != {}, 
  	tofixat = Flatten@Map[Position[Reactlist, #] &, tofix];
  	Reactlist[[tofixat]] = Map[# <> "_RXN" &, Reactlist[[tofixat]]];
  	tofixat = Flatten@Map[Position[Reverselist, #] &, tofix];
  	Reverselist[[tofixat]] = Map[# <> "_RXN" &, Reverselist[[tofixat]]];
	tofixat = Flatten@Map[Position[FluxValues[[All,1]], #] &, tofix];
	FluxValues[[tofixat,1]] = Map[# <> "_RXN" &, FluxValues[[tofixat,1]]];
  	];
 
  Complist = MapThread[StringJoin, {CompIDs, ConstantArray[" ", Length[CompIDs]],Compartlist}];
  CompartIDs=DeleteDuplicates[Compartlist];
  CompartOuters = ConstantArray["Missing", Length[CompartIDs]];
  CompartSizes = ConstantArray["1.0", Length[CompartIDs]];
  CompartUnits = ConstantArray["volume", Length[CompartIDs]];
  CompartConstant = ConstantArray["Missing", Length[CompartIDs]];
  CompartDims = ConstantArray["3", Length[CompartIDs]];
  
  (* Read flux unit specification if given in the file*)
  If[flunstart > 0, 
	fluxread = StringCases[fulltext[[flunstart]], 
		unit__ ~~ ":" ~~ Shortest[flot__] ~~ Whitespace ~~ Shortest[flos__] ~~ Whitespace ~~ 
		Shortest[DWt__] ~~ Whitespace ~~ Shortest[DWs__] ~~ Whitespace ~~ timef__ :> 
		{StringTrim[unit], StringTrim[flot], flos,DWt, DWs,  StringTrim[timef]}];
	If[fluxread != {},
		{{FluxUnits, flotype, floscale, DWtype, DWscale, timefactor}} = fluxread]];
      
  (* Read constraints, checking for consistency with reaction listing *)
  FluxBounds = If[constart > 0, Flatten[StringCases[fulltext[[constart ;; finish[constart]]], 
      react__~~Whitespace...~~"["~~lo__~~","~~hi__~~"]" :>
      	{StringTrim[react], String2Number[lo], String2Number[hi]}], 1], {{}}];
  unique = DeleteDuplicates[FluxBounds, SameQ[#1[[1]], #2[[1]]] &];
  excess = Length[FluxBounds] - Length[unique]; 
  If[excess > 0, excess=Flatten@Cases[Tally[FluxBounds], {multi_, _?(# > 1 &)} -> multi];DialogInput[{
  		"Data error - constraint specifications\n" <> ToString[excess] <> 
   		"\nare duplicated in the input file.
   		\nOnly the first one is retained in all cases.",
   		DefaultButton[DialogReturn[]]}]];
  FluxBounds = unique;
  If[tofix != {},   
  	tofixat = Flatten@Map[Position[FluxBounds[[All,1]], #] &, tofix];
	FluxBounds[[tofixat,1]] = Map[# <> "_RXN" &, FluxBounds[[tofixat,1]]]];
	
  absent = Complement[FluxBounds[[All, 1]], Reactlist];
  If[absent != {}, DialogInput[{
   			"Data error - the constraint entries\n" <> ToString[absent] <> 
    		"\ndo not appear on the reaction list and are omitted.",
   			DefaultButton[DialogReturn[]]}]; 
   		Do[FluxBounds = DeleteCases[FluxBounds, {absent[[i]], _, _}], {i, Length[absent]}]];
  FluxSigns = Map[MapAt[Sign, #, 2] &, FluxValues];
  BoundFluxes = If[Length[FluxBounds] > 0, 
	Map[{#[[1]], Sign[Sign[#[[2]]] + Sign[#[[3]]]]} &, FluxBounds], {}];
  BoundFluxes = DeleteCases[BoundFluxes, {_, 0}];
  FluxSigns = Union[FluxSigns, BoundFluxes];
  If[FluxSigns != {}, Fluxmess = "\n"<>ToString[Length[FluxSigns]]<>
   " Fluxes and/or Bounds were read from the BioOpt input file."];
  
  (* Read explicitly listed external compounds, check for typo's *)
  StoichioExternals = If[exstart>0, 
    ExIDs = Flatten@StringCases[fulltext[[exstart ;; finish[exstart]]], comp__ :>
    	StringTrim[comp]]; 
    ExIDs = DeleteCases[ExIDs, ""];ExIDs = StringReplace[ExIDs, Whitespace -> "_"];
    absent = Complement[ExIDs,CompIDs]; 
    If[absent != {}, 
   DialogInput[{"Data error - metabolites\n" <> ToString[absent] <> 
      "\nlisted as external metabolites do not appear in reaction equations\nnor in the METABOLITES section.
       \nThey are ignored.", DefaultButton[DialogReturn[]]}];ExIDs=Intersection[CompIDs,ExIDs];
   ];
    exat = Map[Position[CompIDs, #] &, ExIDs];Complist[[Flatten[exat]]], {}];
  
  (* Read objectives, with and/or without explicit min/max specifications. 
  	Secondary objective is also read, and a text qualifier identified to distinguish it from
  	the main objective. *)
  oblist = {}; ObjWeights = {}; 
  If[maxstart > 0, AppendTo[oblist, 
    Flatten@StringCases[fulltext[[maxstart ;; finish[maxstart]]], 
      rhs__ :> {"MAXIMIZE_", StringTrim[rhs]}]]]; 
  If[minstart > 0, AppendTo[oblist, 
    Flatten@StringCases[fulltext[[minstart ;; finish[minstart]]], 
      rhs__ :> "MINIMIZE_" <> StringTrim[rhs]]]];
  If[obstart > 0, oblist = Join[oblist, 
     Flatten[StringCases[fulltext[[obstart ;; finish[obstart]]],
     	{"min"|"MIN" ~~ rhs__ :> {"MINIMIZE_", StringTrim[rhs]}, 
         "max"|"MAX" ~~ rhs__ :> {"MAXIMIZE_", StringTrim[rhs]}, 
         rhs__ :> {" ", StringTrim[rhs]}
        }], 1]]];
  If[dobstart > 0, objtoo = 
    StringCases[fulltext[[dobstart - 1]],
    	"-"~~qual:WordCharacter..~~"OBJ"~~___ :> qual][[1]];
   oblist = Join[oblist, Flatten[StringCases[fulltext[[dobstart ;; finish[dobstart]]],
   	{"min"|"MIN" ~~ rhs__ :> {"MINIMIZE_"~~objtoo, StringTrim[rhs]}, 
     "max"|"MAX" ~~ rhs__ :> {"MAXIMIZE_"~~objtoo, StringTrim[rhs]}, 
      rhs__ :> {objtoo, StringTrim[rhs]}
    }], 1]]];
    oblist=DeleteCases[oblist, {_, ""}];
    (* Split the objective expressions into terms using arithmetic signs *)
  If[oblist != {}, prefixterms = Map[StringSplit[#, s:eqsep :> StringTrim[s]]&, oblist]; 
   prefixterms[[All, 2]] = Map[If[First[#]=="+"||First[#]=="-",Null, Prepend[#, "+"]] &, 
     prefixterms[[All, 2]]];
    (* Now isolate coefficients and make them negative where appropriate *)
   Do[ 
   	signs = Cases[prefixterms[[ob, 2]], "+" | "-"] /. {"+" -> 1, "-" -> -1};
    terms = Cases[prefixterms[[ob, 2]], Except["+" | "-"]]; 
    splitterms = Flatten[StringCases[terms,
    	{coef:NumberString ~~ Whitespace ~~ comp__ ->	{comp, coef},
    	Whitespace ... ~~ comp__ -> {comp, "1"}
    	}], 1]; 
    splitterms[[All, 2]] = signs*Map[String2Number[#] &, splitterms[[All, 2]]]; 
    splitterms = Map[Insert[#1, StringTrim@prefixterms[[ob,1,1]], 2]&, splitterms];
    ObjWeights = Join[ObjWeights, splitterms], {ob, Length[prefixterms]}]];
  (* Where reaction IDs were made unique before, this has to be applied to objectives too *)  
  If[tofix != {},   
	tofixat = Flatten@Map[Position[ObjWeights[[All,1]], #] &, tofix];
	ObjWeights[[tofixat,1]] = Map[# <> "_RXN" &, ObjWeights[[tofixat,1]]]];
  (* Check that objectives refer to known reactions *)
  absent = Complement[ObjWeights[[All, 1]], Reactlist];
  If[absent != {}, DialogInput[{
		"Data error - the objective entries\n" <> ToString[absent] <> 
   		"\ndo not appear on the reaction list and are omitted.",
		DefaultButton[DialogReturn[]]}]; 
   		Do[ObjWeights = DeleteCases[ObjWeights, {absent[[i]], _, _}], {i, Length[absent]}]];
  
  Remove[fulltext];
   True ]

TabsepnBlock[Namelist_]:=StringJoin[Riffle[StringJoin[Riffle[#,"\t"]]&/@Partition[Namelist,5,5,1,{}],"\n"]]

Reconcile[AllItems_,WithAttributes_,Default_:None]:=Module[{flag=1,attribute},
	Map[(attribute=Default;
		Which[flag>Length[WithAttributes],
			Null,MatchQ[WithAttributes[[flag]],{#,_}],
			attribute=WithAttributes[[flag,2]];flag++;Exit;,
			True,Null];attribute)&,AllItems]
	]

TSVExport[Title_,Blocklisting_]:=Module[{BlockIntRows,BlockIntCols,BlockRevCols,BlockRows,BlockExtRows,
	BlockStoi,Inflows,Outflows,Crossflows,Reactlisting,Reverselisting,Intlisting,Extlisting,
	Smatlisting,Outfile,i},
If[Length[Blocklisting]>0,
	BlockIntRows=Flatten[Table[Position[Complist,Blocklisting[[i]]],{i,1,Length[Blocklisting]}]];
	BlockIntCols=Nonzerows[Normal[Transpose[S[[BlockIntRows]]]]];
	BlockRevCols=BlockIntCols\[Intersection]Reversibles;
	BlockRows=Nonzerows[Normal[S[[All,BlockIntCols]]]];
	BlockExtRows=Complement[BlockRows,BlockIntRows];
	Inflows=BlockExtRows\[Intersection]InputExternalrows;
	Outflows=BlockExtRows\[Intersection]OutputExternalrows;
	Crossflows=Complement[BlockExtRows,Inflows\[Union]Outflows];
	BlockExtRows=Join[Inflows,Crossflows,Outflows];
	BlockRows=Join[BlockIntRows,BlockExtRows];
	BlockStoi=S[[BlockRows,BlockIntCols]];
	Reactlisting=TabsepnBlock[Reactlist[[BlockIntCols]]];
	Reverselisting=TabsepnBlock[Reactlist[[BlockRevCols]]];
	Intlisting=TabsepnBlock[Complist[[BlockIntRows]]];
	Extlisting=TabsepnBlock[Complist[[BlockExtRows]]];
	Smatlisting=ExportString[SparseArray[BlockStoi],"MTX"]; 
	Outfile=Title[[1]]<>".tsv";
	Export[Outfile,{"<Title>",Title[[1]]<>" : "<>Title[[2]],"<Reactions>",Reactlisting,"<ReversibleReactions>",Reverselisting,
		"<InternalCompounds>",Intlisting,"<ExternalCompounds>",Extlisting,
		"<Stoichiometry>",Smatlisting},"List"];
]]

SBMLreader[filepath_]:=Module[{InputXML,NoteXML,ExIDs,ExNames,Srules,SBMLerror,
	FuncDefXML,UnitDefXML,CompartTypeXML,SpecieTypeXML,ParameterXML,
	InitAssXML,RulesXML,ConstraintsXML,EventsXML,LowLims,UpLims,BoundFluxes},
	
	InputXML=Import[filepath,"XML"];
	{{ModelID, ModelName}} = Cases[InputXML,XMLElement["model", {___,"id" -> modid_, ___,"name" -> modname_, ___}, _] :> {modid, modname}, Infinity];
	NoteXML = Cases[InputXML, XMLElement["model",_,{XMLElement["notes",{},{note_}],__}]:>note, Infinity];
	NoteXML = If[NoteXML == {}, "No notes read from SBML", NoteXML[[1]]];
	ModelNote = If[StringQ[NoteXML], NoteXML, StringJoin@Riffle[Flatten@Cases[NoteXML, XMLElement["p", _, parag_] :> parag, Infinity],"\n"]];
(* Note: global definitions appear at nesting level 5 in InputXML. 
   Using Infinity will pick up local Parameter declarations in reactions (level 11) as well.*)
	FuncDefXML = Cases[InputXML, XMLElement["listOfFunctionDefinitions", _, _], 5];
	FuncDefs = XMLtoText[FuncDefXML];
	UnitDefXML = Cases[InputXML, XMLElement["listOfUnitDefinitions", _, _], 5];
	UnitDefs = XMLtoText[UnitDefXML];
	CompartTypeXML = Cases[InputXML, XMLElement["listOfCompartmentTypes", _, _], 5];
	CompartTypes = XMLtoText[CompartTypeXML];
	SpecieTypeXML = Cases[InputXML, XMLElement["listOfSpeciesTypes", _, _], 5];
	SpecieTypes = XMLtoText[SpecieTypeXML];
	ParameterXML = Cases[InputXML, XMLElement["listOfParameters", _, _], 5];
	Parameters = XMLtoText[ParameterXML];
	InitAssXML = Cases[InputXML, XMLElement["listOfInitialAssignments", _, _], 5];
	InitAss = XMLtoText[InitAssXML];
	RulesXML = Cases[InputXML, XMLElement["listOfRules", _, _], 5];
	Rules = XMLtoText[RulesXML];
	ConstraintsXML = Cases[InputXML, XMLElement["listOfConstraints", _, _], 5];
	Constraints = XMLtoText[ConstraintsXML];
	EventsXML = Cases[InputXML, XMLElement["listOfEvents", _, _], 5];
	Events = XMLtoText[EventsXML];

	CompartIDs=Cases[InputXML,XMLElement["compartment",{___,"id"->compartid_,___},{}]:> compartid,Infinity];
	CompartOuters=Cases[InputXML,XMLElement["compartment",{___,"id"->compartid_,___,"outside"->outid_,___},{}]:> {compartid,outid},Infinity];
	CompartOuters=Reconcile[CompartIDs,CompartOuters,"Missing"];
	CompartDims=Cases[InputXML,XMLElement["compartment",{___,"id"->compartid_,___,"spatialDimensions"->outid_,___},{}]:> {compartid,outid},Infinity];
	CompartDims=Reconcile[CompartIDs,CompartDims,"3"];
	CompartSizes=Cases[InputXML,XMLElement["compartment",{___,"id"->compartid_,___,"size"->outid_,___},{}]:> {compartid,outid},Infinity];
	CompartSizes=Reconcile[CompartIDs,CompartSizes,"Missing"];
	CompartUnits=Cases[InputXML,XMLElement["compartment",{___,"id"->compartid_,___,"units"->outid_,___},{}]:> {compartid,outid},Infinity];
	CompartUnits=Reconcile[CompartIDs,CompartUnits,"Missing"];
	CompartConstant=Cases[InputXML,XMLElement["compartment",{___,"id"->compartid_,___,"constant"->outid_,___},{}]:> {compartid,outid},Infinity];
	CompartConstant=Reconcile[CompartIDs,CompartConstant,"Missing"];

	ReactionXML = Cases[InputXML, XMLElement["reaction", _, _], Infinity];
	Reactlist=Cases[ReactionXML,XMLElement["reaction",{___,"id"->reactid_,___},{___}]:>reactid,Infinity];
	Reverselist=Complement[Reactlist,Cases[ReactionXML,XMLElement["reaction",{___,"id"->reactid_,___,"reversible"->"false",___},{___}]:>reactid,Infinity]];
	ReactNames=Cases[ReactionXML,XMLElement["reaction",{___,"id"->reactid_,___,"name"->reactname_,___},{}]:> {reactid,reactname},Infinity];
	ReactNames=Reconcile[Reactlist,ReactNames,"NoName"];

    FluxValues = Cases[ReactionXML,
		XMLElement["reaction", {___, "id" -> reactid_, ___}, {__, 
     		XMLElement["kineticLaw", {}, {___, 
				XMLElement["listOfParameters", {}, {___, 
					XMLElement["parameter", {___,"id"->"FLUX"|"FLUX_VALUE"|"FLUXVALUE",___,"value"->fluxval_,___},{}],
		 ___}],___}],___}] :> {reactid,String2Number[fluxval]}, 1];
	LowLims = Cases[ReactionXML,
		XMLElement["reaction",{___,"id"->reactid_,___}, {__, 
     		XMLElement["kineticLaw",{},{___, 
				XMLElement["listOfParameters",{},{___, 
					XMLElement["parameter",{___,"id"->"LOWERBOUND"|"LOWER_BOUND",___,"value"->fluxval_,___},{}],
		___}],___}],___}] :> {reactid,String2Number[fluxval]},1];
	UpLims = Cases[ReactionXML,
		XMLElement["reaction",{___,"id"->reactid_,___}, {__, 
     		XMLElement["kineticLaw",{},{___, 
				XMLElement["listOfParameters",{},{___, 
					XMLElement["parameter",{___,"id"->"UPPERBOUND"|"UPPER_BOUND",___,"value"->fluxval_,___},{}],
		___}],___}],___}] :> {reactid,String2Number[fluxval]},1];
	FluxSigns = Map[MapAt[Sign, #, 2] &, DeleteCases[FluxValues, {_, 0}]];
	BoundFluxes = If[Length[UpLims] == Length[LowLims], 
		MapThread[{#1[[1]], Sign[Sign[#1[[2]]] + Sign[#2[[2]]]]} &, {UpLims,LowLims}], {}];
	BoundFluxes = DeleteCases[BoundFluxes, {_, 0}];
	FluxSigns = Union[FluxSigns, BoundFluxes];
	If[FluxSigns != {}, Fluxmess = "\n"<>ToString[Length[FluxSigns]]<>
		" Fluxes and/or Bounds were read from the SBML input file."];
	
	SpeciesXML = Cases[InputXML, XMLElement["species", _, _], Infinity];
	CompIDs=Cases[SpeciesXML,XMLElement["species",{___,"id"->compid_,___},{}]:> compid,Infinity];
	CompNames = Cases[SpeciesXML, XMLElement["species", {___,"id"->compid_,___, "name" -> compname_, ___}, {}] :> {compid,compname}, Infinity];
	CompNames=Reconcile[CompIDs,CompNames,"NoName"];
	Compartlist = Cases[SpeciesXML, XMLElement["species", {___,"id"->compid_,___, "compartment" -> compname_, ___}, {}] :> {compid,compname}, Infinity];
	If[Compartlist!={},Compartlist=Reconcile[CompIDs,Compartlist,DefaultCompartment]];

	If[CompIDs[[1]] == "DUMMY",
		SpeciesXML = Rest[SpeciesXML]; CompIDs = Rest[CompIDs]; CompNames = Rest[CompNames]; 
  		If[Compartlist != {}, Compartlist = Rest[Compartlist];]];
  	
	SBMLerror=Catch[Srules=Flatten[Map[ReactionParse,ReactionXML],1];];
	If[SBMLerror,Return[False]];
	S=SparseArray[Srules];
		
	ExIDs=Cases[InputXML,XMLElement["species",{___,"id"->compid_,___,"boundaryCondition"->"true",___},{}]:> compid,Infinity];
	ExNames = Cases[InputXML, 
		XMLElement["species",{___,"id"->compid_,___,"name"->compname_,___,"boundaryCondition"->"true",___},{}] | 
     	XMLElement["species",{___,"id"->compid_,___,"boundaryCondition"->"true",___,"name"->compname_,___},{}] :>{compid, compname}, Infinity];
	ExNames = Reconcile[ExIDs, ExNames, "NoName"];
	If[CompIDs!= CompNames,
		StoichioExternals=MapThread[StringJoin,{ExIDs,ConstantArray[" ",Length[ExIDs]],ExNames}];
		Complist=MapThread[StringJoin,{CompIDs,ConstantArray[" ",Length[CompIDs]],CompNames}],
		StoichioExternals=ExIDs; Complist=CompIDs];
	Complist=Map[StringReplace[#, {(StartOfString ~~Whitespace) | (Whitespace ~~ EndOfString) :>  "",Whitespace->" "}]&,Complist];
	StoichioExternals=Map[StringReplace[#, {(StartOfString ~~Whitespace) | (Whitespace ~~ EndOfString) :>  "",Whitespace->" "}]&,StoichioExternals];
Remove[InputXML];True]


XMLtoText[XML_] := 
 If[XML == {}, {}, 
  Flatten@Map[
    StringSplit[
      StringReplace[
       ExportString[#, "XML"], {"\r\n" -> "", 
        ">" -> ">\[DownArrow]"}], "\[DownArrow]"] &, XML]]
     
ReactionParse[RxnSBML_]:=Module[{Scol,SrowIn,SrowOut,reaction},
	Scol=RxnSBML/.XMLElement["reaction",{___,"id"->rname_,___},_]:>Position[Reactlist,rname,1][[1,1]];
	reaction=Reactlist[[Scol]];
	SrowIn=Cases[Cases[RxnSBML,XMLElement["listOfReactants",_,_],{2}],
		XMLElement["speciesReference", {___,"species" -> spec_,___,"stoichiometry" -> stoich_, ___}, {}] | 
 		XMLElement["speciesReference", {___, "species" -> spec_, ___}, {}]:>
			{Quiet@Check[Position[CompIDs,spec,1][[1,1]],Squeal[spec,reaction];Throw[True],Part::partw],
				N@If[stoich== "",1,String2Number[stoich]]},{3}];
	SrowOut=Cases[Cases[RxnSBML,XMLElement["listOfProducts",_,_],{2}],
		XMLElement["speciesReference", {___,"species" -> spec_,___,"stoichiometry" -> stoich_, ___}, {}] | 
 		XMLElement["speciesReference", {___, "species" -> spec_, ___}, {}]:>
			 {Quiet@Check[Position[CompIDs,spec,1][[1,1]],Squeal[spec,reaction];Throw[True],Part::partw],N@If[stoich== "",1,String2Number[stoich]]},{3}];
	Join[Map[{#[[1]],Scol}->-#[[2]]&,SrowIn],Map[{#[[1]],Scol}->#[[2]]&,SrowOut]]]

Squeal[spec_,rname_]:=DialogInput[DialogNotebook[{TextCell[Style[
		"Inconsistent data - reading abandoned.\n\nMetabolite "<>spec<>"\nin reaction "<>rname<>
		"\nnot found in SBML listofSpecies section ! ",Purple,"Panel"]],
		Button["OK",DialogReturn[1]]},WindowTitle->"SBML Error"]];

ReactantSBML[i_]:=Module[{speciecount},
	speciecount=Length[ Reactions[[i,2]]];
	If[speciecount>0,Join[
		{"<listOfReactants>"},
		Table["<speciesReference species=\"" ~~Reactions[[i,2,j,1]]~~ "\""~~ 
			" stoichiometry=\""~~ ToString[Reactions[[i,2,j,2]]]~~ "\" />",{j,1,speciecount}],
		{"</listOfReactants>"}],
		{}]];

ProductSBML[i_]:=Module[{speciecount},
	speciecount=Length[ Reactions[[i,3]]];
	If[speciecount>0,Join[
		{"<listOfProducts>"},
		Table["<speciesReference species=\"" ~~Reactions[[i,3,j,1]]~~ "\""~~ 
			" stoichiometry=\""~~ ToString[Reactions[[i,3,j,2]]]~~ "\" />",{j,1,speciecount}],
		{"</listOfProducts>"}],
		{}]];

KineticParamSBML[reactid_,FluxValues_,FluxBounds_,ObjWeights_] := 
  Module[{flux,fluxval,bounds, low, high, objective = True, weight, defs, headsbml, 
		boundsbml, objsbml,fluxsbml, prefix,prefsplit,mul=1, tailsbml}, 
  fluxval = Cases[FluxValues, {reactid, val_?NumberQ} -> val];
  If[Length[fluxval]==0,fluxval=0;flux=False,fluxval=fluxval[[1]];flux=True];
  bounds = Cases[FluxBounds, {reactid, lo_?NumberQ, hi_?NumberQ} -> {lo, hi}];
  If[ArrayQ[bounds, 2], {{low, high}} = bounds; bounds = True, 
  						{{low, high}} = {{0, 1000}}; bounds = False];
  objective = Cases[ObjWeights, {reactid, pref__, w_?NumberQ} -> {pref, w}]; 
  objective = If[Length[objective] == 0, False, {prefix, weight}=Transpose[objective]; True];
  If[!bounds && !objective && !flux, Return[{}]];
  defs={};
  If[flux,defs={"<ci>FLUX_VALUE</ci>"}];
  If[bounds,defs=Join[defs,{"<ci>LOWER_BOUND</ci>","<ci>UPPER_BOUND</ci>"}]];
  If[objective,defs=Join[defs,Table[
  	
  	(* These 3 lines assume that the optimizer minimizes, so it reverses objective weights for maximization*)
  	If[prefix[[ob]]!="",prefsplit = StringSplit[prefix[[ob]], "_"]; 
  		mul = If[prefsplit[[1]] == "MAXIMIZE", -1, 1];
		prefix[[ob]] = If[Length[prefsplit] > 1, prefsplit[[2]], ""]];
	(* If they are commented out, the MAXIMIZE or MINIMIZE part of the prefix is put in the SBML   *)
	
	"<ci>" <> prefix[[ob]] <> 
      "OBJECTIVE_COEFFICIENT</ci>",{ob, Length[prefix]}]]];
  headsbml = Join[{"<kineticLaw>", 
    		  "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">"},defs, 
    		   {"</math>", "<listOfParameters>"}];
  boundsbml = If[bounds, 
  	{"<parameter id=\"LOWER_BOUND\" value=\"" <> ToString[low] <> 
      "\" units=\""<>FluxUnits<>"\"/>",
     "<parameter id=\"UPPER_BOUND\" value=\"" <> ToString[high] <> 
      "\" units=\""<>FluxUnits<>"\"/>"}, {}];
  objsbml = If[objective, Table[ "<parameter id=\"" <> prefix[[ob]] <> 
      "OBJECTIVE_COEFFICIENT\" value=\"" <> ToString[mul*weight[[ob]]] <> "\" units=\"dimensionless\"/>",
       						{ob, Length[weight]}], {}]; 
  fluxsbml=If[flux,
  	{"<parameter id=\"FLUX_VALUE\" value=\""<>ToString[fluxval]<>"\" units=\""<>FluxUnits<>"\"/>"},{}];
  tailsbml = {"</listOfParameters>", "</kineticLaw>"};
  Join[headsbml, boundsbml, objsbml, fluxsbml, tailsbml]
  ]

PutAttribute[xml_, element_, attribute_, ValueSet_] := 
 Module[{XML = DeleteCases[xml, (attribute -> _), 3]}, 
  Do[XML[[i]]=XML[[i]] /. 
     XMLElement[element, {tags___}, specs___] -> 
      XMLElement[element, {tags, Rule[attribute, ValueSet[[i]]]}, specs],
      	{i, Length[XML]}]; XML]

SBMLfixID[compids_] := Module[{fixed, duplis},
  (*SBML rules only allow letter,digit or understroke,no digit at start of ID*)
  fixed = StringReplace[compids, 
  	{"-" -> "_", "+" -> "x", Except[WordCharacter] .. -> "_",
      StartOfString ~~ s : Except[LetterCharacter | "_"] -> "_" ~~ s}];
  duplis = Length[fixed] - Length@DeleteDuplicates[fixed];
  If[duplis > 0, 
   DialogInput[DialogNotebook[{TextCell[
       "Forced conformation to SBML ID naming rules \ncaused duplicate ID's:\n" <> 
        ToString@Take[Commonest[fixed], Min[3, duplis]] <> 
       "\nManual editing of ID's required to fix this."], 
      Button["OK", DialogReturn[1]]}]]];
  fixed]
   
SBMLfixName[name_]:=(StringReplace[name, {"&" -> "&amp;", "\"" -> "&quot;","'" -> "&apos;"}])   
   
SBMLExport[Title_,Blocklisting_]:=Module[{BlockIntRows,BlockIntCols,BlockRows,BlockExtRows,
	Inflows,Outflows,Crossflows,BlockStoi,Header,IntNames,ExtNames,
	BlockCompartIDs,BlockCompartOuters,BlockCompartDims,BlockCompartConstant,BlockCompartSize,BlockCompartUnit,
	CompartmentSBL,Speclist,IntIDs, ExtIDs,BlockBC,BlockSpeclist, SpeciesSBML, BlockReactlist,BlockReactNames,
	SBMLReverselist,tofix,tofixat,BlockRev,Stoicol,Reacts,Prods,ReactionSBML,SBMLfluxes,SBMLbounds,SBMLobjectives,Outfile,i},
	Outfile=StringJoin[StringSplit[Title][[1;;2]],".sbml"];
	Outfile=Title[[1]]<>".sbml";
	If[Length[Blocklisting]>0,
		BlockIntRows=Flatten[Table[Position[Complist,Blocklisting[[i]]],{i,1,Length[Blocklisting]}]];
		BlockIntCols=Nonzerows[Normal[Transpose[S[[BlockIntRows]]]]];
		BlockRows=Nonzerows[Normal[S[[All,BlockIntCols]]]];
		BlockExtRows=Complement[BlockRows,BlockIntRows];
		Inflows=BlockExtRows\[Intersection]InputExternalrows;
		Outflows=BlockExtRows\[Intersection]OutputExternalrows;
		Crossflows=Complement[BlockExtRows,Inflows\[Union]Outflows];
		BlockExtRows=Join[Inflows,Crossflows,Outflows];
		BlockRows=Join[BlockIntRows,BlockExtRows];
		BlockStoi=S[[BlockRows,BlockIntCols]];
		Header={"<?xml version=\"1.0\" encoding=\"UTF-8\"?> <sbml ",
			"xmlns=\"http://www.sbml.org/sbml/level2/version4\" version=\"4\" level=\"2\"",
			"xmlns:math=\"http://www.w3.org/1998/Math/MathML\"",
			"xmlns:html=\"http://www.w3.org/1999/xhtml\"   >",
			"<model id=\"" ~~Title[[1]]~~ "\" name=\"" ~~Title[[2]] ~~"\" >",
			"<notes>  <body xmlns=\"http://www.w3.org/1999/xhtml\"> <pre>",
			"Metabolic subnetwork extracted by NetSplitter on "~~DateString[]~~
			".\nStoichiometry file processed is: " ~~Infile ,
			"</pre>  </body>  </notes>"};
		If[UnitDefs=={}, UnitDefs={"<listOfUnitDefinitions>",
			"<unitDefinition id=\""<>FluxUnits<>"\" >","<listOfUnits>",
			"<unit kind=\""<>flotype<>"\" scale=\""<>floscale<>"\" multiplier=\"1\" />",
			"<unit kind=\""<>DWtype<>"\" scale=\""<>DWscale<>"\" exponent=\"-1\" multiplier=\"1\" />",
			"<unit kind=\"second\" multiplier=\""<>timefactor<>"\" exponent=\"-1\" />",
			"</listOfUnits>","</unitDefinition>","</listOfUnitDefinitions>"}];	
		BlockCompartIDs=Union[Compartlist[[BlockRows]]];
		BlockCompartOuters=CompartOuters[[Flatten[Table[Position[CompartIDs,
			BlockCompartIDs[[i]]],{i,1,Length[BlockCompartIDs]}]]]];
		BlockCompartDims=CompartDims[[Flatten[Table[Position[CompartIDs,
			BlockCompartIDs[[i]]],{i,1,Length[BlockCompartIDs]}]]]];
		BlockCompartSize=CompartSizes[[Flatten[Table[Position[CompartIDs,
			BlockCompartIDs[[i]]],{i,1,Length[BlockCompartIDs]}]]]];
		BlockCompartUnit=CompartUnits[[Flatten[Table[Position[CompartIDs,
			BlockCompartIDs[[i]]],{i,1,Length[BlockCompartIDs]}]]]];
		BlockCompartConstant=CompartConstant[[Flatten[Table[Position[CompartIDs,
			BlockCompartIDs[[i]]],{i,1,Length[BlockCompartIDs]}]]]];
		CompartmentSBL=Table["<compartment id=\""~~ BlockCompartIDs[[i]]~~
			"\" spatialDimensions=\""~~BlockCompartDims[[i]]~~"\""~~
			If[BlockCompartOuters[[i]]== "Missing",""," outside=\""~~ BlockCompartOuters[[i]]~~"\""]~~
			If[BlockCompartSize[[i]]== "Missing",""," size=\""~~BlockCompartSize[[i]]~~"\""]~~
			If[BlockCompartUnit[[i]]== "Missing",""," units=\""~~BlockCompartUnit[[i]]~~"\""]~~
			If[BlockCompartConstant[[i]]== "Missing",""," constant=\""~~BlockCompartConstant[[i]]~~"\""]~~
			" />",{i,1,Length[BlockCompartIDs]}];

	  If[SpeciesXML == {},
		IntIDs = SBMLfixID[CompIDs[[BlockIntRows]]]; ExtIDs = SBMLfixID[CompIDs[[BlockExtRows]]];
		IntNames = SBMLfixName[CompNames[[BlockIntRows]]]; ExtNames = SBMLfixName[CompNames[[BlockExtRows]]];
		Speclist = Join[IntIDs, ExtIDs];
(*Dummy is a workaround for CellnetAnalyzer bug, 
	name field in first record not read and upsets internal/external recognition *)
		SpeciesSBML= Join[{"<species id=\"DUMMY\" compartment=\""<>DefaultCompartment<>"\" name=\"Dummy\" />"}, 
			Table["<species id=\""~~IntIDs[[i]]~~
				"\" compartment=\""~~Compartlist[[BlockIntRows[[i]]]]~~
				"\" boundaryCondition=\"false\""~~ 
				If[IntNames[[i]]=="NoName", "", 
					" name=\""~~IntNames[[i]]~~"\""]~~
					" />", {i, 1, Length[BlockIntRows]}], 
			Table["<species id=\""~~ExtIDs[[i]]~~
				"\" compartment=\""~~Compartlist[[BlockExtRows[[i]]]]~~
				"\" boundaryCondition=\"true\""~~ 
				If[ExtNames[[i]]=="NoName", "", 
					" name=\""~~ExtNames[[i]]~~"\""]~~
					" />", {i, 1, Length[BlockExtRows]}]];,

		(*Else, the full SBML specification of each species as read from an SBML input file is
			transmitted, except that BoundaryCondition is adjusted according to 
			external/internal status used by Netsplitter*)
		BlockBC=Join[ConstantArray["false",Length[BlockIntRows]],ConstantArray["true",Length[BlockExtRows]]];
		BlockSpeclist=PutAttribute[SpeciesXML[[BlockRows]], "species", "boundaryCondition", BlockBC];
		SpeciesSBML=XMLtoText[BlockSpeclist];
		  ];
		  
	  If[ReactionXML == {}, 
	  	BlockReactNames=SBMLfixName[ReactNames[[BlockIntCols]]];
 		BlockReactlist = SBMLfixID[Reactlist[[BlockIntCols]]]; 
 		SBMLReverselist = SBMLfixID[Reverselist];
 		SBMLfluxes = FluxValues; SBMLfluxes[[All, 1]] = SBMLfixID[SBMLfluxes[[All, 1]]];
 		SBMLbounds = FluxBounds; SBMLbounds[[All, 1]] = SBMLfixID[SBMLbounds[[All, 1]]];
 		SBMLobjectives = ObjWeights; SBMLobjectives[[All, 1]] = SBMLfixID[SBMLobjectives[[All, 1]]];
 		
 		tofix = Intersection[BlockReactlist, Speclist];
		If[tofix != {}, tofixat = Flatten@Map[Position[BlockReactlist, #] &, tofix];
			BlockReactlist[[tofixat]] = Map[# <> "_RXN" &, BlockReactlist[[tofixat]]];
			tofixat = Flatten@Map[Position[SBMLReverselist, #] &, tofix];
			SBMLReverselist[[tofixat]] = Map[# <> "_RXN" &, SBMLReverselist[[tofixat]]];
			tofixat = Flatten@Map[Position[SBMLfluxes[[All,1]], #] &, tofix];
			SBMLfluxes[[tofixat,1]] = Map[# <> "_RXN" &, SBMLfluxes[[tofixat,1]]];
			tofixat = Flatten@Map[Position[SBMLbounds[[All,1]], #] &, tofix];
			SBMLbounds[[tofixat,1]] = Map[# <> "_RXN" &, SBMLbounds[[tofixat,1]]];
			tofixat = Flatten@Map[Position[SBMLobjectives[[All,1]], #] &, tofix];
			SBMLobjectives[[tofixat,1]] = Map[# <> "_RXN" &, SBMLobjectives[[tofixat,1]]];
		  ];		
 		

		BlockReactlist = Table[{BlockReactlist[[i]], 
			If[Count[SBMLReverselist, BlockReactlist[[i]]] > 0, True, False]},
			 {i, 1, Length[BlockReactlist]}];
 		Stoicol=(Transpose@Normal@BlockStoi);
		Reacts=Map[Flatten[Position[#,x_/;x<0]]&,Stoicol];
		Prods=Map[Flatten[Position[#,x_/;x>0]]&,Stoicol];
		Reactions=Table[{BlockReactlist[[i]],
			Transpose[{Speclist[[Reacts[[i]]]],-Stoicol[[i,Reacts[[i]]]]}],
			Transpose[{Speclist[[Prods[[i]]]],Stoicol[[i,Prods[[i]]]]}]},{i,1,Length[BlockReactlist]}];
		ReactionSBML=Flatten@Table[Join[{"<reaction id=\""~~ Reactions[[i,1,1]]~~"\"" ~~
			If[ Reactions[[i,1,2]],""," reversible=\"false\" "]~~
			If[BlockReactNames[[i]]=="NoName", ""," name=\"" ~~ BlockReactNames[[i]]~~"\""] ~~
			" >"},
			ReactantSBML[i],ProductSBML[i],
			KineticParamSBML[BlockReactlist[[i,1]],SBMLfluxes,SBMLbounds,SBMLobjectives],{"</reaction>"}],
				{i,1,Length[Reactions]}];,
		
	(*Else, the full SBML specification of each reaction as read from an SBML input file is
	transmitted, except that reversibility is adjusted according to that used by Netsplitter*)
		BlockRev=Table[If[Count[Reversibles,BlockIntCols[[i]]]>0,"true","false"],{i,1,Length[BlockIntCols]}];
		BlockReactlist=PutAttribute[ReactionXML[[BlockIntCols]], "reaction", "reversible", BlockRev];
		ReactionSBML=XMLtoText[BlockReactlist];
			
	   ];
	Export[Outfile,Join[Header,FuncDefs,UnitDefs,CompartTypes,SpecieTypes,
		{"<listOfCompartments>"},CompartmentSBL,{"</listOfCompartments>",
		"<listOfSpecies>"},SpeciesSBML,{"</listOfSpecies>"},
		Parameters,InitAss,Rules,Constraints,
		{"<listOfReactions>"},ReactionSBML,{"</listOfReactions>"},Events,{"</model>","</sbml>"}],"List"];

]]

End[] (* End Private Context *)

EndPackage[]
