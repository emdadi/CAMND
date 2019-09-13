(* Mathematica Package *)

(* COPYRIGHT (C) 2010 Wynand Verwoerd 
under the terms of the GNU General Public License V3
This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
*)

(* This package run by init.m declares and/or initialises variables semi-globally for 
   use in all packages used by NetSplitter.
   Remember to add the context {"NetSplitter`"} as a second argument of BeginPackage
   in all packages that use any of them. *)

BeginPackage["NetSplitter`"]
(* SymbolName::usage="text" adds a hover message *) 

Foundtinge::usage="Background colour of blocks that were found "
Chosentinge::usage="Background colour of selected blocks"
Targetinge::usage="Background colour of target metabolites"
PrintColumns::usage="No of columns in metabolite listings"
Effpower::usage="Power law exponent for calculating split efficacy"

(* Configurations *)
NSPLVersion=" Version 1.3.2  December 2011  ";
Targetinge = Hue[.925, 0.4, 1]; Foundtinge = 
 Hue[.5833, 0.5, 1]; Chosentinge = Hue[.4166, 0.5, 1];
PrintColumns=3; PrintMultHorPages=True; SafeDraw=False;
Attention = Sound[SoundNote[#, .4, "Whistle"] & /@ {7, 16, 12}];
Crisis = Sound[SoundNote[5, 1, "Bird"]] ;
Effpower := 0.25 (CNodes - Length[StructuralExternalrows])^0.5; (* OR 2<const<8 *)
DefaultCompartment="Cytosol";
SetAttributes[StringSplit,Listable];

(* Main dialog choices *)
(* Edit the 1st StartupDirectory line for a custom path, then swap with 2nd  *)
(* NOTE: Use double backslash \\ in explicit directory path for Windows *)
StartupDirectory="H:\\My Documents\\Research\\Projects\\MetabolicModels\\";
StartupDirectory=NotebookDirectory[];  (* DEFAULT *)
DistChoice=SokalSneathDissimilarity;LinkChoice="Single";ReverseThem=False;
Manygreys=1; GreyMax=0.5;GreyLimit;GroupLevel;GroupLevelRange={0.3,0.7};ConnexMax=8;

(* Startup state specification *)
StartTime=0;CalculationTime=0;Iterationcount;
SplitDone = False; Convert = False; Dropstring = "";
Refusals={};Stepcount=0;RepeatMerge=False;OldMerge={};NoMatch={};
Status = "Waiting for action";DecisionTrace;MergeTrace;

(* Files and data *)
Infile; Inpath;
Exfile = StartupDirectory; Targetfile = StartupDirectory; Fluxfile = StartupDirectory; 
Outfile="ExternalMetabolites.txt";
FluxUnits="mmol_per_gDW_per_min";
flotype="mole";floscale="-3";DWtype="gram";DWscale="0";timefactor="0.016667";
S;CR;RC;Complist;CompNames;CompIDs;
Reactlist;ReactNames;Reverselist;Reversibles;CNodes;RNodes;
Externals;Externalrows;StExrows;InputExternalrows={};OutputExternalrows={};
StructuralExternalrows;ConnectionExternalrows;
BasicExternals;StoichioExternals;
PreMergeExternals;ChosenExternals;
Internalrows;InternalComplist;InternalCNodes;Connex;
DropComps;OrphanRows;
TargetCompounds;TargetIDs;
Compartlist;CompartIDs;CompartOuters;CompartDims;CompartConstant;CompartSizes;CompartUnits;


(* Blocking variables *)
Nameseq;Sinks;Sources;CrossBlocks;IsolationBlock;
SelectedBlocks;
ShuffleDAG;ConsolBlocking;
FinalBlocks;FinalBlockCount;BlockCount;
IntrowsperBlock;ColsperBlock;DAG;
Sourceseq;Sinkseq;
NewExternals;
BlockSinks;BlockChoice;BlockCompart;BlockHue;
StartEff;PreMergeEff;PostMergeEff;

EndPackage[]

