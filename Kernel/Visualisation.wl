(* ::Package:: *)

(* ::Chapter:: *)
(*Visualisation*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[PlotSection, PlotVelocity]


PlotSection::usage = 
"PlotSection[horNHsorted]"


PlotVelocity::usage = 
"PlotVelocity[velModel]"


(* ::Section:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


PlotSection[horNHsorted_] := 
ListLinePlot[horNHsorted,
	GridLinesStyle -> Directive[Thick, Blue],
	FrameStyle -> Directive[Black, 18], 
	Filling -> Bottom, Frame -> True, ImageSize -> 800,
	PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horNHsorted]] - 1)],
	PlotLabel -> "Depth Model", 
	LabelStyle -> Directive[18, Bold, Gray]
]


PlotVelocity[velModel_] :=
ListContourPlot[Flatten[velModel, 2], 
	PlotTheme -> "Detailed"
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]