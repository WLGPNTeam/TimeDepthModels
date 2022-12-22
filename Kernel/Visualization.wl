(* ::Package:: *)

(* ::Chapter:: *)
(*BuildSection*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[PlotSection]


PlotSection::usage = 
"PlotSection[horNHsorted]"


(* ::Section:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


PlotSection[horNHsorted_, numOfLayers_Integer] := 
ListLinePlot[horNHsorted, 
	GridLinesStyle -> Directive[Thick, Blue], 
	FrameStyle -> Directive[Black, 18], 
	Filling -> Bottom, Frame -> True, ImageSize -> 800, 
	PlotLabels -> Map["Hor " <> ToString[#] &, (Range[numOfLayers + 1] - 1)], 
	PlotLabel -> "Depth Model", 
	LabelStyle -> Directive[18, Bold, Gray]
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
