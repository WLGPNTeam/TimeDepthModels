(* ::Package:: *)

(* ::Chapter:: *)
(*Visualisation*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[PlotDepthSection, PlotVelocity, PlotTimeSection, PlotDepthSectionWithWells]


PlotDepthSection::usage = 
"PlotDepthSection[horNHsorted]"


PlotVelocity::usage = 
"PlotVelocity[velModel]"


PlotTimeSection::usage = 
"PlotTimeSection[timeNH]"


PlotDepthSectionWithWells::usage = 
"PlotDepthSectionWithWells[horNHsorted, wells]"


(* ::Section:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


PlotDepthSection[horNHsorted_] := 
ListLinePlot[horNHsorted,
	GridLinesStyle -> Directive[Thick, Blue],
	FrameStyle -> Directive[Black, 18], 
	Filling -> Bottom, Frame -> True, ImageSize -> 800,
	PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horNHsorted]] - 1)],
	PlotLabel -> "Depth Section", 
	LabelStyle -> Directive[18, Bold, Gray]
]


PlotVelocity[velModel_] :=
ListContourPlot[Flatten[velModel, 2], 
	PlotTheme -> "Detailed"
]


PlotTimeSection[timeNH_] := 
ListLinePlot[timeNH,
	GridLinesStyle -> Directive[Thick, Blue],
	FrameStyle -> Directive[Black, 18], 
	Filling -> Bottom, Frame -> True, ImageSize -> 800,
	PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[timeNH]] - 1)],
	PlotLabel -> "Time Section", 
	LabelStyle -> Directive[18, Bold, Gray]
]


PlotDepthSectionWithWells[horNHsorted_, wells_] := 
ListLinePlot[horNHsorted,
	GridLines->{wells[[All, 2]], None},
	GridLinesStyle -> Directive[Thick, Blue],
	FrameStyle -> Directive[Black, 18], 
	Filling -> Bottom, Frame -> True, ImageSize -> 800,
	PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horNHsorted]] - 1)],
	PlotLabel -> "Depth Section with wells", 
	LabelStyle -> Directive[18, Bold, Gray]
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
