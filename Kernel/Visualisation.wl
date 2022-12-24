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
"PlotVelocity[velModel, horNHsorted]"


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
					FrameStyle -> Directive[Black, 18], 
					Filling -> Bottom, Frame -> True, ImageSize -> 800,
					PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horNHsorted]] - 1)],
					PlotLabel -> "Depth Section", 
					LabelStyle -> Directive[18, Bold, Gray]
]


PlotVelocity[velModel_, horNHsorted_] :=
Show[ListContourPlot[Flatten[velModel, 2], 
										ColorFunction -> ColorData[{"RedBlueTones", "Reverse"}],
										PlotLegends -> BarLegend[Automatic, LegendLabel -> "Vel, m/s"],
										PlotLabel -> "Velocity Distribution"
										],
			ListLinePlot[horNHsorted, 
									PlotStyle -> {Directive[Thickness[0.005], Black]},
									PlotLabel->"Velocity Distribution",
									PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horNHsorted]] - 1)],
									LabelStyle -> Directive[14, Bold, Gray],  
									ImageSize -> 800
									], 
			Frame -> True,
			FrameStyle -> Directive[Black, 12],
			PlotRangePadding -> {{Scaled[0.05], Scaled[0.2]}, {Scaled[0.05], Scaled[0.05]}}
] 


PlotTimeSection[timeNH_] := 
ListLinePlot[timeNH,
					FrameStyle -> Directive[Black, 18], 
					Filling -> Bottom, Frame -> True, ImageSize -> 800,
					PlotLabels -> Map["t " <> ToString[#] &, (Range[Length[timeNH]] - 1)],
					PlotLabel -> "Time Section", 
					LabelStyle -> Directive[18, Bold, Gray]
]


PlotDepthSectionWithWells[horNHsorted_, wells_] := 
ListLinePlot[horNHsorted,
					GridLines->{wells[[All, 2]], None},
					GridLinesStyle -> Directive[Thick, Grey],
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
