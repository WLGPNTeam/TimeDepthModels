(* ::Package:: *)

(* ::Chapter:: *)
(*Visualisation*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[PlotDepthSection, PlotVelocity, PlotTimeSection, LinearRegressionPlots]


PlotDepthSection::usage = 
"PlotDepthSection[horizons, opts:OptionsPattern[ListLinePlot]]"


PlotVelocity::usage = 
"PlotVelocity[model, horizons]"


PlotTimeSection::usage = 
"PlotTimeSection[time, opts:OptionsPattern[ListLinePlot]]"


LinearRegressionPlots::usage =
"LinearRegressionPlots[points, lmSet, variable, plotName, axesName]"


(* ::Section:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


PlotDepthSection[horizons_, opts:OptionsPattern[ListLinePlot]] := 
ListLinePlot[horizons, opts]


PlotDepthSection[horizons_, datasetWells_, opts:OptionsPattern[ListLinePlot]] := 
ListLinePlot[horizons, opts]


PlotVelocity[model_, horizons_] :=
Show[ListContourPlot[Flatten[model, 2][[All, 2 ;; 4]], 
						ColorFunction -> ColorData[{"RedBlueTones", "Reverse"}],
						PlotLegends -> BarLegend[Automatic, LegendLabel -> "v, m/s"],
						PlotLabel -> "Velocity Distribution",
						Contours -> {Automatic, 20}],
			ListLinePlot[horizons, 
							PlotStyle -> {Directive[Thickness[0.005], Black]},
							PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horizons]] - 1)],
							ImageSize -> 500], 
			Frame -> True,
			FrameStyle -> Directive[14, Black],
			PlotRangePadding -> {{Scaled[0.05], Scaled[0.2]}, {Scaled[0.05], Scaled[0.05]}}
] 


PlotTimeSection[time_, opts:OptionsPattern[ListLinePlot]] :=
Module[{
				timeReverseT, (*dont know how to reverse positive axis t on the plot, so there is such a solation*)
				i,
				j
},
				timeReverseT = Table[{time[[i, j, 1]], -time[[i, j, 2]]}, {i, Length[time]}, {j, Length[time[[i]]]}];
				ListLinePlot[timeReverseT, opts]
]


LinearRegressionPlots[points_, lmSet_, variable_, plotName_, axesName_]:= 
Module[{
				plots,
				min,
				max,
				i
},
				min = Table[Min[points[[i]][[All, 1]]], {i, Length[points]}];
				max = Table[Max[points[[i]][[All, 1]]], {i, Length[points]}];
				plots = Table[Show[ListPlot[points[[i]], ImageSize -> 500,
												PlotLabel -> StringJoin["Horizon ", ToString[i],". ", plotName],
												LabelStyle -> Directive[14, Gray],
												GridLines -> {points[[i]][[All, 1]], points[[i]][[All, 2]]}
									],
									Plot[lmSet[[i]][variable], {variable, min[[i]], max[[i]]}, PlotLegends -> NumberForm[Normal[lmSet[[i]]], 3], 
									
									Frame -> True,
									FrameStyle -> Directive[14, Black],
									FrameLabel -> axesName]], {i, Length[points]}									
							];
							
							Return[plots]

]




(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
