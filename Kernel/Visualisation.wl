(* ::Package:: *)

(* ::Chapter:: *)
(*Visualisation*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[PlotDepthSection, PlotVelocity, PlotTimeSection, LinearRegressionPlots, PlotPredictionModel2D]


PlotDepthSection::usage = 
"PlotDepthSection[horizons, opts:OptionsPattern[ListLinePlot]]"


PlotVelocity::usage = 
"PlotVelocity[model, horizons]"


PlotTimeSection::usage = 
"PlotTimeSection[time, opts:OptionsPattern[ListLinePlot]]"


LinearRegressionPlots::usage =
"LinearRegressionPlots[points, lmSet, variable, plotName, axesName]"


PlotPredictionModel2D::usage = 
"PlotPredictionModel2D[ds, allHorizons, opts:OptionsPattern[ListLinePlot]]"


PlotPredictionModel2D::usage = 
"PlotPredictionModel2D[allHorizons, opts:OptionsPattern[ListLinePlot]]"


(* ::Section:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


PlotDepthSection[horizons_, opts:OptionsPattern[ListLinePlot]] := 
ListLinePlot[horizons, opts]


PlotDepthSection[horizons_, ds_, opts:OptionsPattern[ListLinePlot]] := 
ListLinePlot[horizons, opts]


PlotVelocity[model_, horizons_ (*there are must be 2 OptionsPatterns*)] :=
Show[ListContourPlot[Flatten[model, 2][[All, 2 ;; 4]], 
						{ColorFunction -> ColorData[{"RedBlueTones", "Reverse"}],
						PlotLegends -> BarLegend[Automatic, LegendLabel -> "v, \:043c/\:0441"],
						PlotLabel -> "\:0421\:043a\:043e\:0440\:043e\:0441\:0442\:043d\:043e\:0439 \:0440\:0430\:0437\:0440\:0435\:0437",
						Contours -> {Automatic, 20},
						ContourStyle -> None}],
			ListLinePlot[horizons, {PlotStyle -> {Directive[Thickness[0.001], Black]}, PlotLabels -> Placed[Map["Hor " <> ToString[#] &, (Range[Length[horizons]] - 1)], {Above}],ImageSize -> 500}], {Frame -> True,(*
			FrameStyle -> Directive[14, Black],*) FrameLabel->{"x, \:043c","z, \:043c"},
			PlotRangePadding -> {{Scaled[0.05], Scaled[0.2]}, {Scaled[0.05], Scaled[0.05]}}}
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


LinearRegressionPlots[points_, lmSet_, variable_, plotName_, axesName_ (*need 3 OptionsPatterns*)]:= 
Module[{
				plots,
				min,
				max,
				i
},
				
				min = Table[Min[points[[i]][[All, 1]]], {i, Length[points]}];
				max = Table[Max[points[[i]][[All, 1]]], {i, Length[points]}];
				plots = Table[Show[ListPlot[points[[i]], PlotMarkers->{Automatic, Scaled[0.01]}, ImageSize -> 500,
												PlotLabel -> StringJoin[ "\:041c\:0435\:0442\:043e\:0434:", plotName],
												LabelStyle -> Directive[14, Gray](*,
												GridLines -> {points[[i]][[All, 1]], points[[i]][[All, 2]]}*)
									],
									Plot[lmSet[[i]][variable], {variable, min[[i]], max[[i]]}, PlotLegends -> Normal[lmSet[[i]]]], 
									
									Frame -> True,
									FrameStyle -> Directive[14, Black],
									FrameLabel -> axesName], {i, Length[points]}									
							];
							
							Return[plots]
]


LinearRegressionPlots[points_, lmSet_, variable_, plotName_, axesName_, hnames_(*need 3 OptionsPatterns*)]:= 
Module[{
				plots,
				min,
				max,
				i
},
				
				min = Table[Min[points[[i]][[All, 1]]], {i, Length[points]}];
				max = Table[Max[points[[i]][[All, 1]]], {i, Length[points]}];
				plots = Table[Show[ListPlot[points[[i]], PlotMarkers->{Automatic, Scaled[0.015]}, ImageSize -> 500,
												PlotLabel -> StringJoin[ "\:0413\:043e\:0440\:0438\:0437\:043e\:043d\:0442: ", ToString[hnames[[i]]], ". \:041c\:0435\:0442\:043e\:0434: ", plotName],
												LabelStyle -> Directive[14, Gray](*,
												GridLines -> {points[[i]][[All, 1]], points[[i]][[All, 2]]}*)
									],
									Plot[lmSet[[i]][variable], {variable, min[[i]], max[[i]]}, PlotLegends -> Normal[lmSet[[i]]]], 
									
									Frame -> True,
									FrameStyle -> Directive[14, Black],
									FrameLabel -> axesName], {i, Length[points]}									
							];
							
							Return[plots]
]


PlotPredictionModel2D[allHorizons_, opts:OptionsPattern[ListLinePlot]]:=
Show[Table[ListLinePlot[allHorizons[[i]], opts[[i]]], {i, Length[allHorizons]}]]


PlotPredictionModel2D[ds_, allHorizons_, opts:OptionsPattern[ListLinePlot]]:=
Show[Table[ListLinePlot[allHorizons[[i]], opts[[i]]], {i, Length[allHorizons]}]]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
