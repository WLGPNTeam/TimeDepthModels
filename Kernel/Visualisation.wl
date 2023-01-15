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


PlotHT::usage = 
"PlotHT[dataWellsHT, lmSet, t]"


PlotVT::usage = 
"PlotVT[setVT, lmSet, t]"


(* ::Section:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


PlotDepthSection[horNHsorted_,opts:OptionsPattern[ListLinePlot]] := 
ListLinePlot[horNHsorted,
				opts,
				FrameStyle -> Directive[14, Black], 
				Filling -> Bottom, Frame -> True, ImageSize -> 500,
				PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horNHsorted]] - 1)],
				PlotLabel -> "Depth Section", 
				LabelStyle -> Directive[14, Gray]
]


PlotDepthSection[horNHsorted_,datasetWells_, opts:OptionsPattern[ListLinePlot]] := 
ListLinePlot[horNHsorted,
				opts,
				GridLines->{DeleteDuplicates[Normal[datasetWells[All, "x"]]], None},
				GridLinesStyle -> Directive[Thick, Gray],
				FrameStyle -> Directive[14, Black], 
				Filling -> Bottom, Frame -> True, ImageSize -> 500,
				PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horNHsorted]] - 1)],
				PlotLabel -> "Depth Section with wells", 
				LabelStyle -> Directive[14, Gray]
]


PlotVelocity[velModel_, horNHsorted_,opts:OptionsPattern[{ListLinePlot,ListContourPlot}]] :=
Show[ListContourPlot[Flatten[velModel, 2][[All, 2 ;; 4]], 
						FilterRules[opts,Options[ListContourPlot]],
						ColorFunction -> ColorData[{"RedBlueTones", "Reverse"}],
						PlotLegends -> BarLegend[Automatic, LegendLabel -> "Vel, m/s"],
						PlotLabel -> "Velocity Distribution"
									],
			ListLinePlot[horNHsorted, 
							FilterRules[opts,Options[ListLinePlot]],
							PlotStyle -> {Directive[Thickness[0.005], Black]},
							PlotLabel->"Velocity Distribution",
							PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horNHsorted]] - 1)],
							LabelStyle -> Directive[14, Gray],  
							ImageSize -> 500
									], 
			Frame -> True,
			FrameStyle -> Directive[14, Black],
			PlotRangePadding -> {{Scaled[0.05], Scaled[0.2]}, {Scaled[0.05], Scaled[0.05]}}
] 


PlotTimeSection[timeNH_] :=
Module[{
				timeNHreverseT, (*dont know how to reverse positive axis t on the plot, so there is such a solation*)
				i,
				j
},
				timeNHreverseT = Table[{timeNH[[i, j, 1]], -timeNH[[i, j, 2]]}, {i, Length[timeNH]}, {j, Length[timeNH[[i]]]}];
				ListLinePlot[timeNHreverseT,
								FrameStyle -> Directive[14, Black], 
								Filling -> Bottom, Frame -> True, ImageSize -> 500,
								PlotLabels -> Map["t " <> ToString[#] &, (Range[Length[timeNH]] - 1)],
								PlotLabel -> "Time Section", 
								LabelStyle -> Directive[14, Gray]
							]
]


PlotDepthSectionWithWells[horNHsorted_, datasetWells_] := 
ListLinePlot[horNHsorted,
				GridLines->{DeleteDuplicates[Normal[datasetWells[All, "x"]]], None},
				GridLinesStyle -> Directive[Thick, Gray],
				FrameStyle -> Directive[14, Black], 
				Filling -> Bottom, Frame -> True, ImageSize -> 500,
				PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horNHsorted]] - 1)],
				PlotLabel -> "Depth Section with wells", 
				LabelStyle -> Directive[14, Gray]
]


PlotHT[dataWellsHT_, lmSet_, t_]:= 
Module[{
				plots,
				tmin,
				tmax,
				i
},
				tmin = Table[Min[dataWellsHT[[i]][[All, 1]]], {i, Length[dataWellsHT]}];
				tmax = Table[Max[dataWellsHT[[i]][[All, 1]]], {i, Length[dataWellsHT]}];
				plots = Table[Show[ListPlot[dataWellsHT[[i]], ImageSize -> 500,
												PlotLabel -> StringJoin["Horizon ",ToString[i - 1],". h = f(t)"],
												LabelStyle -> Directive[14, Gray],
												GridLines -> {dataWellsHT[[i]][[All, 1]],dataWellsHT[[i]][[All, 2]]}
									],
									Plot[lmSet[[i]][t], {t, tmin[[i]], tmax[[i]]}], 
									
									Frame -> True,
									FrameStyle -> Directive[14, Black],
									FrameLabel -> {"t, s", "h, m"}], {i, 2, Length[dataWellsHT]}
							];

				Return[plots]
	
]


PlotVT[setVT_, lmSet_, t_]:= 
Module[{
				plots,
				tmin,
				tmax,
				i
},
				tmin = Table[Min[setVT[[i]][[All,2]]], {i, Length[setVT]}];
				tmax = Table[Max[setVT[[i]][[All,2]]], {i, Length[setVT]}];
				plots = Table[Show[ListPlot[setVT[[i]][[All, 2;;3]], ImageSize -> 500,
												PlotLabel -> StringJoin["Horizon ",ToString[i - 1],". v = f(t)"],
												LabelStyle -> Directive[14, Gray],
												GridLines -> {setVT[[i]][[All,2]], setVT[[i]][[All, 3]]}
									],
									Plot[lmSet[[i]][t], {t, tmin[[i]], tmax[[i]]}], 
									
									Frame -> True,
									FrameStyle -> Directive[14, Black],
									FrameLabel -> {"t, s", "v, m"}], {i, 2, Length[setVT]}
							];

				Return[plots]
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
