(* ::Package:: *)

(* ::Chapter:: *)
(*Visualisation*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[PlotDepthSection, PlotVelocity, PlotTimeSection, PlotHT, PlotVT, PlotTH]


PlotDepthSection::usage = 
"PlotDepthSection[horNH, opts:OptionsPattern[ListLinePlot]]"


PlotVelocity::usage = 
"PlotVelocity[velModel, horNH, opts:OptionsPattern[{ListLinePlot, ListContourPlot}]]"


PlotTimeSection::usage = 
"PlotTimeSection[timeNH, opts:OptionsPattern[ListLinePlot]]"


PlotHT::usage = 
"PlotHT[dataWellsHT, lmSet, t]"


PlotVT::usage = 
"PlotVT[setVT, lmSet, t]"


PlotTH::usage = 
"PlotTH[dataWellsTH, lmSet, h]"


(* ::Section:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


PlotDepthSection[horNH_, opts:OptionsPattern[ListLinePlot]] := 
ListLinePlot[horNH, opts]


PlotDepthSection[horNH_, datasetWells_, opts:OptionsPattern[ListLinePlot]] := 
ListLinePlot[horNH, opts]


PlotVelocity[velModel_, horNH_, opts:OptionsPattern[{ListLinePlot, ListContourPlot}]] :=
Show[ListContourPlot[Flatten[velModel, 2][[All, 2 ;; 4]], 
						FilterRules[opts, Options[ListContourPlot]],
						ColorFunction -> ColorData[{"RedBlueTones", "Reverse"}],
						PlotLegends -> BarLegend[Automatic, LegendLabel -> "Vel, m/s"],
						PlotLabel -> "Velocity Distribution"
									],
			ListLinePlot[horNH, 
							FilterRules[opts,Options[ListLinePlot]],
							PlotStyle -> {Directive[Thickness[0.005], Black]},
							PlotLabel->"Velocity Distribution",
							PlotLabels -> Map["Hor " <> ToString[#] &, (Range[Length[horNH]] - 1)],
							LabelStyle -> Directive[14, Gray],  
							ImageSize -> 500
									], 
			Frame -> True,
			FrameStyle -> Directive[14, Black],
			PlotRangePadding -> {{Scaled[0.05], Scaled[0.2]}, {Scaled[0.05], Scaled[0.05]}}
] 


PlotTimeSection[timeNH_, opts:OptionsPattern[ListLinePlot]] :=
Module[{
				timeNHreverseT, (*dont know how to reverse positive axis t on the plot, so there is such a solation*)
				i,
				j
},
				timeNHreverseT = Table[{timeNH[[i, j, 1]], -timeNH[[i, j, 2]]}, {i, Length[timeNH]}, {j, Length[timeNH[[i]]]}];
				ListLinePlot[timeNHreverseT, opts]
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
												PlotLabel -> StringJoin["Horizon ",ToString[i],". h = f(t)"],
												LabelStyle -> Directive[14, Gray],
												GridLines -> {dataWellsHT[[i]][[All, 1]], dataWellsHT[[i]][[All, 2]]}
									],
									Plot[lmSet[[i]][t], {t, tmin[[i]], tmax[[i]]}], 
									
									Frame -> True,
									FrameStyle -> Directive[14, Black],
									FrameLabel -> {"t, s", "h, m"}], {i, Length[dataWellsHT]}
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
									FrameLabel -> {"t, s", "v, m"}], {i, Length[setVT]}
							];

				Return[plots]
]


PlotTH[dataWellsTH_, lmSet_, h_]:= 
Module[{
				plots,
				hmin,
				hmax,
				i
},
				hmin = Table[Min[dataWellsTH[[i]][[All, 1]]], {i, Length[dataWellsTH]}];
				hmax = Table[Max[dataWellsTH[[i]][[All, 1]]], {i, Length[dataWellsTH]}];
				plots = Table[Show[ListPlot[dataWellsTH[[i]], ImageSize -> 500,
												PlotLabel -> StringJoin["Horizon ",ToString[i],". t = f(h)"],
												LabelStyle -> Directive[14, Gray],
												GridLines -> {dataWellsTH[[i]][[All, 2]], dataWellsTH[[i]][[All, 1]]}
									],
									Plot[lmSet[[i]][h], {h, hmin[[i]], hmax[[i]]}], 
									
									Frame -> True,
									FrameStyle -> Directive[14, Black],
									FrameLabel -> {"h, m", "t, s"}], {i, Length[dataWellsTH]}
							];

				Return[plots]
	
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
