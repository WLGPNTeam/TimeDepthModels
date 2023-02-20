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
"PlotDepthSection[horizons, opts:OptionsPattern[ListLinePlot]]"


PlotVelocity::usage = 
"PlotVelocity[model, horizons]"


PlotTimeSection::usage = 
"PlotTimeSection[time, opts:OptionsPattern[ListLinePlot]]"


PlotHT::usage = 
"PlotHT[wellValues, lmSet, t]"


PlotVT::usage = 
"PlotVT[wellValues, lmSet, t]"


PlotTH::usage = 
"PlotTH[wellValues, lmSet, h]"


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


PlotHT[wellValues_, lmSet_, t_]:= 
Module[{
				plots,
				tmin,
				tmax,
				i
},
				tmin = Table[Min[wellValues[[i]][[All, 1]]], {i, Length[wellValues]}];
				tmax = Table[Max[wellValues[[i]][[All, 1]]], {i, Length[wellValues]}];
				plots = Table[Show[ListPlot[wellValues[[i]], ImageSize -> 500,
												PlotLabel -> StringJoin["Horizon ",ToString[i],". h = f(t)"],
												LabelStyle -> Directive[14, Gray],
												GridLines -> {wellValues[[i]][[All, 1]], wellValues[[i]][[All, 2]]}
									],
									Plot[lmSet[[i]][t], {t, tmin[[i]], tmax[[i]]}], 
									
									Frame -> True,
									FrameStyle -> Directive[14, Black],
									FrameLabel -> {"t, s", "h, m"}], {i, Length[wellValues]}
							];

				Return[plots]
	
]


PlotVT[wellValues_, lmSet_, t_]:= 
Module[{
				plots,
				tmin,
				tmax,
				i
},
				tmin = Table[Min[wellValues[[i]][[All,2]]], {i, Length[wellValues]}];
				tmax = Table[Max[wellValues[[i]][[All,2]]], {i, Length[wellValues]}];
				plots = Table[Show[ListPlot[wellValues[[i]][[All, 2;;3]], ImageSize -> 500,
												PlotLabel -> StringJoin["Horizon ",ToString[i],". v = f(t)"],
												LabelStyle -> Directive[14, Gray],
												GridLines -> {wellValues[[i]][[All,2]], wellValues[[i]][[All, 3]]}
									],
									Plot[lmSet[[i]][t], {t, tmin[[i]], tmax[[i]]}], 
									
									Frame -> True,
									FrameStyle -> Directive[14, Black],
									FrameLabel -> {"t, s", "v, m"}], {i, Length[wellValues]}
							];

				Return[plots]
]


PlotTH[wellValues_, lmSet_, h_]:= 
Module[{
				plots,
				hmin,
				hmax,
				i
},
				hmin = Table[Min[wellValues[[i]][[All, 1]]], {i, Length[wellValues]}];
				hmax = Table[Max[wellValues[[i]][[All, 1]]], {i, Length[wellValues]}];
				plots = Table[Show[ListPlot[wellValues[[i]], ImageSize -> 500,
												PlotLabel -> StringJoin["Horizon ", ToString[i],". t = f(h)"],
												LabelStyle -> Directive[14, Gray],
												GridLines -> {wellValues[[i]][[All, 1]], wellValues[[i]][[All, 2]]}
									],
									Plot[lmSet[[i]][h], {h, hmin[[i]], hmax[[i]]}], 
									
									Frame -> True,
									FrameStyle -> Directive[14, Black],
									FrameLabel -> {"h, m", "t, s"}], {i, Length[wellValues]}
							];

				Return[plots]
	
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
