(* ::Package:: *)

(* ::Chapter:: *)
(*SupportFunctions*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[subsets, CheckVelocity, ImportData, ExportPredictedDepths, ExportPredictionModelsParameters, SelectHorizons, ExportDepthSection]


subsets::usage = 
"subsets[set]"


CheckVelocity::usage = 
"CheckVelocity[ds]"


ImportData::usage = 
"ImportData[file]"


ExportPredictedDepths::usage =
"ExportPredictedDepths[f, horizonsNames, wellNames, methodsNames, depthObjective, depthPredicted]"


ExportPredictionModelsParameters::usage = 
"ExportPredictionModelsParameters[f, res, RSquared, RMS, methodsNames]"


SelectHorizons::usage = 
"SelectHorizons[ds, set]"


ExportDepthSection::usage = 
"ExportDepthSection[f, horizons]"


(* ::Section::Closed:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


subsets[set_]:= Table[Join[{set[[1]]}, Cases[set[[2 ;; -1]], {__, i}][[All, 1]]], {i, 4}];


CheckVelocity[ds_]:= 
Module[{
                i,
                j,
                dt,
                dh,
                table,
                wellValuesSelected,
                wellCount,
                horizonsCount,
                result                
},
                wellCount = Normal[Values[ds]][[-1, 1]]; (*count wells*)
                horizonsCount = Max[Normal[Values[ds]][[All, 3]]]; (*count horizons*)
                table = Table[{0, 0, 0, 0, 0}, {i, 1, wellCount}, {j, 1, horizonsCount - 1}]; (*there will be results*)

                For[i = 1, i <= wellCount, i++, (*take i well*)
                    wellValuesSelected = Select[Normal[Values[ds]], #[[1]] == i& ]; (*take all data from this well (x, horizons, depths, times)*)
                    For[j = 1, j <= horizonsCount - 1, j++, (*go through layers*)
                        dh = -(wellValuesSelected [[j + 1, 4]] - wellValuesSelected [[j, 4]]);  (*in turn, calculate the layer thickness,*)  
                        dt = wellValuesSelected [[j + 1, 5]] - wellValuesSelected [[j, 5]];                    (*time*)
                        table[[i]][[j]] = {i, j, dh, dt, Abs[2 dh/dt]}                                    (*and interval speed*)
                        ] (*go to next layer*)
                    ];
                
                result = Dataset[Map[<|"well" -> #[[1]], "layer" -> #[[2]], "dh" -> #[[3]], "dt" -> #[[4]], "Vint" -> #[[5]]|>&, Flatten[table, 1]]]; (*make dataset*)
				
                Return[<|"ds" -> result|>]
]


ImportData[f_] := 
Module[{
                dataset,
                data
},
                data = Import[f, "Table"];
                dataset = Dataset[Map[<|Table[data[[1, i]] -> #[[i]],  {i, Length[data[[1]]]}]|>&, data[[2 ;; -1]]]];
 
                Return[dataset]
]


ExportPredictedDepths[f_, horizonsNames_, wellNames_, methodsNames_, depthObjective_, depthPredicted_] :=
Module[{
                data,
				ds,
				cols
},                         
				data = Flatten[Table[Flatten[{wellNames[[i, j]], horizonsNames[[i]], depthObjective[[i, j]], depthPredicted[[All, i,  j, 2]]}, 1], {k, 1}, {i, Length[wellNames]}, {j, Length[wellNames[[i]]]}], 2];
				cols = Flatten[{"well", "horizon", "z, \:043c", methodsNames}, 1];
				ds = Dataset[Map[<|Table[cols[[i]] -> #[[i]],  {i, Length[data[[1]]]}]|>&, data]];
				Export[f, ds];
 
Return[ds]
]


ExportPredictionModelsParameters[f_, res_, RSquared_, RMS_, methodsNames_] :=
Module[{
                 ds,
                 k
},                         
				(*ds = Dataset[Table[<|methodsNames[[i]]->Table[<|"fit"->Normal[lmSet[[i,j]]],"RSquared"->RSquared[[i,j]]|>, {j, Length[lmSet[[i]]]}]|>, {i, Length[methodsNames]}]];*)
				k = {{"P1art", "B"},{"P1art", "B"},{"P1art", "B"}}; (*!!!*)
				ds = Dataset[Table[<|"method" -> methodsNames[[i]], "horizon" -> (*k[[i,j]]*) "B", "fit / average velocity, \:043c/\:0441    " -> Normal[res[[i, j]]], "RSquared" -> RSquared[[i, j]], "RMS, \:043c" -> RMS[[i, j]]|>, {i, Length[methodsNames]}, {j, Length[res[[i]]]}]];
				Export[f, ds, "Table"];
 
				Return[ds]
]


SelectHorizons[ds_, set_]:= Select[ds, MemberQ[set, #horizon]&] (*S s !!!*)


ExportDepthSection[f_, horizons_]:=
Module[{
			ds,
			i,
			j
},
			ds = Dataset[Map[<|"horizon" -> #[[1]], "x,m" -> #[[2]], "h,m" -> #[[3]]|>&, Flatten[Table[{i, horizons[[i, j, 1]], horizons[[i, j, 2]]}, {i, Length[horizons]}, {j, Length[horizons[[i]]]}], 1]]];
			
			Export[f, ds]
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
