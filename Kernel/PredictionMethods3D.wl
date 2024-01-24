(* ::Package:: *)

(* ::Chapter:: *)
(*PredictionMethods3D*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[HTMethod3D, VTMethod3D, THMethod3D, VaveMethod3D, dHdTMethod3D, dVdTMethod3D, dTdHMethod3D, dVaveMethod3D]


HTMethod3D::usage = 
"HTMethod3D[ds, timeMaps, var]"


VTMethod3D::usage = 
"VTMethod3D[ds, timeMaps, var]"


THMethod3D::usage = 
"THMethod3D[ds, timeMaps, var]"


VaveMethod3D::usage = 
"VaveMethod3D[ds, timeMaps]"


dHdTMethod3D::usage = 
"dHdTMethod3D[ref, ds, timeMaps, var]"


dVdTMethod3D::usage = 
"dVdTMethod3D[ref, ds, timeMaps, var]"


dTdHMethod3D::usage = 
"dTdHMethod3D[ref, ds, timeMaps, var]"


dVaveMethod3D::usage = 
"dVaveMethod3D[ref, ds, timeMaps, var]"


(* ::Section::Closed:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


HTMethod3D[ds_, timeMaps_, var_]:= 
Module[{
                i,
                j,
                xy,
                values,
                lmSet,
                fits,
                lmParametres,
                depthObjective,
                depthPredicted,
                errors,
                iErrors,
                interpolationErrors,
                result
},               
               values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "time", "depth"}]]]], {i, Normal[DeleteDuplicates[ds[[All, "horizon"]]]]}];  (*time, depth*)
				values[[All, All, 2]] /= 2000;   (*divide time*)     
                If[values[[1, 1, 3]] < 0, values[[All, All, 3]] *= -1];   (*transform depths to positive values for easy calculations and presentable plots*) 
                
                lmSet = Table[LinearModelFit[values[[i]][[All, 2 ;; 3]], var, var], {i,  Length[values]}];	(*evaluate linear models h(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
	  
	            fits  = Table[Table[{timeMaps[[i, j, 1]], timeMaps[[i, j, 2]], -(lmParametres[[i, 2]] * timeMaps[[i, j, 3]]/2000 + lmParametres[[i, 1]])}, {j, Length[timeMaps[[i]]]}], {i, Length[timeMaps]}];	           

				xy = Table[Values[Normal[ds][[All, {"x", "y"}]]], {i, Length[DeleteDuplicates[ds[[All, "horizon"]]]]}]; (*x y of wells*)

                (*find errors between wells values and predicted by model*)
				
				If[values[[1, 1, 3]] < 0, depthObjective = values[[All, All, {1, 3}]], depthObjective = Table[Map[{#[[1]], -#[[3]]}&, values[[i]]], {i, Length[values]}]];
                depthPredicted = Table[Table[{values[[i, j, 1]], -(lmParametres[[i, 2]] * values[[i, j, 2]]+ lmParametres[[i, 1]])}, {j, Length[values[[i]]]}], {i, Length[values]}];
				errors = Table[Table[{values[[i, j, 1]], xy[[i, j, 1]], xy[[i, j, 2]], (depthObjective[[i, j, 2]] - depthPredicted[[i, j, 2]])}, {j,  Length[values[[i]]]}], {i, Length[values]}];	
                
                  (*interpolatioan and adding errors*)
		(*\:0434\:043e\:0431\:0430\:0432\:0438\:0442\:044c \:0434\:0432\:0435 \:0441\:0442\:0440\:043e\:043a\:0438 \:043d\:0438\:0436\:043d\:0438\:0435 \:0432 \:043e\:0441\:0442\:0430\:043b\:044c\:043d\:044b\:0435 \:0431\:043b\:043e\:043a\:0438*)
                iErrors = Table[Interpolation[DeleteDuplicates[errors[[i, All, 2 ;; -1]]], InterpolationOrder -> 1], {i, Length[errors]}];
                result = Table[{#[[1]]*1., #[[2]]*1., #[[3]] + iErrors[[i]][#[[1]], #[[2]]]}&/@fits[[i]], {i, Length[fits]}];
                
                Return[<|"depthObjective" -> depthObjective,
                            "depthPredicted" -> depthPredicted,
                            "wellValues" -> values,
                            "result" -> result,
                            "lmSet" -> lmSet,
                            "lmParametres" -> lmParametres,
                            "fits" -> fits,
                            "errors" -> errors,
                            "RMSError" -> Map[StandardDeviation[#[[All, 4]]]&, errors],
                            "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


VTMethod3D[ds_, timeMaps_, var_]:= 
Module[{
                i,
                j,
                xy,
                table,                
                lmParametres,
                values,
                lmSet,
                fits,
                depthObjective,                
                errors,
                interpolationErrors,
                depthPredicted,
                result                
},               	          
	            values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "time", "depth"}]]]], {i, Normal[DeleteDuplicates[ds[[All, "horizon"]]]]}];  (*time, depth*)
				values[[All, All, 2]] /= 2000; (*divide time*)
	            If[values[[1, 1, 3]] < 0, values[[All, All, 3]] *= -1];   (*transform depths to positive values for easy calculations and presentable plots*) 
                
                table = Table[Table[{values[[i, j, 2]], values[[i, j, 3]]/values[[i, j, 2]]}, {j, Length[values[[i]]]}], {i, Length[values]}]; (*{time, velocity}*)
           
		        lmSet = Table[LinearModelFit[table[[i]], var, var], {i, Length[table]}]; (*evaluate linear models v(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
                
                fits = Table[Table[{timeMaps[[i, j, 1]], timeMaps[[i, j, 2]], -(lmParametres[[i, 2]] * (timeMaps[[i, j, 3]]/2000)^2  + lmParametres[[i, 1]] * timeMaps[[i, j, 3]]/2000)}, {j, Length[timeMaps[[i]]]}], {i, Length[timeMaps]}];(*evaluate fits (x, h = at^2+bt). make data suitable for ListLinePlot *)
                        
				xy = Table[Values[Normal[ds][[All, {"x", "y"}]]], {i, Length[DeleteDuplicates[ds[[All, "horizon"]]]]}]; (*x y of wells*)

                (*find errors between wells values and predicted by model*)
				
				If[values[[1, 1, 3]] < 0, depthObjective = values[[All, All, {1, 3}]], depthObjective = Table[Map[{#[[1]], -#[[3]]}&, values[[i]]], {i, Length[values]}]];
                depthPredicted = Table[Table[{values[[i, j, 1]],-(lmParametres[[i,2]] * values[[i, j, 2]]^2 + lmParametres[[i, 1]] * values[[i, j, 2]])}, {j, Length[values[[i]]]}], {i, Length[values]}];
				errors = Table[Table[{values[[i, j, 1]], xy[[i, j, 1 ;; 2]], (depthObjective - depthPredicted)[[i, j, 2]]}, {j,  Length[values[[i]]]}], {i, Length[values]}];	
                
                (*interpolatioan and adding error*)


                Return[<|"depthObjective" -> depthObjective,
                            "depthPredicted" -> depthPredicted,
                            "wellValues" -> values,
                            "result" -> result,
                            "plotData" -> table,
                            "lmSet" -> lmSet,
                            "lmParametres" -> lmParametres,
                            "fits" -> fits,
                            "errors" -> errors,
                            "RMSError" -> 4,
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


THMethod3D[ds_, timeMaps_, var_]:= 
Module[{
                i,
                j,
                xy,
                values,
                lmSet,
                fits,
                lmParametres,
                depthObjective,
                depthPredicted,
                errors,
                interpolationErrors,
                result                
},               
				               
               
	            values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "depth", "time"}]]]], {i, Normal[DeleteDuplicates[ds[[All, "horizon"]]]]}];   
				values[[All, All, 3]] /= 2000;  (*divide time*)
                If[values[[1, 1, 2]] < 0, values[[All, All, 2]] *= -1];  (*transform depths to positive values for easy calculations and presentable plots*)           
                                          
                lmSet = Table[LinearModelFit[values[[i]][[All, 2 ;; 3]], var, var], {i, Length[values]}]; (*evaluate linear models v(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)

                fits = Table[Table[{timeMaps[[i, j, 1]], timeMaps[[i, j, 2]], -(timeMaps[[i, j, 3]] / lmParametres[[i, 2]] / 2000 - lmParametres[[i, 1]] / lmParametres[[i, 2]])}, {j, Length[timeMaps[[i]]]}], {i, Length[timeMaps]}]; (*evaluate fits (x, h = t/a - b/a). make data suitable for ListLinePlot*)
                        
				xy = Table[Values[Normal[ds][[All, {"x", "y"}]]], {i, Length[DeleteDuplicates[ds[[All, "horizon"]]]]}]; (*x y of wells*)
                
                (*find errors between wells values and predicted by model*)
				
                If[values[[1, 1, 2]] < 0, depthObjective = values[[All, All, {1, 2}]], depthObjective = Table[Map[{#[[1]], -#[[2]]}&, values[[i]]], {i, Length[values]}]];
                depthPredicted = Table[Table[{values[[i, j, 1]], -(values[[i, j, 3]] / lmParametres[[i, 2]] - lmParametres[[i, 1]] / lmParametres[[i, 2]])}, {j, Length[values[[i]]]}], {i, Length[values]}];
				errors = Table[Table[{values[[i, j, 1]], xy[[i, j, 1 ;; 2]], (depthObjective - depthPredicted)[[i, j, 2]]}, {j,  Length[values[[i]]]}], {i, Length[values]}];	
                
                (*interpolatioan and adding error*)

              
                Return[<|"depthObjective" -> depthObjective,
                            "depthPredicted" -> depthPredicted,
                            "wellValues" -> values,
                            "fits"-> fits,
                            "result" -> result,     
                            "lmSet" -> lmSet,
                            "lmParametres" -> lmParametres,
                            "errors" -> errors[[All, 2]],
                            "RMSError"->Map[StandardDeviation[#[[All, 3]]]&, errors],
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


VaveMethod3D[ds_,  timeMaps_]:= 
Module[{                
                i,
                j,
                xy,
                values,                
                errors,
                interpolationErrors,
                result,
                vAveTable,
                fits,
                depthObjective,
                depthPredicted
},               			
                
                values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "depth", "time"}]]]], {i, Normal[DeleteDuplicates[ds[[All, "horizon"]]]]}];  
				values[[All, All, 3]] /= 2000;  (*divide time*)
                If[values[[1, 1, 2]] < 0, values[[All, All, 2]] *= -1];  (*transform depths to positive values for easy calculations and presentable plots*)           
                                 
                vAveTable = Table[Mean[Table[(values[[i, j, 2]])/values[[i, j, 3]], {j, Length[values[[i]]]}]], {i, Length[values]}]; (*define average velocity for each horizon for each position on section*)
                fits = Table[Table[{timeMaps[[i, j, 1]], timeMaps[[i, j, 2]], -vAveTable[[i]] * timeMaps[[i, j, 3]]/2000}, {j, Length[timeMaps[[i]]]}], {i, Length[timeMaps]}]; (*evaluate fits (x, h = vAve*t/2). make data suitable for ListLinePlot*)
				
               xy = Table[Values[Normal[ds][[All, {"x", "y"}]]], {i,  Length[DeleteDuplicates[ds[[All,"horizon"]]]]}]; (*x y of wells*)  (**drop Table! everywhere!!!!!!!!**)
				
				(*find errors between wells values and predicted by model*)

				If[values[[1, 1, 2]] < 0, depthObjective = values[[All, All, {1, 2}]], depthObjective = Table[Map[{#[[1]], -#[[2]]}&, values[[i]]], {i, Length[values]}]];
                depthPredicted = Table[Table[{values[[i, j, 1]], -vAveTable[[i]] * values[[i, j, 3]]} , {j, Length[values[[i]]]}], {i, Length[vAveTable]}];
				errors = Table[Table[{values[[i,j, 1]], xy[[i, j, 1;;2]], (depthObjective - depthPredicted)[[i, j, 2]]}, {j,  Length[values[[i]]]}], {i, Length[vAveTable]}];	

                (*interpolatioan and adding error*)
				
				
                Return[<|"depthObjective" -> depthObjective,
                            "depthPredicted" -> depthPredicted,
                            "wellValues" -> values,
                            "fits"-> fits,
                            "result" -> result,
                            "vAveTable" -> vAveTable,
                            "errors" -> errors,
                            "RMSError" -> Map[StandardDeviation[#[[All, 3]]]&, errors]|>]
]    


dHdTMethod3D[ref_, ds_, timeMaps_, var_]:= 
Module[{                
                i,
                j,
                xy,
                values,
                plotData,
                lmSet,                
                lmParametres,
                fits,
                depthObjective,
                depthPredicted,
                errors,
                interpolationErrors,
                result				
},                              
                
                values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "time", "depth"}]]]], {i, Normal[DeleteDuplicates[ds[[All, "horizon"]]]]}];  
                values[[All, All, 2]] /= 2000;   (*divide time*)     
	            
	            plotData = Table[Table[values[[i, j, {2, 3}]] - values[[1, j, {2, 3}]], {j, Length[values[[i]]]}], {i, 2, Length[values]}];
                If[values[[1, 1, 3]] < 0, values[[All, All, 3]] *= -1];   (*transform depths to positive values for easy calculations and presentable plots*) 
                
                lmSet = Table[LinearModelFit[(values[[i]][[All, 2 ;; 3]] - values[[1]][[All, 2 ;; 3]]), var, var], {i,  2, Length[values]}];	(*evaluate linear models h(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
	  
	            fits  = Table[Table[{timeMaps[[i, j, 1]], timeMaps[[i, j, 2]], ref[[j, 3]] - (lmParametres[[i, 2]] * (timeMaps[[i + 1, j, 3]] - timeMaps[[1, j, 3]])/2000 + lmParametres[[i, 1]])}, {j, Length[timeMaps[[i]]]}], {i, Length[lmSet]}];	           

				xy = Table[Values[Normal[ds][[All, {"x", "y"}]]], {i,  Length[DeleteDuplicates[ds[[All,"horizon"]]]]}]; (*x y of wells*)

                (*find errors between wells values and predicted by model*)
				
				If[values[[1, 1, 3]] < 0, depthObjective = values[[2 ;; -1, All, {1, 3}]], depthObjective = Table[Map[{#[[1]], -#[[3]]}&, values[[i + 1]]], {i, Length[values] - 1}]];
                depthPredicted = Table[Table[{values[[1, j, 1]], -values[[1, j, 3]] - (lmParametres[[i, 2]] * (values[[i + 1, j, 2]] - values[[1, j, 2]]) + lmParametres[[i, 1]])}, {j, Length[values[[i]]]}], {i, Length[values] - 1}];
				errors = Table[Table[{values[[1, j, 1]], xy[[1, j, 1 ;; 2]], (depthObjective[[i, j, 2]] - depthPredicted[[i, j, 2]])}, {j,  Length[values[[i]]]}], {i, Length[values] - 1}];			  
				
                
                (***********interpolatioan and adding errors********)
                
                (************************************************)

                Return[<|"depthObjective" -> depthObjective,
                            "depthPredicted" -> depthPredicted,
                            "wellValues" -> values,
							"plotData" -> plotData,
                            "result" -> result,
                            "lmSet" -> lmSet,
                            "lmParametres" -> lmParametres,
                            "fits" -> fits,
                            "errors" -> errors,
                            "RMSError" -> Map[StandardDeviation[#[[All, 3]]]&, errors],
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


dVdTMethod3D[ref_, ds_, timeMaps_, var_]:= 
Module[{
                i,
                j,
                xy,                
                plotData,                
                lmParametres,
                values,
                lmSet,
                fits,
                depthObjective,                
                errors,
                interpolationErrors,
                depthPredicted,
                result                  
},               
	          
	            values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "time", "depth"}]]]], {i, Normal[DeleteDuplicates[ds[[All, "horizon"]]]]}];            
				values[[All, All, 2]] /= 2000; (*divide time*)
	            
	            If[values[[1, 1, 3]] < 0, values[[All, All, 3]] *= -1];   (*transform depths to positive values for easy calculations and presentable plots*) 
                
                plotData = Table[Table[{values[[i + 1, j, 2]] - values[[1, j, 2]], (values[[i + 1, j, 3]] - values[[1, j, 3]])/(values[[i + 1, j, 2]] - values[[1, j, 2]])}, {j, Length[values[[i]]]}], {i, Length[values] - 1}]; (*{time, velocity}*)
           
		        lmSet = Table[LinearModelFit[plotData[[i]], var, var], {i, Length[plotData]}]; (*evaluate linear models v(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
                
                fits = Table[Table[{timeMaps[[i, j, 1]], timeMaps[[i, j, 2]], ref[[j, 3]] - (lmParametres[[i, 2]] * ((timeMaps[[i + 1, j, 3]] - timeMaps[[1, j, 3]])/2000)^2  + lmParametres[[i, 1]] * (timeMaps[[i + 1, j, 3]] - timeMaps[[1, j, 3]])/2000)}, {j, Length[timeMaps[[i]]]}], {i, Length[lmSet]}];(*evaluate fits (x, h = at^2+bt). make data suitable for ListLinePlot *)
                        
				xy = Table[Values[Normal[ds][[All,{ "x", "y"}]]], {i,  Length[DeleteDuplicates[ds[[All, "horizon"]]]]}]; (*x y of wells*)

                (*find errors between wells values and predicted by model*)
				
				If[values[[1, 1, 3]] < 0, depthObjective = values[[2 ;; -1, All, {1, 3}]], depthObjective = Table[Map[{#[[1]], -#[[3]]}&, values[[i + 1]]], {i, Length[values] - 1}]];
                depthPredicted = Table[Table[{values[[1, j, 1]], -values[[1, j, 3]] - (lmParametres[[i, 2]] * (values[[i + 1, j, 2]]- values[[1, j, 2]])^2 + lmParametres[[i, 1]] * (values[[i + 1, j, 2]]- values[[1, j, 2]]))}, {j, Length[values[[i]]]}], {i, Length[values]-1}];
				errors = Table[Table[{values[[1, j, 1]], xy[[1, j, 1 ;; 2]], (depthObjective[[i, j, 2]] - depthPredicted[[i, j, 2]])}, {j,  Length[values[[i]]]}], {i, Length[values] - 1}];	
                
                (*interpolatioan and adding error*)


                Return[<|"depthObjective" -> depthObjective,
                            "depthPredicted" -> depthPredicted,
                            "wellValues" -> values,
                            "result" -> result,
                            "plotData" -> plotData,
                            "lmSet" -> lmSet,
                            "lmParametres" -> lmParametres,
                            "fits" -> fits,
                            "errors" -> errors,
                            "RMSError" -> Map[StandardDeviation[#[[All, 3]]]&, errors],
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


dTdHMethod3D[ref_, ds_, timeMaps_, var_]:= 
Module[{
                i,
                j,
                xy,
                values,
                plotData,
                lmSet,
                lmParametres,
                depthObjective,
                depthPredicted,                
                fits,
                errors,
                interpolationErrors,
                result	              
},               
				values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "depth", "time"}]]]], {i, Normal[DeleteDuplicates[ds[[All, "horizon"]]]]}];  
				values[[All, All, 3]] /= 2000;  (*divide time*)
	            plotData = Table[Table[values[[i, j, {2, 3}]] - values[[1, j, {2, 3}]], {j, Length[values[[i]]]}], {i, 2, Length[values]}];

                If[values[[1, 1, 2]] < 0, values[[All, All, 2]] *= -1];  (*transform depths to positive values for easy calculations and presentable plots*)           
                                          
                lmSet = Table[LinearModelFit[values[[i]][[All, 2 ;; 3]] - values[[1]][[All, 2 ;; 3]], var, var], {i, 2, Length[values]}]; (*evaluate linear models v(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)

                fits = Table[Table[{timeMaps[[i, j, 1]], timeMaps[[i, j, 2]], ref[[j, 3]]-((timeMaps[[i + 1, j, 3]] - timeMaps[[1, j, 3]])/ lmParametres[[i, 2]]/2000 - lmParametres[[i, 1]]/lmParametres[[i, 2]])}, {j, Length[timeMaps[[i]]]}], {i, Length[lmSet]}]; (*evaluate fits (x, h = t/a - b/a). make data suitable for ListLinePlot*)
                        
				xy = Table[Values[Normal[ds][[All,{ "x", "y"}]]], {i,  Length[DeleteDuplicates[ds[[All, "horizon"]]]]}]; (*x y of wells*)
                
                (*find errors between wells values and predicted by model*)
				
                If[values[[1, 1, 2]] < 0, depthObjective = values[[2 ;; -1, All, {1, 2}]], depthObjective = Table[Map[{#[[1]], -#[[2]]}&, values[[i + 1]]], {i, Length[values] - 1}]];
                depthPredicted = Table[Table[{values[[i, j, 1]], -values[[1, j, 2]] - ((values[[i + 1, j, 3]] - values[[1, j, 3]]) / lmParametres[[i, 2]] - lmParametres[[i, 1]] / lmParametres[[i, 2]])}, {j, Length[values[[i]]]}], {i, Length[values] - 1}];
				errors = Table[Table[{values[[1, j, 1]], xy[[1, j, 1 ;; 2]], (depthObjective[[i, j, 2]] - depthPredicted[[i, j, 2]])}, {j,  Length[values[[i]]]}], {i, Length[values] - 1}];	
                
                (*interpolatioan and adding error*)

              
                Return[<|"depthObjective" -> depthObjective,
                            "depthPredicted" -> depthPredicted,
                            "wellValues" -> values,
							"plotData" -> plotData,
                            "fits"-> fits,
                            "result" -> result,     
                            "lmSet" -> lmSet,
                            "lmParametres" -> lmParametres,
                            "errors" -> errors,
                            "RMSError" -> Map[StandardDeviation[#[[All, 3]]]&, errors],
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


dVaveMethod3D[ref_, ds_, timeMaps_]:= 
Module[{                
                i,
                j,
                xy,
                values,                
                errors,
                interpolationErrors,
                result,
                vAveTable,
                fits,
                depthObjective,
                depthPredicted
},               			
                
                values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "depth", "time"}]]]], {i, Normal[DeleteDuplicates[ds[[All, "horizon"]]]]}];  
                values[[All, All, 3]] /= 2000;  (*divide time*)
                If[values[[1, 1, 2]] < 0, values[[All, All, 2]] *= -1];  (*transform depths to positive values for easy calculations and presentable plots*)           
                                 
                vAveTable = Table[Mean[Table[(values[[i + 1, j, 2]] - values[[1, j, 2]])/(values[[i + 1, j, 3]] - values[[1, j, 3]]), {j, Length[values[[i]]]}]], {i, Length[values] - 1}]; (*define average velocity for each horizon for each position on section*)
                fits = Table[Table[{timeMaps[[i, j, 1]], timeMaps[[i, j, 2]], ref[[j, 3]] - vAveTable[[i]] * (timeMaps[[i + 1, j, 3]]  -  timeMaps[[1,j, 3]])/ 2000}, {j, Length[timeMaps[[i]]]}], {i, Length[timeMaps] - 1}]; (*evaluate fits (x, h = vAve*t/2). make data suitable for ListLinePlot*)
				
                xy = Table[Values[Normal[ds][[All,{ "x", "y"}]]], {i,  Length[DeleteDuplicates[ds[[All, "horizon"]]]]}]; (*x y of wells*)
               
				(*find errors between wells values and predicted by model*)

				If[values[[1, 1, 2]] < 0, depthObjective = values[[2 ;; -1, All, {1, 2}]], depthObjective =Table[Map[{#[[1]], -#[[2]]}&, values[[i + 1]]], {i, Length[values] - 1}]];
                depthPredicted = Table[Table[{values[[1, j, 1]], -values[[1, j, 2]] - vAveTable[[i]] * (values[[i + 1, j, 3]]- values[[1, j, 3]])}, {j, Length[values[[i]]]}], {i, Length[values] - 1}];
				errors = Table[Table[{values[[1, j, 1]], xy[[1, j, 1 ;; 2]], (depthObjective[[i, j, 2]] - depthPredicted[[i, j, 2]])}, {j,  Length[values[[i]]]}], {i, Length[values] - 1}];	

                (*interpolatioan and adding error*)
				
				
                Return[<|"depthObjective" -> depthObjective,
                            "depthPredicted" -> depthPredicted,
                            "wellValues" -> values,
                            "fits"-> fits,
                            "result" -> result,
                            "vAveTable" -> vAveTable,
                            "errors" -> errors,
                            "RMSError" -> Map[StandardDeviation[#[[All, 3]]]&, errors]|>]
]    


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
