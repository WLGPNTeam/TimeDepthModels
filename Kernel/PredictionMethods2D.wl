(* ::Package:: *)

(* ::Chapter:: *)
(*PredictionMethods2D*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[HTMethod2D, VTMethod2D, THMethod2D, VaveMethod2D, dHdTMethod2D, dVdTMethod2D, dTdHMethod2D, dVaveMethod2D, VaveVOGTMethod2D]


HTMethod2D::usage = 
"HTMethod2D[ds, timeSection, var]"


VTMethod2D::usage = 
"VTMethod2D[ds, timeSection, var]"


THMethod2D::usage = 
"THMethod2D[ds, timeSection, var]"


VaveMethod2D::usage = 
"VaveMethod2D[ds, timeSection]"


dHdTMethod2D::usage = 
"dHdTMethod2D[ref, ds, timeSection, var]"


dVdTMethod2D::usage = 
"dVdTMethod2D[ref, ds, timeSection, var]"


dTdHMethod2D::usage = 
"dTdHMethod2D[ref, ds, timeSection, var]"


dVaveMethod2D::usage = 
"dVaveMethod2D[ref, ds, timeSection]"


VaveRegressionMethod2D::usage = 
"VaveRegressionMethod2D[ds, horizons, timeSection, velSection, var]"


(* ::Section::Closed:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


HTMethod2D[ds_, timeSection_, var_]:= 
Module[{
                
                i,
                j,
                dx,
                len,
                surface,
                positions,
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
				dx = timeSection[[1, 2, 1]]; (*section step*)
	            len = timeSection[[1, -1, 1]]; (*length of section*)
	            
                values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "time", "depth"}]]]], {i,  DeleteDuplicates[Normal[ds[All, "horizon"]]][[2 ;; -1]]}]; (*time, depth*)
                values[[All, All, 2]] /= 2;   (*divide time*)     
                values[[All, All, 3]] *= -1;   (*transform depths to positive values for easy calculations and presentable plots*) 
                
                lmSet = Table[LinearModelFit[values[[i]][[All, 2 ;; 3]], var, var], {i,  Length[values]}];	(*evaluate linear models h(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
                fits = Table[Table[{(j - 1) dx, -(lmParametres[[i, 2]] * timeSection[[i + 1, j, 2]] / 2 + lmParametres[[i, 1]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*make data suitable for ListLinePlot. (x, h = at+b)*)
								
					(*find errors between wells values and predicted by model*)
					positions = DeleteDuplicates[Normal[ds[All, "x"]]]; (*x positions of wells*)
					depthObjective = -values[[All, All, 3]];
					depthPredicted = Table[Flatten[Cases[fits[[i]], {#, __}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
					errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, Length[positions]}], {i, Length[depthObjective]}];
					interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
					
					result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
					
				
                Return[<|"depthObjective" -> depthObjective,
                "depthPredicted" -> depthPredicted,
                "wellValues" -> values,
                "result" -> result,
                "lmSet" -> lmSet,
                "lmParametres" -> lmParametres,
                "fits" -> fits,
                "errors" -> errors[[All, All,2]],
                "RMSError" -> Map[StandardDeviation[#[[All, 2]]]&, errors],
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


dHdTMethod2D[subset_, ref_, ds_, timeSection_, var_]:= 
Module[{
                i,
                j,
                dx,
                len,
                values,
                lmSet ,
                lmParametres,
                fits,
                positions,
                depthObjective,
                depthPredicted,
                errors,
                interpolationErrors,
                result,
                plotData
},               
		
            dx = timeSection[[1, 2, 1]]; (*x step on section*)
            len = timeSection[[1, -1, 1]]; (*length of section*)
            positions = DeleteDuplicates[Normal[ds[All, "x"]]]; (*x positions of wells*)
            
            values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"time", "depth"}]]]], {i, subset}] ;(*time, depth*) 
            values[[All, All, 1]] /= 2; (*divide time*)
            values[[All, All, 2]] *= -1; (*transform depths to positive values for easy calculations and presentable plots*)    

            plotData = Table[values[[i]] - values[[1]], {i, 2, Length[values]}];
            
            lmSet = Table[LinearModelFit[values[[i]] - values[[1]], var, var], {i, 2, Length[values]}]; (*evaluate linear models dh(dt) = adt + b for each horizon*) 
            lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
            fits = Table[Table[{(j - 1) dx, ref[[j, 2]] - (lmParametres[[i, 2]] * (timeSection[[subset[[i + 1]], j, 2]]  - timeSection[[subset[[1]], j, 2]])/ 2 + lmParametres[[i, 1]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*evaluate fits (x, h = href + dt/a - b/a). make data suitable for ListLinePlot*)

	        (*find errors between wells values and predicted by model*)
            depthObjective = -values[[2 ;; -1, All, 2]];
            depthPredicted = Table[Flatten[Cases[fits[[i]], {#, __}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
            errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}], {i, Length[depthObjective]}];
            interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)


            result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
				
				
            Return[<|"depthObjective" -> depthObjective,
            "depthPredicted" -> depthPredicted,
            "wellValues" -> values,
            "plotData" -> plotData,
            "result" -> result,
            "lmSet" -> lmSet,
            "lmParametres" -> lmParametres,
            "fits" -> fits,
            "errors" -> errors[[All, All, 2]],
            "RMSError" -> Map[StandardDeviation[#[[All, 2]]]&, errors],
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]			
         
]


VTMethod2D[ds_, timeSection_, var_]:= 
Module[{
               
                i,
                j,
                x,
                h,
                v,
                dx,
                len,
                table,                
                lmParametres,
                values,
                lmSet,
                fits,
                depthObjective,                
                positions,
                errors,
                interpolationErrors,
                depthPredicted,
                result                
},               
	            dx = timeSection[[1, 2, 1]]; (*section step*)
	            len = timeSection[[1, -1, 1]]; (*length of section*)	  
	                      
	            values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "time", "depth"}]]]], {i,  DeleteDuplicates[Normal[ds[All, "horizon"]]][[2 ;; -1]]}]; (*well, time, depth*)                
	            values[[All, All, 2]] /= 2; (*divide time*)
	            values[[All, All, 3]] *= -1; (*transform depths to positive values for easy calculations and presentable plots*) 
                
                table = Table[Table[{i, values[[i, j, 2]], values[[i, j, 3]] / values[[i, j, 2]]}, {j, Length[values[[i]]]}], {i, Length[values]}]; (*{horizon, time, vel}*)                
 
		        lmSet = Table[LinearModelFit[table[[i]][[All, 2 ;; 3]], var, var], {i, Length[table]}]; (*evaluate linear models v(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
                fits = Table[Table[{(j - 1) dx, -1/2*(lmParametres[[i, 2]] * timeSection[[i + 1, j, 2]]^2/2  + lmParametres[[i, 1]] * timeSection[[i + 1, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}];(*evaluate fits (x, h = at^2+bt). make data suitable for ListLinePlot *)
                
	
				(*find errors between wells values and predicted by model*)
				positions = DeleteDuplicates[Normal[ds[All, "x"]]]; (*x positions of wells*)
				depthObjective = -values[[All, All, 3]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#,__}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, Length[positions]}], {i, Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}];(*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
								
                Return[<|"depthObjective" -> depthObjective,
                "depthPredicted" -> depthPredicted,
                "wellValues" -> table,
                "lmSet" -> lmSet,
                "lmParametres" -> lmParametres,
                "result" -> result,
                "fits" -> fits, 
                "errors" -> errors[[All,All, 2]],
                "RMSError" -> Map[StandardDeviation[#[[All,2]]]&, errors],
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


dVdTMethod2D[subset_, ref_, ds_, timeSection_, var_]:= 
Module[{
               
                i,
                j,                
                dx,
                len,
                table,                
                lmParametres,
                values,
                lmSet,
                fits,                
                depthObjective,                
                positions,
                errors,
                interpolationErrors,
                depthPredicted,
                result	                           
},               
	            dx = timeSection[[1, 2, 1]]; (*x step on section*)
	            len = timeSection[[1, -1, 1]]; (*length of section*)
	            
	                     
	            values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"x", "time", "depth"}]]]], {i, subset}]; (*x, time, depth*)
	            values[[All, All, 2]] /= 2; (*divide time*)
	            values[[All, All, 3]] *= -1; (*transform depths to positive values for easy calculations and presentable plots*)    

		        positions = DeleteDuplicates[Normal[ds[All, "x"]]]; (*x positions of wells*)
		        
                table = Table[Table[{i, values[[i + 1, j, 2]] - values[[1, j, 2]], (values[[i + 1, j, 3]] - values[[1, j, 3]])/(values[[i + 1, j, 2]] - values[[1, j, 2]]) }, {j, Length[values[[i]]]}], {i, Length[values] - 1}]; (*{horizon, dt[i,j] = t[i,j] - tref[i,j] , 0}*)
               
 
		        lmSet = Table[LinearModelFit[table[[i]][[All, 2 ;; 3]], var, var], {i, Length[table]}]; (*evaluate linear models v(dt) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
                fits = Table[Table[{(j - 1) dx, ref[[j, 2]] - 1/2 * (lmParametres[[i, 2]] * (timeSection[[subset[[i + 1]], j, 2]] - timeSection[[subset[[1]], j, 2]])^2/2  + lmParametres[[i, 1]] * (timeSection[[subset[[i + 1]], j, 2]] - timeSection[[subset[[1]], j, 2]]))}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*evaluate fits (x, h = href + adt^2 + bdt). make data suitable for ListLinePlot *)
                 
	
				(*find errors between wells values and predicted by model*)
				depthObjective = -values[[2 ;; -1, All, 3]]; 
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#, __}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, Length[positions]}], {i, Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
				
				
                Return[<|"depthObjective" -> depthObjective,
                "positions" -> positions,
                "depthPredicted" -> depthPredicted,
                "wellValues" -> values,                
                "plotData" -> table,
                "lmSet" -> lmSet,
                "lmParametres" -> lmParametres,
                "result" -> result,
                "fits" -> fits,
                "errors" -> errors[[All, All, 2]],
                "RMSError" -> Map[StandardDeviation[#[[All,2]]]&, errors],
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


THMethod2D[ds_, timeSection_, var_]:= 
Module[{
                i,
                j,
                dx,
                len,
                values,
                lmSet,
                fits,
                lmParametres,
                positions,
                depthObjective,
                depthPredicted,
                errors,
                interpolationErrors,
                result                
},               
				dx = timeSection[[1, 2, 1]]; (*x step on section*)
	            len = timeSection[[1, -1, 1]]; (*length of section*)
                
                
	            values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "depth", "time"}]]]], {i,  DeleteDuplicates[Normal[ds[All, "horizon"]]][[2 ;; -1]]}]; (*h, t*)
                values[[All, All, 3]] /= 2;  (*divide time*)
                values[[All, All, 2]] *= -1;  (*transform depths to positive values for easy calculations and presentable plots*)           
                
                lmSet = Table[LinearModelFit[values[[i]][[All, 2 ;; 3]], var, var], {i, Length[values]}]; (*evaluate linear models v(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
                fits = Table[Table[{(j - 1) dx, -(timeSection[[i + 1, j, 2]] / lmParametres[[i, 2]] / 2 - lmParametres[[i, 1]] / lmParametres[[i, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*evaluate fits (x, h = t/a - b/a). make data suitable for ListLinePlot*)
								
				(*find errors between wells values and predicted by model*)
				positions = DeleteDuplicates[Normal[ds[All, "x"]]]; (*x positions of wells*)
				depthObjective = -values[[All, All, 2]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#, __}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, Length[positions]}],{i,  Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
				
                Return[<|"depthObjective" -> depthObjective,
                "depthPredicted" -> depthPredicted,
                "wellValues" -> values,
                "fits"-> fits,
                "result" -> result,
                "lmSet" -> lmSet,
                "lmParametres" -> lmParametres,
                "errors" -> errors[[All, All,2]],
                "RMSError" -> Map[StandardDeviation[#[[All, 2]]]&, errors],
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


dTdHMethod2D[subset_, ref_, ds_, timeSection_, var_]:= 
Module[{
                i,
                j,
                dx,
                len,
                values,
                lmSet,
                fits,
                lmParametres,
                positions,
                depthObjective,
                depthPredicted,
                errors,
                interpolationErrors,
                result,
                plotData             
},               
				dx = timeSection[[1, 2, 1]]; (*x step on section*)
	            len = timeSection[[1, -1, 1]]; (*length of section*)
                
                values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {4, 5}]]]], {i, subset}]; (*h, t*)
                values[[All, All, 2]] /= 2;   (*divide time*)
                values[[All, All, 1]] *= -1;   (*transform depths to positive values for easy calculations and presentable plots*)             
                
                positions = DeleteDuplicates[Normal[ds[All, "x"]]]; (*x positions of wells*)
                
                plotData = Table[values[[i]] - values[[1]], {i, 2, Length[values]}];
                
                lmSet = Table[LinearModelFit[values[[i]] - values[[1]], var, var], {i, 2, Length[values]}];	(*evaluate linear models dt(dh) = adh + b for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)                
                fits = Table[Table[{(j - 1) dx,  ref[[j, 2]] - ((timeSection[[subset[[i + 1]], j, 2]] - timeSection[[subset[[1]], j, 2]])/ lmParametres[[i, 2]] / 2 - lmParametres[[i, 1]] / lmParametres[[i, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*evaluate fits (x, h = href + dt/a - b/a). make data suitable for ListLinePlot*)
								
				(*find errors between wells values and predicted by model*)
				depthObjective = -values[[2 ;; -1, All, 1]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#, __}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}],{i,  Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
								
                Return[<|"depthObjective" -> depthObjective,
                "depthPredicted" -> depthPredicted,
                "wellValues" -> values,
                "plotData" -> plotData,
                "result" -> result,
                "lmSet" -> lmSet,
                "lmParametres" -> lmParametres,
                "fits" -> fits,
                "errors" -> errors[[All, All,2]],
                "RMSError" -> Map[StandardDeviation[#[[All, 2]]]&, errors],
                 "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


VaveMethod2D[ds_,  timeSection_]:= 
Module[{`                
                i,
                j,
                dx,
                len,
                values,
                positions,
                errors,
                interpolationErrors,
                result,
                vAveTable,
                fits,                
                depthObjective,
                depthPredicted
},               
				dx = timeSection[[1, 2, 1]]; (*x step on section*)
	            len = timeSection[[1, -1, 1]]; (*length of section*)
                                
	            values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "depth", "time"}]]]], {i,  DeleteDuplicates[Normal[ds[All, "horizon"]]][[2 ;; -1]]}]; (*h, t*)
                values[[All, All, 3]] /= 2;   (*divide time*)
                values[[All, All, 2]] *= -1;   (*transform depths to positive values for easy calculations and presentable plots*)     
                
                vAveTable = Table[Mean[Table[(values[[i, j, 2]])/values[[i, j, 3]], {j, Length[values[[i]]]}]], {i, Length[values]}]; (*define average velocity for each horizon for each position on section*)
                fits = Table[Table[{(j - 1) dx, -vAveTable[[i]] * timeSection[[i + 1, j, 2]] / 2}, {j, len/dx + 1}], {i, Length[vAveTable]}]; (*evaluate fits (x, h = vAve*t/2). make data suitable for ListLinePlot*)
								
				(*find errors between wells values and predicted by model*)
				positions = DeleteDuplicates[Normal[ds[All, "x"]]]; (*x positions of wells*)
				depthObjective = -values[[All, All, 2]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#, __}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, Length[positions]}], {i, Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[fits]}]; (*add errors to fits*)
				
                Return[<|"depthObjective" -> depthObjective,
                "depthPredicted" -> depthPredicted,
                "wellValues" -> values,
                "fits"-> fits,
                "result" -> result,
                "vAveTable" -> vAveTable,
                "errors" -> errors[[All, All,2]],
                "RMSError" -> Map[StandardDeviation[#[[All, 2]]]&, errors]|>]
]


dVaveMethod2D[subset_, ref_, ds_, timeSection_]:= 
Module[{                
                i,
                j,
                dx,
                len,
                values,
                positions,
                depthObjective,
                errors,
                interpolationErrors,
                result,
                vAveTable,
                fits,
                depthPredicted,
	            referenceValuesInterpolated
},               
				dx = timeSection[[1, 2, 1]]; (*x step on section*)
	            len = timeSection[[1, -1, 1]]; (*length of section*)
                
                values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"depth", "time"}]]]], {i, subset}]; (*pairs (depth, time)*)
                values[[All, All, 2]] /= 2;  (*divide time*)
                values[[All, All, 1]] *= -1;  (*transform depths to positive values for easy calculations and presentable plots*)        
                
                positions = DeleteDuplicates[Normal[ds[All, "x"]]]; (*x positions of wells*)
                
                vAveTable = Table[Mean[Table[Abs[values[[i, j, 1]] - values[[1, j, 1]]]/(values[[i, j, 2]] - values[[1, j, 2]]), {j, Length[values[[i]]]}]], {i, 2, Length[values] }]; (*define average velocity for each horizon for each position on section. wellValues[[1]] - well values of reference horizon *)
                
                fits = Table[Table[{(j - 1) dx, ref[[j, 2]] - vAveTable[[i]] * (timeSection[[subset[[i + 1]], j, 2]] - timeSection[[subset[[1]], j, 2]])/ 2}, {j, len/dx + 1}], {i, Length[vAveTable]}]; (*evaluate fits (x, h = href + vAve*dt/2). make data suitable for ListLinePlot*)
								 
				(*find errors between wells values and predicted by model*)
				depthObjective = -values[[2 ;; -1, All, 1]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#, __}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, Length[positions]}], {i, Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[fits]}]; (*add errors to fits*)
								
                Return[<|"depthObjective" -> depthObjective,
                "depthPredicted" -> depthPredicted,
                "fits" -> fits,
                "result" -> result,
                "vAveTable" -> vAveTable,
                "errors" -> errors[[All,All, 2]],
                "RMSError" -> Map[StandardDeviation[#[[All, 2]]]&, errors]|>]
]


VaveRegressionMethod2D[ds_, horizons_, timeSection_, velSection_, var_]:= 
Module[{                
                i,
                j,
                dx,
                len,
                values,
                positions,                
                errors,
                interpolationErrors,
                result,
                vAveTable,
				vAveTableFullHor,
				table,
				lmSet,
				lmParametres,
                fits,
                depthObjective,
                depthPredicted,
                vOGTTable
},               
				dx = timeSection[[1, 2, 1]]; (*x step on section*)
	            len = timeSection[[1, -1, 1]]; (*length of section*)
                
                values = Table[Values[Normal[ds[Select[#horizon == i&]][[All, {"well", "depth", "time"}]]]], {i,  DeleteDuplicates[Normal[ds[All, "horizon"]]][[2 ;; -1]]}]; (*h, t*)                
                values[[All, All, 2]] *= -1;   (*transform depths to positive values for easy calculations and presentable plots*)     
                
                positions = DeleteDuplicates[Normal[ds[All, "x"]]];
                
                vAveTableFullHor = Table[Table[Abs[(horizons[[i, j, 2]] - horizons[[1, j, 2]]) / timeSection[[i, j, 2]]], {j, Length[horizons[[i]]]}], {i, 2, Length[horizons]}];
				vAveTable = Extract[#, Partition[positions/dx, 1]]&/@vAveTableFullHor; (*define average velocity for each horizon for each position on section*)
				vOGTTable = RandomReal[{0.5, 1.5}] * #&/@vAveTable;(*WILL BE CHANGED = vOGT*)
				table = Table[Table[{vOGTTable[[i, j]], vAveTable[[i, j]]}, {j, Length[vAveTable[[i]]]}], {i, Length[vAveTable]}];
				
				lmSet = Table[LinearModelFit[table[[i]], var, var], {i, Length[table]}];	(*evaluate linear model for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)  

                 fits = Table[Table[{(j - 1) dx, -(vAveTableFullHor[[i, j]] * lmParametres[[i, 2]] + lmParametres[[i, 1]]) * timeSection[[i + 1, j, 2]] / 2}, {j, len/dx + 1}], {i, Length[vAveTableFullHor]}];(*evaluate fits (x, h = vAve*t/2). make data suitable for ListLinePlot*)
								
				(*find errors between wells values and predicted by model*)
				 (*x positions of wells*)
				depthObjective = -values[[All,All,2]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#, __}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, Length[positions]}],{i, Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[fits]}]; (*add errors to fits*)
				
                Return[<|"depthObjective" -> depthObjective,
                "depthPredicted" -> depthPredicted,
                "wellValues" -> values,
                "result" -> result,
                "lmSet" -> lmSet,
                "lmParametres" -> lmParametres,
                "plotData" -> table,
                "fits" -> fits,
                "errors" -> errors[[All, All,2]],
                "RMSError" -> Map[StandardDeviation[#[[All, All,2]]]&, errors],
                "RSquared" -> Table[lmSet[[i]]["RSquared"], {i, Length[lmSet]}]|>]
]


(*VOGTMethod[wellDataset_, horizons_,  time_, vOGT_, v_]:= 
Module[{                
                i,
                j,
                dx,
                len,
                wellValues,
                surface,
                positions,                
                errors,
                interpolationErrors,
                result,
                vAveTable,
				vAveTableFullHor,
				table,
				lmSet,
				lmParametres,
                fits,
                depthObjective,
                depthPredicted
},               
				dx = time[[1, 2, 1]]; (*x step on section*)
	            len = time[[1, -1, 1]]; (*length of section*)
                
                wellValues = Table[Values[Normal[wellDataset[Select[#horizon == i&]][[All, 4]]]], {i, 2, Max[wellDataset[All, 3]]}]; (*depth*)
                wellValues[[All, All, 1]] *= -1;   (*transform depths to positive values for easy calculations and presentable plots*)     
                
                positions = DeleteDuplicates[Normal[wellDataset[All, "x"]]];
                surface = Values[Normal[wellDataset[Select[#horizon == 1&]][[All, {2, 4}]]]]; (*(x, depth) of first horizon*)
             
		    	fits = Table[Table[{(j - 1) dx, -(vOGT[[i(*!*), j]] * time[[i + 1, j, 2]] / 2)}, {j, len/dx + 1}], {i, Length[time] - 1 }];(*evaluate fits (x, h = vAve*t/2). make data suitable for ListLinePlot*)
								
				(*find errors between wells values and predicted by model*)
				(*x positions of wells*)
				depthObjective = -wellValues;
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#, __}]&/@positions, 1], {i, 1, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}],{i, 1, Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[fits]}]; (*add errors to fits*)
				If[Length[surface] != 0, result = Join[{surface}, result]];
				
                Return[<|"depthObjective" -> depthObjective,
                "depthPredicted" -> depthPredicted,
                "fits" -> fits,
                "result" -> result,
                "vAveTable" -> vAveTable,
                "fits" -> fits,
                "errors" -> errors[[All, 2]],
                "RMSError" -> Map[StandardDeviation[#[[All, 2]]]&, errors]|>]
]*)


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
