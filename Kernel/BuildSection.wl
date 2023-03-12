(* ::Package:: *)

(* ::Chapter:: *)
(*BuildSection*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[BuildDepthSection, BuildVelocitySection, BuildTimeSection, BuildWellSet, CheckVelocity, HTMethod, VTMethod, THMethod, VaveMethod, dHdTMethod, dVdTMethod, dTdHMethod, dVaveMethod, AllMethods]


BuildDepthSection::usage = 
"BuildDepthSection[listH, slopes, len, dx, type, hDispersion, hTapering, gaussRadius]"


BuildVelocitySection::usage =
"BuildVelocitySection[horizons, listV, trends, anomalies]"


BuildTimeSection::usage = 
"BuildTimeSection[horizons, model]"


BuildWellSet::usage = 
"BuildWellSet[horizons, time, wellCount, wellType]"


CheckVelocity::usage = 
"CheckVelocity[wellDataset]"


HTMethod::usage = 
"HTMethod[wellDataset, time, t]"


VTMethod::usage = 
"VTMethod[wellDataset, time, t]"


THMethod::usage = 
"THMethod[wellDataset, time, t]"


VaveMethod::usage = 
"VaveMethod[wellDataset, time]"


dHdTMethod::usage = 
"dHdTMethod[setHT, reference, wellDataset, time, dt]"


dVdTMethod::usage = 
"dVdTMethod[setVT, reference, wellDataset, time, dt]"


dTdHMethod::usage = 
"dTdHMethod[setTH, reference, wellDataset, time, dh]"


dVaveMethod::usage = 
"dVaveMethod[setVave, reference, wellDataset, time]"


AllMethods::usage = 
"AllMethods[set, horizons, wellDataset, time, dt, dh]"


(* ::Section::Closed:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


BuildDepthSection[listH_, slopes_, len_, dx_, type_, hDispersion_, hTapering_, gaussRadius_] := 
Module[{
                i,
                j,
                horizons,
                radius,
                parametres,
                data,
                anomaly,
                surface, 
                anomalyFiltered,
                listSumH,
                table,
                totalH,
                dh1,
                dh2,
                max,
                sums
},
				If[hTapering >= 1/2 ((len/dx) - radius) listH[[1]]/totalH , Return["hTapering is too big"]]; 
 
				radius = gaussRadius; (*GaussianFilter radius for layer with the highest dispersion*)
				parametres = {0.1, 1}; (*parametres for ARProcess*)
                
                totalH = Total[listH]; 
                                
                (*flat layers made*)       
                If[MatchQ[slopes, {} (*Or {0,0,0,...}*)], 
                  table = Join[{Table[0., len/dx + 1]}, Table[-Sum[listH[[i]], {i, 1, i}], {i, Length[listH]}, {j, len/dx + 1}]],
                   (*else*)        
                  table = Table[Table[0, {j, len/dx + 1}], {i, Length[listH] + 1}]; (*add dipping if slope != 0*)
                  sums = Table[Sum[listH[[i]], {i, 1, i}], {i, Length[listH]}];
				  For[i = 1, i <= Length[listH] + 1 , i++, 
					If[i == 1,
						For[j = 1, j <= len/dx + 1, j++, table[[1, j]] = -N[Tan[slopes[[1]]]*(j - 1) * dx]],
						For[j = 1, j <= len/dx + 1, j++, table[[i, j]] = -N[Tan[slopes[[i]]]*(j - 1) * dx]];
						dh1 = Abs[Max[Table[table[[i, j]] - table[[i - 1, j]], {j, len/dx + 1}]]];
						dh2 = sums[[i - 1]];
						If[dh1 >= dh2, max = dh1, max = dh2];
					    table[[i, All]] = Function[x, x - max*1.1]/@table[[i, All]]
                    ] 
                  ]
                ];
                (*make the highest anomaly*)
                If[hDispersion != 0,
                anomaly = Table[0, len/dx + 1]; (*init anomaly*);
                data = RandomFunction[ARProcess[{parametres[[1]]}, parametres[[2]]], {1, len + 1}]; 
                anomaly = hDispersion GaussianFilter[data["SliceData", Range[1, len + 1, dx]], radius][[1]];
                
                listSumH = Join[{1}, Table[Sum[listH[[i]], {i, 1, i}], {i, 1, Length[listH]}]]; (*needed for increasing Gaussian Filter radius in While below*)
                                
                i = 1;
                While[i < Length[table] + 1, 
                            table[[-type i]] += anomaly; (*"type" defines which (last or first) hor has the highest anomaly; add anomaly*)
                               anomalyFiltered = GaussianFilter[anomaly, (radius + hTapering totalH/listSumH[[-i]])];(*changing local anomaly for each horizon*) 
                               max = Abs[Max[Table[type*(anomaly[[j]] - anomalyFiltered[[j]]), {j, 1, Length[anomaly]}]]]; (*finding max difference between previous and current anomaly; need to know so that the inclined do not intersect*)
                               anomaly = Table[anomalyFiltered[[j]] + type * max, {j, 1, Length[anomalyFiltered]}]; (*evaluate anomaly for the next horizon*)
                               i++;
                ]
    ];
                (*make data suitable for ListLinePlot and add first horizon as {0,0...0} if it does not exist*)
                If[MatchQ[table[[1]], Table[0., {j, len/dx + 1}]], 
                horizons = Table[{(j - 1) dx, Round[table[[i, j]]]}, {i, Length[listH] + 1}, {j, len/dx + 1}],
                max = Max[table[[1]]]; (*need to shift horizons under {0,0...0}*)
                surface = Table[{(j - 1) dx, 0}, {j, len/dx + 1}];
                horizons = Join[{surface}, Table[{(j - 1) dx, (Round[table[[i, j]]] - 1.1 * max)}, {i, Length[listH] + 1}, {j, len/dx + 1}]]                  
                  ];
                                    
                Return[<|"horizons" -> horizons, "table" -> table|>]
]


BuildVelocitySection[horizons_, listV_, trendsQ_, anomaliesQ_] := 
Module[{
                i, 
                j, 
                dh,
                dx,
                f,
                maxH,
                vAve,
                trendQ,
                anomalyQ,
                horizonsJoined,
                layerThickness,
                dataset,
                model,
                trend,
                anomaly
                
},
                dx = horizons[[1, 2, 1]] - horizons[[1, 1, 1]]; (*section step*)
                maxH = -Min[horizons[[-1]]]; (*find max depth*)
                
                f[vAve_, trendQ_, anomalyQ_, dh_] := vAve + trendQ * 100 Sqrt[dh] + anomalyQ * RandomInteger[{-vAve/100, vAve/100}];(*vAve - velocity from listV for each layer, trendsQ - depth trend exist or not, anomaliesQ - anomalies exist or not*)
                                              
                model = Table[{0, 0, 0, 0}, {i, 1, Length[horizons], 1}, {j, 1, Length[horizons[[i]]]}, {dh, 1, maxH}]; 
                (*init model [horizons count x pickets count x layer thickness]; use maxH here beccause thickness is different at each (i,j) position; zeros will be deleted after)*)
                horizonsJoined = Join[Part[horizons, All, All, -1], {Table[-(maxH + 1), {j, 1, Length[horizons[[1]]]}]}]; (*modificated horizons*)
               
                For[i = 1, i <= Length[horizons], i++, (*go into layer*)
                    vAve = listV[[i]];(*init vAve for the layer*)
                    trendQ = trendsQ[[i]];(*trendQ = 1 or 0*)
                    anomalyQ = anomaliesQ[[i]]; (*anomalyQ = 1 or 0*)
                    For[j = 1, j <= Length[horizons[[i]]], j++, (*go to j position on section*)
                        layerThickness = Ceiling[Abs[horizonsJoined[[i + 1]][[j]] - horizonsJoined[[i]][[j]]]]; (*define i layer thickness on j position*)
                        For[dh = 0, dh <= layerThickness, dh++, (*go vertical and evaluate velocity on each dh step*)
                            model[[i]][[j]][[dh + 1]] = {i, (j - 1) dx, Round[horizonsJoined[[i, j]]] - dh, Round[f[vAve, trendQ, anomalyQ, dh]]}; (* = {layer, x coordinate, h cootdinate, velocity}*)
                        ]
                    ]
                ];  
                
                model = Table[DeleteCases[model[[i, j]], {0, 0, 0, 0}], {i, 1, Length[horizons]}, {j, 1, Length[horizons[[i]]]}]; (*delete (0, 0, 0, 0) cases*)
                dataset = Dataset[Map[<|"layer" -> #[[1]], "x" -> #[[2]], "h" -> #[[3]], "v" -> #[[4]]|>&, Flatten[model, 2]]]; (*make dataset*)
                
                <|"model" -> model, "dataset" -> dataset|>
]


BuildTimeSection[horizons_, model_] := 
Module[{
                time,
                dx,
                i,
                j,
                vave,
                dh,
				dt
},
                dx = horizons[[1, 2, 1]] - horizons[[1, 1, 1]]; (*section step*)
                time = Table[{0, 0}, {i, 1, Length[horizons]}, {j, 1, Length[horizons[[1]]]}]; (*init time [num of inclined x num of pickets] *)
        
       
                For[i = 1, i <= Length[horizons], i++, (*go into layer*)
                For[j = 1, j <= Length[horizons[[i]]], j++, (*go to j position on section*)
                    If[i == 1, time[[i]][[j]] = {(j - 1) dx, 0}, (*first horizon is surface so time = 0*)
                        vave = Table[Mean[Flatten[model[[All, j]][[k, All, 4]]]], {k, 2, i}]; (*find average velocity in layers laying above i horizon *)
                        dh = Table[Abs[horizons[[k, j, 2]] - horizons[[k - 1, j, 2]]], {k, 2, i}]; (*find thicknesess of layers laying above i horizon*)
                        If[Length[dh] == Length[vave], dt = Abs[N[2 dh/vave]], Return["mistake"]]; (*make table of dt*)
                        
                        time[[i]][[j]] = {(j - 1) dx, Total[dt]} (*for each horizon for each picket on section find time as sum of dt*)
                        ]
                    ]
                ];
                
                Return[<|"time" -> time|>]
]


BuildWellSet[horizons_ ,time_, wellCount_, wellType_] := 
Module[{
                i,
                k,
                dx,
                len,
                table,
                positions,                           
                dataset
},
                
                dx = horizons[[1, 2, 1]] - horizons[[1, 1, 1]]; (*section step*)
                len = horizons[[1, -1, 1]]; (*length of section*)
                
                If[wellType == "max",
                    positions = Last @ SortBy[Length] @ Table[Round[FindPeaks[horizons[[i]][[All, 2]], 2]][[All, 1]], {i, Length[horizons]}];(*find max values of the horizon with the highet dispersion*)
                    If[positions == {}, Return["No peaks!"], If[Length[positions] > wellCount, positions = Sort[RandomSample[positions, wellCount]]]], (*choose part of found values determined by wellCount*)
                    If[wellType == "min", 
                        positions = Last @ SortBy[Length] @ Table[Round[FindPeaks[-horizons[[i]][[All, 2]], 2]][[All, 1]], {i, Length[horizons]}]; (*find min values of the horizon with the highet dispersion*)
                        If[positions == {}, Return["No peaks!"], If[Length[positions] > wellCount, positions = Sort[RandomSample[positions, wellCount]]]], (*choose part of found values determined by wellCount*)
                        If[wellType == "regular",
                            positions = Range[Round[len/dx/wellCount/2], len/dx, Round[len/dx/wellCount]], (*wells on the regular net*)
                            If[wellType == "random",
                                positions = Sort[RandomSample[Range[1, len/dx, 1], wellCount]], (*wells positions are random*)
                                Print["wellType undefined"]]] (*or else*)
                ]
                ];
                               
                table = Flatten[Table[Table[{k, dx (positions[[k]] - 1), i, N[Interpolation[horizons[[i]], dx *(positions[[k]] - 1)]], N[Interpolation[time[[i]], dx * (positions[[k]] - 1)]]}, {i, Length[horizons]}],{k, Length[positions]}], 1]; (*table (wellNum, position on section, horizon, depth, time)*)
				dataset = Dataset[Map[<|"well" -> #[[1]], "x" -> #[[2]], "horizon" -> #[[3]], "depth" -> #[[4]], "time" -> #[[5]]|>&, table]]; (*make dataset of wells*)
				
                Return[<|"table" -> table, "dataset" -> dataset, "positions" -> dx positions|>]
]


CheckVelocity[wellDataset_]:= 
Module[{
                i,
                j,
                dt,
                dh,
                table,
                wellValuesSelected,
                wellCount,
                horizonsCount,
                dataset
                
},

                wellCount = Normal[Values[wellDataset]][[-1, 1]]; (*count wells*)
                horizonsCount = Max[Normal[Values[wellDataset]][[All, 3]]]; (*count horizons*)
                table = Table[{0, 0, 0, 0, 0}, {i, 1, wellCount}, {j, 1, horizonsCount - 1}]; (*there will be results*)

                For[i = 1, i <= wellCount, i++, (*take i well*)
                    wellValuesSelected = Select[Normal[Values[wellDataset]], #[[1]] == i& ]; (*take all data from this well (x, horizons, depths, times)*)
                    For[j = 1, j <= horizonsCount - 1, j++, (*go through layers*)
                        dh = -(wellValuesSelected [[j + 1, 4]] - wellValuesSelected [[j, 4]]);  (*in turn, calculate the layer thickness,*)  
                        dt = wellValuesSelected [[j + 1, 5]] - wellValuesSelected [[j, 5]];                    (*time*)
                        table[[i]][[j]] = {i, j, dh, dt, Abs[2 dh/dt]}                                    (*and interval speed*)
                        ] (*go to next layer*)
                    ];
                
                dataset = Dataset[Map[<|"well" -> #[[1]], "layer" -> #[[2]], "dh" -> #[[3]], "dt" -> #[[4]], "Vint" -> #[[5]]|>&, Flatten[table, 1]]]; (*make dataset*)
				
                Return[<|"dataset" -> dataset|>]
]


HTMethod[wellDataset_, time_, t_]:= 
Module[{
                
                i,
                j,
                dx,
                len,
                surface,
                positions,
                wellValues,
                lmSet,
                fits,
                lmParametres,
                depthObjective,
                depthPredicted,
                errors,
                interpolationErrors,
                result
},               
				dx = time[[1, 2, 1]]; (*section step*)
	            len = time[[1, -1, 1]]; (*length of section*)
	            
                surface = Values[Normal[wellDataset[Select[#horizon == 1&]][[All, {2, 4}]]]]; (*first horizon take as surface*)
                
                wellValues = Table[Values[Normal[wellDataset[Select[#horizon == i&]][[All, {5, 4}]]]], {i, 2, Max[wellDataset[All, 3]]}]; (*time, depth*)
                wellValues[[All, All, 1]] /= 2;   (*divide time*)     
                wellValues[[All, All, 2]] *= -1;   (*transform depths to positive values for easy calculations and presentable plots*) 
                
                lmSet = Table[LinearModelFit[wellValues[[i]], t, t], {i,  Length[wellValues]}];	(*evaluate linear models h(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
                fits = Table[Table[{(j - 1) dx, -(lmParametres[[i, 2]] * time[[i + 1, j, 2]] / 2 + lmParametres[[i, 1]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*make data suitable for ListLinePlot. (x, h = at+b)*)
								
					(*find errors between wells values and predicted by model*)
					positions = DeleteDuplicates[Normal[wellDataset[All, "x"]]]; (*x positions of wells*)
					depthObjective = -wellValues[[All, All, 2]];
					depthPredicted = Table[Flatten[Cases[fits[[i]], {#,__}]&/@positions, 1], {i, 1, Length[fits]}][[All, All, 2]];
					errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}], {i, Length[depthObjective]}];
					interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
					
					result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
					If[Length[surface] != 0, result = Join[{surface}, result]];
				
                Return[<|"wellValues" -> wellValues, "result" -> result, "lmSet" -> lmSet, "lmParametres" -> lmParametres|>]
]


dHdTMethod[setHT_, reference_, wellDataset_, time_, dt_]:= 
Module[{
                dx,
                len,
                wellValues,
                lmSet ,
                lmParametres,
                fits ,
                positions,
                depthObjective,
                depthPredicted,
                errors,
                interpolationErrors,
                result,
                plotData
},               
		
            dx = time[[1, 2, 1]]; (*x step on section*)
            len = time[[1, -1, 1]]; (*length of section*)
            positions = DeleteDuplicates[Normal[wellDataset[All, "x"]]]; (*x positions of wells*)
            
            wellValues = Table[Values[Normal[wellDataset[Select[#horizon == i&]][[All, {5, 4}]]]], {i, setHT}] ;(*time, depth*) 
            wellValues[[All, All, 1]] /= 2; (*divide time*)
            wellValues[[All, All, 2]] *= -1; (*transform depths to positive values for easy calculations and presentable plots*)    

            plotData = Table[wellValues[[i]]-wellValues[[1]], {i, 2, Length[wellValues]}];
            
            lmSet = Table[LinearModelFit[wellValues[[i]] - wellValues[[1]](*!!!*), dt, dt], {i,  2, Length[wellValues] }]; (*evaluate linear models dh(dt) = adt + b for each horizon*) 
            lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
            fits = Table[Table[{(j - 1) dx, reference[[j, 2]] - (lmParametres[[i, 2]] * (time[[setHT[[i + 1]], j, 2]]  - time[[setHT[[1]], j, 2]])/ 2 + lmParametres[[i, 1]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*evaluate fits (x, h = href + dt/a - b/a). make data suitable for ListLinePlot*)

	        (*find errors between wells values and predicted by model*)
            depthObjective = -wellValues[[2 ;; -1, All, 2]];
            depthPredicted = Table[Flatten[Cases[fits[[i]], {#,__}]&/@positions, 1], {i, 1, Length[fits]}][[All, All, 2]];
            errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}], {i, Length[depthObjective]}];
            interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)


            result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
				
				
            Return[<|"plotData" -> plotData, "result" -> result, "lmSet" -> lmSet, "lmParametres" -> lmParametres, "fits" -> fits|>]
				
         
]


VTMethod[wellDataset_, time_, t_]:= 
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
                wellValues,
                lmSet,
                fits,
                surface,
                depthObjective,                
                positions,
                errors,
                interpolationErrors,
                depthPredicted,
                result
                
},               
	            dx = time[[1, 2, 1]]; (*section step*)
	            len = time[[1, -1, 1]]; (*length of section*)
	            
	            surface = Values[Normal[wellDataset[Select[#horizon == 1&]][[All, {2, 4}]]]]; (*first horizon take as surface*)
	            
	            wellValues = Table[Values[Normal[wellDataset[Select[#horizon == i&]][[All, {2, 5, 4}]]]], {i, 2, Max[wellDataset[All, 3]]}]; (*well, time, depth*)
	            wellValues[[All, All, 2]] /= 2; (*divide time*)
	            wellValues[[All, All, 3]] *= -1; (*transform depths to positive values for easy calculations and presentable plots*) 
                
                table = Table[Table[{i, wellValues[[i, j, 2]], 0}, {j, Length[wellValues[[i]]]}], {i, Length[wellValues]}]; (*{horizon, time, 0}*)
                For[i = 1, i <= Length[wellValues], i++,
                    For[j = 1, j <= Length[wellValues[[i]]], j++,                                              
                        table[[i, j, 3]] = wellValues[[i, j, 3]] / wellValues[[i, j, 2]] (*evaluate velocities*)
                        ]
                    ];
 
		        lmSet = Table[LinearModelFit[table[[i]][[All, 2 ;; 3]], t, t], {i, 1, Length[table]}]; (*evaluate linear models v(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
                fits = Table[Table[{(j - 1) dx, -1/2*(lmParametres[[i, 2]] * time[[i + 1, j, 2]]^2/2  + lmParametres[[i, 1]] * time[[i + 1, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}];(*evaluate fits (x, h = at^2+bt). make data suitable for ListLinePlot *)
                
	
				(*find errors between wells values and predicted by model*)
				positions = DeleteDuplicates[Normal[wellDataset[All, "x"]]]; (*x positions of wells*)
				depthObjective = -wellValues[[All, All, 3]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#,__}]&/@positions, 1], {i, 1, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}], {i, 1, Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}];(*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
				If[Length[surface] != 0, result = Join[{surface}, result]];
				
                Return[<|"depthObjective"->depthObjective, "depthPredicted"->depthPredicted, "wellValues" -> table, "lmSet" -> lmSet, "lmParametres" -> lmParametres, "fits"-> fits, "result" -> result, "fits" -> fits,  "errors" -> errors[[All,2]],  "RMSError"->Map[StandardDeviation[#[[All,2]]]&,errors] |>]
]


dVdTMethod[setVT_, reference_, wellDataset_, time_, dt_]:= 
Module[{
               
                i,
                j,                
                dx,
                len,
                table,                
                lmParametres,
                wellValues,
                lmSet,
                fits,                
                depthObjective,                
                positions,
                errors,
                interpolationErrors,
                depthPredicted,
                result
	                           
},               
	            dx = time[[1, 2, 1]]; (*x step on section*)
	            len = time[[1, -1, 1]]; (*length of section*)
	            
	                     
	            wellValues = Table[Values[Normal[wellDataset[Select[#horizon == i&]][[All, {2, 5, 4}]]]], {i, setVT}]; (*well, time, depth*)
	            wellValues[[All, All, 2]] /= 2; (*divide time*)
	            wellValues[[All, All, 3]] *= -1; (*transform depths to positive values for easy calculations and presentable plots*)    

		        positions = DeleteDuplicates[Normal[wellDataset[All, "x"]]]; (*x positions of wells*)
		        
                table = Table[Table[{i, wellValues[[i +  1, j, 2]] - wellValues[[1, j, 2]], 0}, {j, Length[wellValues[[i]]]}], {i, Length[wellValues] - 1}]; (*{horizon, dt[i,j] = t[i,j] - tref[i,j] , 0}*)
                For[i = 1, i <= Length[wellValues] - 1, i++,
                    For[j = 1, j <= Length[wellValues[[i]]], j++,                                              
                        table[[i, j, 3]] = (wellValues[[i + 1, j, 3]] - wellValues[[1, j, 3]])  / (wellValues[[i + 1, j, 2]] - wellValues[[1, j, 3]] ) (*evaluate velocities dh/dt, dh[i,j] = href[i,j] - h[i,j]*)
                        ]
                    ];
 
		        lmSet = Table[LinearModelFit[table[[i]][[All, 2 ;; 3]], dt, dt], {i, Length[table]}]; (*evaluate linear models v(dt) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
                fits = Table[Table[{(j - 1) dx, reference[[j, 2]] - 1/2*(lmParametres[[i, 2]] * (time[[setVT[[i + 1]], j, 2]] - time[[setVT[[1]], j, 2]])^2/2  + lmParametres[[i, 1]] * (time[[setVT[[i + 1]], j, 2]] - time[[setVT[[1]], j, 2]]))}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*evaluate fits (x, h = href + adt^2 + bdt). make data suitable for ListLinePlot *)
                 
	
				(*find errors between wells values and predicted by model*)
				depthObjective = -wellValues[[2 ;; -1, All, 3]]; 
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#,__}]&/@positions, 1], {i, 1, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}], {i, 1, Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
				
				
                Return[<|"plotData" -> table, "lmSet" -> lmSet, "lmParametres" -> lmParametres, "result" -> result, "fits" -> fits |>]
]


THMethod[wellDataset_, time_, h_]:= 
Module[{
                i,
                j,
                dx,
                len,
                wellValues,
                lmSet,
                fits,
                lmParametres,
                surface,
                positions,
                depthObjective,
                depthPredicted,
                errors,
                interpolationErrors,
                result                
},               
				dx = time[[1, 2, 1]]; (*x step on section*)
	            len = time[[1, -1, 1]]; (*length of section*)
                
                wellValues = Table[Values[Normal[wellDataset[Select[#horizon == i&]][[All, {4, 5}]]]], {i, 2, Max[wellDataset[All, 3]]}]; (*h, t*)
                wellValues[[All, All, 2]] /= 2;  (*divide time*)
                wellValues[[All, All, 1]] *= -1;  (*transform depths to positive values for easy calculations and presentable plots*)           
                
                surface = Values[Normal[wellDataset[Select[#horizon == 1&]][[All, {2, 4}]]]]; (*{horizon, time, 0}*)
                
                lmSet = Table[LinearModelFit[wellValues[[i]], h, h], {i, Length[wellValues]}]; (*evaluate linear models v(t) for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)
                fits = Table[Table[{(j - 1) dx, -(time[[i + 1, j, 2]] / lmParametres[[i, 2]] / 2 - lmParametres[[i, 1]] / lmParametres[[i, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*evaluate fits (x, h = t/a - b/a). make data suitable for ListLinePlot*)
								
				(*find errors between wells values and predicted by model*)
				positions = DeleteDuplicates[Normal[wellDataset[All, "x"]]]; (*x positions of wells*)
				depthObjective = -wellValues[[All, All, 1]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#,__}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}],{i,  Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
				If[Length[surface] != 0, result = Join[{surface}, result]];
				
                Return[<|"depthObjective"->depthObjective, "depthPredicted"->depthPredicted, "wellValues" -> wellValues, "fits"-> fits, "result" -> result, "lmSet" -> lmSet, "lmParametres" -> lmParametres, "errors" -> errors[[All,2]],  "RMSError"->Map[StandardDeviation[#[[All,2]]]&,errors]|>]
]


dTdHMethod[setTH_, reference_, wellDataset_, time_, dh_]:= 
Module[{
                i,
                j,
                dx,
                len,
                wellValues,
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
				dx = time[[1, 2, 1]]; (*x step on section*)
	            len = time[[1, -1, 1]]; (*length of section*)
                
                wellValues = Table[Values[Normal[wellDataset[Select[#horizon == i&]][[All, {4, 5}]]]], {i, setTH}]; (*h, t*)
                wellValues[[All, All, 2]] /= 2;   (*divide time*)
                wellValues[[All, All, 1]] *= -1;   (*transform depths to positive values for easy calculations and presentable plots*)             
                
                positions = DeleteDuplicates[Normal[wellDataset[All, "x"]]]; (*x positions of wells*)
                
                plotData = Table[wellValues[[i]]-wellValues[[1]], {i, 2, Length[wellValues]}];
                
                lmSet = Table[LinearModelFit[wellValues[[i]] - wellValues[[1]], dh, dh], {i, 2, Length[wellValues]}];	(*evaluate linear models dt(dh) = adh + b for each horizon*)
                lmParametres = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}]; (*linear model fits parameters*)                
                fits = Table[Table[{(j - 1) dx,  reference[[j, 2]] - ((time[[setTH[[i + 1]], j, 2]] - time[[setTH[[1]], j, 2]])/ lmParametres[[i, 2]] / 2 - lmParametres[[i, 1]] / lmParametres[[i, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*evaluate fits (x, h = href + dt/a - b/a). make data suitable for ListLinePlot*)
								
				(*find errors between wells values and predicted by model*)
				depthObjective = -wellValues[[2 ;; -1, All, 1]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#, __}]&/@positions, 1], {i, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}],{i,  Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[lmParametres]}]; (*add errors to fits*)
								
                Return[<|"plotData" -> plotData, "result" -> result, "lmSet" -> lmSet, "lmParametres" -> lmParametres, "fits" -> fits|>]
]


VaveMethod[wellDataset_,  time_]:= 
Module[{
                
                i,
                j,
                dx,
                len,
                wellValues,
                surface,
                positions,
                depthObjective,
                errors,
                interpolationErrors,
                result,
                vAveTable,
                fits,
                depthPredicted
},               
				dx = time[[1, 2, 1]]; (*x step on section*)
	            len = time[[1, -1, 1]]; (*length of section*)
                
                wellValues = Table[Values[Normal[wellDataset[Select[#horizon == i&]][[All, {4, 5}]]]], {i, 2, Max[wellDataset[All, 3]]}]; (*pairs (depth, time)*)
                wellValues[[All, All, 2]] /= 2;   (*divide time*)
                wellValues[[All, All, 1]] *= -1;   (*transform depths to positive values for easy calculations and presentable plots*)     
                
                surface = Values[Normal[wellDataset[Select[#horizon == 1&]][[All, {2, 4}]]]]; (*(x, depth) of first horizon*)
                
                vAveTable = Table[Mean[Table[-(-wellValues[[i, j, 1]] - surface[[j, 2]])/wellValues[[i, j, 2]], {j, Length[wellValues[[i]]]}]], {i, Length[wellValues]}]; (*define average velocity for each horizon for each position on section*)
                fits = Table[Table[{(j - 1) dx, -vAveTable[[i]] * time[[i + 1, j, 2]] / 2}, {j, len/dx + 1}], {i, Length[vAveTable]}]; (*evaluate fits (x, h = vAve*t/2). make data suitable for ListLinePlot*)
								
				(*find errors between wells values and predicted by model*)
				positions = DeleteDuplicates[Normal[wellDataset[All, "x"]]]; (*x positions of wells*)
				depthObjective = -wellValues[[All, All, 1]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#,__}]&/@positions, 1], {i, 1, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}],{i, 1, Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[fits]}]; (*add errors to fits*)
				If[Length[surface] != 0, result = Join[{surface}, result]];
				
                Return[<|"depthObjective"->depthObjective, "depthPredicted"->depthPredicted, "fits"-> fits, "result" -> result, "vAveTable" -> vAveTable, "fits" -> fits , "errors" -> errors[[All,2]],  "RMSError"->Map[StandardDeviation[#[[All,2]]]&,errors]|>]
]


dVaveMethod[setVave_, reference_, wellDataset_, time_]:= 
Module[{
                
                i,
                j,
                dx,
                len,
                wellValues,
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
				dx = time[[1, 2, 1]]; (*x step on section*)
	            len = time[[1, -1, 1]]; (*length of section*)
                
                wellValues = Table[Values[Normal[wellDataset[Select[#horizon == i&]][[All, {4, 5}]]]], {i, setVave}]; (*pairs (depth, time)*)
                wellValues[[All, All, 2]] /= 2;  (*divide time*)
                wellValues[[All, All, 1]] *= -1;  (*transform depths to positive values for easy calculations and presentable plots*)        
                
                positions = DeleteDuplicates[Normal[wellDataset[All, "x"]]]; (*x positions of wells*)
                
                vAveTable = Table[Mean[Table[Abs[wellValues[[i, j, 1]] - wellValues[[1, j, 1]]]/(wellValues[[i, j, 2]] -  wellValues[[1, j, 2]]), {j, Length[wellValues[[i]]]}]], {i, 2 , Length[wellValues] }]; (*define average velocity for each horizon for each position on section. wellValues[[1]] - well values of reference horizon *)
                
                fits = Table[Table[{(j - 1) dx, reference[[j, 2]] - vAveTable[[i]] * (time[[setVave[[i + 1]], j, 2]] - time[[setVave[[1]], j, 2]])/ 2}, {j, len/dx + 1}], {i, Length[vAveTable]}]; (*evaluate fits (x, h = href + vAve*dt/2). make data suitable for ListLinePlot*)
								 
				(*find errors between wells values and predicted by model*)
				depthObjective = -wellValues[[2;;-1, All, 1]];
				depthPredicted = Table[Flatten[Cases[fits[[i]], {#,__}]&/@positions, 1], {i, 1, Length[fits]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (depthObjective - depthPredicted)[[i, j]]}, {j, 1, Length[positions]}],{i, 1, Length[depthObjective]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx + 1}], {i, Length[errors]}]; (*interpolate errors*)
				
				result = Table[Table[{fits[[i, j, 1]], (fits[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx + 1}], {i, Length[fits]}]; (*add errors to fits*)
								
                Return[<|"result" -> result, "vAveTable" -> vAveTable, "fits" -> fits|>]
]


AllMethods[set_, horizons_, wellDataset_, time_, dt_, dh_]:= 
Module[{
                reference,
                setHT,
                setVT ,
                setTH,
                setVave
                               
},               	                              
                reference = horizons[[set[[1]]]];
                setHT = Join[{set[[1]]}, Cases[set[[2 ;; -1]], {__, 1}][[All, 1]]]; (*make set of horizons (reference and predicted) for dHdT method*)
                setVT = Join[{set[[1]]}, Cases[set[[2 ;; -1]], {__, 2}][[All, 1]]]; (*make set of horizons (reference and predicted) for dVdT method*)
                setTH = Join[{set[[1]]}, Cases[set[[2 ;; -1]], {__, 3}][[All, 1]]]; (*make set of horizons (reference and predicted) for dTdH method*)
                setVave = Join[{set[[1]]}, Cases[set[[2 ;; -1]], {__, 4}][[All, 1]]]; (*make set of horizons (reference and predicted) for dVave method*)
                
                (*get results*)
                
                Return[<|"plotDataTH" -> dTdHMethod[setTH, reference, wellDataset, time, dh][["plotData"]], 
                         "plotDataVT" -> dVdTMethod[setVT, reference, wellDataset, time, dt][["plotData"]],
                         "plotDataHT" -> dHdTMethod[setHT, reference, wellDataset, time, dt][["plotData"]],
                         "resultHT" -> dHdTMethod[setHT, reference, wellDataset, time, dt][["result"]], 
                         "resultVT" -> dVdTMethod[setVT, reference, wellDataset, time, dt][["result"]],
                         "resultTH" -> dTdHMethod[setTH, reference, wellDataset, time, dh][["result"]], 
                         "lmSetHT" -> dHdTMethod[setHT, reference, wellDataset, time, dt][["lmSet"]],
                         "lmParametresHT" -> dHdTMethod[setHT, reference, wellDataset, time, dt][["lmParametres"]],
                         "lmSetVT" -> dVdTMethod[setVT, reference, wellDataset, time, dt][["lmSet"]],
                         "lmParametresVT" -> dVdTMethod[setVT, reference, wellDataset, time, dt][["lmParametres"]],
                         "lmSetTH" -> dTdHMethod[setTH, reference, wellDataset, time, dh][["lmSet"]],
                         "lmParametresTH" -> dTdHMethod[setTH, reference, wellDataset, time, dh][["lmParametres"]],
                         "resultVave" -> dVaveMethod[setVave, reference, wellDataset, time][["result"]],
                         "vAveTable" -> dVaveMethod[setVave, reference, wellDataset, time][["vAveTable"]]|>]
	 
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
