(* ::Package:: *)

(* ::Chapter:: *)
(*BuildSection*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[BuildDepthSection, BuildVelocitySection, BuildTimeSection, BuildWellSet, CheckVelocity, HTMethod, VTMethod, THMethod, VaveMethod]


BuildDepthSection::usage = 
"BuildDepthSection[listH, alpha, len, dx, anoMaxDisp, hTapering, type]"


BuildVelocitySection::usage =
"BuildVelocitySection[horNH, listV]"


BuildTimeSection::usage = 
"BuildTimeSection[horNH, velModel]"


BuildWellSet::usage = 
"BuildWellSet[horNH, timeNH, numOfWells, wellsLocation]"


CheckVelocity::usage = 
"CheckVelocity[datasetWells]"


HTMethod::usage = 
"HTMethod[datasetWells, timeNH, t]"


VTMethod::usage = 
"VTMethod[datasetWells, datasetVelModel, timeNH, t]"


THMethod::usage = 
"THMethod[datasetWells, timeNH, t]"


VaveMethod::usage = 
"VaveMethod[datasetWells, timeNH]"


(* ::Section::Closed:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


BuildDepthSection[listH_, alpha_, len_, dx_, anoMaxDisp_, hTapering_, type_] := 
Module[{
                i,
                horNH,
                horizons,
                anoStartFiltRad,
                anoARParam,
                daySurface,    
                rfAno,
                localAnomaly, 
                localAnomalyFiltered,
                listSumH,
                horizontal,
                horizonts,
                TotalH,
                max,
                numOfhor
},
				If[hTapering >= 1/2 ((len/dx)- anoStartFiltRad) listH[[1]]/TotalH , Return["hTapering is too big"]]; 
 
				anoStartFiltRad = 5; (*GaussianFilter radius for layer with the highest dispersion*)
				anoARParam = {0.1, 1}; (*parametres for ARProcess*)
                daySurface = Table[0, len/dx + 1];  (*day surface as it is (without taking in account relief); len/dx = num of pickets*)  
                localAnomaly = Table[0, len/dx + 1]; (*init localAnomaly*)
                TotalH = Total[listH]; 
                
                (*horizontal layers made*)       
                horizontal = Join[{daySurface}, Table[-Sum[listH[[i]], {i, 1, i}], {i, Length[listH]}, {j, len/dx + 1}]]; 
                
                (*add dipping if alpha != 0*)
                horizonts = Table[(horizontal[[i]][[j]] - (len - (j - 1) dx) Tan[alpha]), {i, Length[listH] + 1}, {j, len/dx + 1}]; 
                
                (*make the highest anomaly*)
                rfAno = RandomFunction[ARProcess[{anoARParam[[1]]}, anoARParam[[2]]], {1, len + 1}]; 
                localAnomaly = anoMaxDisp GaussianFilter[rfAno["SliceData", Range[1, len + 1, dx]], anoStartFiltRad][[1]];
                
                listSumH = Join[{1}, Table[Sum[listH[[i]], {i, 1, i}], {i, 1, Length[listH]}]]; (*needed for increasing Gaussian Filter radius in While below*)
                
                
                i = 1;
                While[i < Length[horizonts] + 1, 
                            horizonts[[-type i]] += localAnomaly; (*"type" defines which (last or first) hor has the highest anomaly; add anomaly*)
                            
                               localAnomalyFiltered = anoMaxDisp GaussianFilter[rfAno["SliceData",Range[1, len + 1, dx]],(anoStartFiltRad + hTapering TotalH/listSumH[[-i]])][[1]];(*changing local anomaly for each horizon*)
                               max = Abs[Max[Table[type*(localAnomaly[[j]] - localAnomalyFiltered[[j]]), {j, 1, Length[localAnomaly]}]]]; (*finding max difference between previous and current anomaly; need to know so that the horizons do not intersect*)
                               localAnomaly = Table[localAnomalyFiltered[[j]] + type * max, {j, 1, Length[localAnomalyFiltered]}]; (*evaluate anomaly for the next horizon*)
                               i++;
                ];
    
                horNH = Table[{(j - 1) dx, Round[horizonts[[i]][[j]]]}, {i, Length[listH] + 1}, {j, len/dx + 1}]; (*make a list suitable for the Plot*)
    
                Return[<|"horNH" -> horNH|>]
]


BuildVelocitySection[horNH_, listV_] := 
Module[{
                i, 
                j, 
                dh,
                dx,
                f,
                maxH,
                vAve,
                horNHJoined,
                layerthickness,
                assocVelModel,
                datasetVelModel,
                velModel
                
},
                dx = horNH[[1, 2, 1]] - horNH[[1, 1, 1]]; (*section step*)
                maxH = -Min[horNH[[-1]]]; (*find max depth*)
                                                              
                f[vAve_, dh_] := vAve + 100 Sqrt[dh] + RandomInteger[{-vAve/100, vAve/100}];(*vAve - velocity from listV for each layer, trend - velocity depth trend, anom - anomaly*)
                                              
                velModel = Table[{0, 0, 0, 0}, {i, 1, Length[horNH], 1}, {j, 1, Length[horNH[[i]]]}, {dh, 1, maxH}]; 
                (*init velModel [num of horizons x num of pickets x layer thickness]; use maxH here beccause thickness is different at each (i,j) position; zeros will be deleted after)*)
                horNHJoined = Join[Part[horNH, All, All, -1], {Table[-(maxH + 1), {i, 1, Length[horNH[[1]]]}]}]; (*modificated horNH*)
               
                For[i = 1, i <= Length[horNH], i++, (*go into layer*)
                    vAve = listV[[i]]; (*init vAve*)
                    For[j = 1, j <= Length[horNH[[i]]], j++, (*go to j position on section*)
                        layerthickness = Ceiling[Abs[horNHJoined[[i + 1]][[j]] - horNHJoined[[i]][[j]]]]; (*evaluate i layer thickness on j position*)
                        For[dh = 0, dh <= layerthickness, dh++, (*go vertical and evaluate velocity on each dh step*)
                            velModel[[i]][[j]][[dh + 1]] = {i, (j - 1) dx, Round[horNHJoined[[i, j]]] - dh, Round[f[vAve, dh]]}; (*layer, x coordinate, h cootdinate, velocity*)
                        ]
                    ]
                ];  
                
                velModel = Table[DeleteCases[velModel[[i, j]], {0, 0, 0, 0}], {i, 1, Length[horNH]}, {j, 1, Length[horNH[[i]]]}]; (*delete (0,0,0,0) cases*)
                assocVelModel = Map[<|"layer" -> #[[1]], "x" -> #[[2]], "h" -> #[[3]], "v" -> #[[4]]|>&, Flatten[velModel, 2]];(*make association*)
				datasetVelModel = Dataset[assocVelModel]; (*make dataset*)
                
                <|"velModel" -> velModel, "datasetVelModel" -> datasetVelModel|>
]


BuildTimeSection[horNH_, velModel_] := 
Module[{
                timeNH,
                dx,
                i,
                j,
                dv,
                dh,
				dt
},
                dx = horNH[[1, 2, 1]] - horNH[[1, 1, 1]]; (*section step*)
                timeNH = Table[{0, 0}, {i, 1, Length[horNH]}, {j, 1, Length[horNH[[1]]]}]; (*init timeNH [num of horizons x num of pickets] *)
        
       
                For[i = 1, i <= Length[horNH], i++, 
                For[j = 1, j <= Length[horNH[[i]]], j++,
                    If[i == 1, timeNH[[i]][[j]] = {(j - 1) dx, 0}, (*first horizon is relief so time = 0*)
                        dv = Table[Mean[Flatten[Part[velModel[[All]][[All, j]]][[All]][[k, All, 4]]]], {k, i}]; (*find average velocities in layers laying above i horizon *)
                        dh = Table[Abs[horNH[[k, j, 2]]],{k, i}]; (*find thicknesess of layers laying above i horizon*)
                        If[Length[dh] == Length[dv], dt = Abs[N[2 dh/dv]], Return["mistake"]]; (*make table of dt*)
                        timeNH[[i]][[j]] = {(j - 1) dx, Total[dt]} (*for each horizon for each picket on section find time as sum of dt*)
                        ]
                    ]
                ];
                
                Return[<|"timeNH" -> timeNH|>]
]


BuildWellSet[horNH_ ,timeNH_, numOfWells_, wellsLocation_] := 
Module[{
                wells,
                positions,
                dx,
                len,
                i,
                k,                
                datasetWells
},
                
                dx = horNH[[1, 2, 1]] - horNH[[1, 1, 1]]; (*section step*)
                len = horNH[[1, -1, 1]]; (*length of section*)
                
                If[wellsLocation == "max",
                    positions = Last @ SortBy[Length] @ Table[Round[FindPeaks[horNH[[i]][[All, 2]], 2]][[All, 1]], {i, Length[horNH]}];(*find max values of the horizon with the highet dispersion*)
                    If[positions == {}, Return["No peaks!"], If[Length[positions] > numOfWells, positions = Sort[RandomSample[positions, numOfWells]]]], (*choose part of found values determined by num of wells*)
                    If[wellsLocation == "min", 
                        positions = Last @ SortBy[Length] @ Table[Round[FindPeaks[-horNH[[i]][[All, 2]], 2]][[All, 1]], {i, Length[horNH]}]; (*find min values of the horizon with the highet dispersion*)
                        If[positions == {}, Return["No peaks!"], If[Length[positions] > numOfWells, positions = Sort[RandomSample[positions, numOfWells]]]], (*choose part of found values determined by num of wells*)
                        If[wellsLocation == "regular",
                            positions = Range[Round[len/dx/numOfWells/2], len/dx, Round[len/dx/numOfWells]], (*wells on the regular net*)
                            If[wellsLocation == "random",
                                positions = Sort[RandomSample[Range[1, len/dx, 1], numOfWells]], (*wells positions are random*)
                                Print["wellsLocation undefined"]]]
                ]
                ];
                               
                wells = Flatten[Table[Table[{k, dx (positions[[k]] - 1), i, N[Interpolation[horNH[[i]], dx *(positions[[k]] - 1)]], N[Interpolation[timeNH[[i]], dx *(positions[[k]]-1)]]}, {i, Length[horNH]}],{k, Length[positions]}], 1]; (*table (wellNum, position on section, horizon, depth, time)*)
				datasetWells = Dataset[Map[<|"well" -> #[[1]], "x" -> #[[2]], "horizon" -> #[[3]], "depth" -> #[[4]], "time" -> #[[5]]|>&, wells]]; (*make dataset of wells*)
				
                Return[<|"wells" -> wells, "datasetWells" -> datasetWells, "positions" -> dx positions|>]
]


CheckVelocity[datasetWells_]:= 
Module[{
                dt,
                dh,
                listVint,
                valuesWells,
                valuesWellsPart,
                keysWells,
                i,
                j,
                numOfWells,
                firstWell,
                firstHor,
                numOfHor,
                assocVint,
                datasetVint
                
},

                valuesWells = Normal[Values[datasetWells]];
                
                numOfWells = valuesWells[[-1, 1]];
                firstWell = valuesWells[[1, 1]];
                firstHor = valuesWells[[1, 3]];
                numOfHor = Max[valuesWells[[All, 3]]];
                listVint = Table[{0, 0, 0, 0, 0}, {i, 1, numOfWells}, {j, 1, numOfHor - 1}];

                For[i = firstWell, i <= numOfWells, i++,
                    valuesWellsPart = Select[valuesWells, #[[1]] == i& ];
                    For[j = firstHor, j <= numOfHor - 1, j++,
                        dh = -(valuesWellsPart [[j + 1, 4]] - valuesWellsPart [[j, 4]]); 
                        dt = valuesWellsPart [[j + 1, 5]] - valuesWellsPart [[j, 5]];
                        listVint[[i]][[j]] = {i, j, dh, dt, Abs[2 dh/dt]}
                        ]
                    ];
                
                listVint = Flatten[listVint, 1];
				assocVint = Map[<|"well" -> #[[1]], "layer" -> #[[2]], "dh" -> #[[3]], "dt" -> #[[4]], "Vint" -> #[[5]]|>&, listVint];
				datasetVint = Dataset[assocVint];
				
                Return[<|"datasetVint" -> datasetVint|>]
]


HTMethod[datasetWells_, timeNH_, t_]:= 
Module[{
                dataWellsHT,
                lmSet,
                i,
                j,
                fitsHT,
                fitsHTParams,
                dx,
                len,
                relief,
                positions,
                wellsSet,
                htSet,
                errors,
                interpolationErrors,
                resultHT
},               
				dx = timeNH[[1, 2, 1]];
	            len = timeNH[[1, -1, 1]];
                
                dataWellsHT = Table[Values[Normal[datasetWells[Select[#horizon == i&]][[All, {5, 4}]]]], {i, 2, Max[datasetWells[All, 3]]}];               
                relief = Values[Normal[datasetWells[Select[#horizon == 1&]][[All, {2, 4}]]]];
                lmSet = Table[LinearModelFit[dataWellsHT[[i]], t, t], {i, 1, Length[dataWellsHT]}];	
                fitsHTParams = Table[lmSet[[i]]["BestFitParameters"], {i, 1, Length[lmSet]}];
                fitsHT = Table[Table[{(j - 1) dx, (fitsHTParams[[i, 2]] * 2 timeNH[[i + 1, j, 2]] + fitsHTParams[[i, 1]])}, {j, len/dx}], {i, Length[fitsHTParams]}];
								
				(*find errors*)
				positions = DeleteDuplicates[Normal[datasetWells[All, "x"]]];
				wellsSet = dataWellsHT[[All, All, 2]];
				htSet = Table[Flatten[Cases[fitsHT[[i]], {#,__}]&/@positions, 1], {i, 1, Length[fitsHT]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (wellsSet - htSet)[[i, j]]},{j, 1, Length[positions]}],{i, 1, Length[wellsSet]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx}], {i, Length[errors]}];
				
				resultHT = Join[{relief}, Table[Table[{fitsHT[[i, j, 1]], (fitsHT[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx}], {i, Length[fitsHTParams]}]];
				
                Return[<|"resultHT" -> resultHT, "lmSet" -> lmSet, "fitsHTParams" -> fitsHTParams|>]
]


VTMethod[datasetWells_, datasetVelModel_, timeNH_, t_]:= 
Module[{
               
                i,
                j,
                x,
                h,
                v,
                dx,
                len,
                setVT,
                datasetVT,
                fitParams,
                wellsXTH,
                lmSet,
                fitsVT,
                relief,
                wellsSet,                
                positions,
                errors,
                interpolationErrors,
                vtSet,
                resultVT
                
},               
	            dx = timeNH[[1, 2, 1]];
	            len = timeNH[[1, -1, 1]];
	            relief = Values[Normal[datasetWells[Select[#horizon == 1&]][[All, {2, 4}]]]];
	            
	            wellsXTH = Table[Values[Normal[datasetWells[Select[#horizon == i&]][[All, {2, 5, 4}]]]], {i, 2, Max[datasetWells[All, 3]]}];
                setVT = Table[Table[{i, wellsXTH[[i, j, 2]], 0},{j, Length[wellsXTH[[i]]]}], {i, Length[wellsXTH]}];
                For[i = 1, i <= Length[wellsXTH], i++,
                    For[j = 1, j <= Length[wellsXTH[[i]]], j++,
                        x = wellsXTH[[i, j, 1]];
                        h = wellsXTH[[i, j, 3]];
                        v = Flatten[Values[Normal[datasetVelModel[Select[(#layer == i + 1 && #x == Round[x] && #h == Round[h]) &]]]], 2];
                        setVT[[i, j, 3]]= v[[4]]
                        ]
                    ];
 
		        lmSet = Table[LinearModelFit[setVT[[i]][[All, 2;;3]], t, t], {i, 1, Length[setVT]}];
                fitParams = Table[lmSet[[i]]["BestFitParameters"], {i, Length[lmSet]}];
                fitsVT = Table[Table[{(j - 1) dx, -(fitParams[[i, 2]] * (timeNH[[i + 1, j, 2]])^2/8 + fitParams[[i, 1]] * timeNH[[i + 1, j, 2]]/4)}, {j, len/dx}], {i, Length[fitParams]}];
                
	
				(*find errors*)
				positions = DeleteDuplicates[Normal[datasetWells[All, "x"]]];
				wellsSet = wellsXTH[[All, All, 3]];
				vtSet = Table[Flatten[Cases[fitsVT[[i]], {#,__}]&/@positions, 1], {i, 1, Length[fitsVT]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (wellsSet - vtSet)[[i, j]]},{j, 1, Length[positions]}],{i, 1, Length[wellsSet]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx}], {i, Length[errors]}];
				
				resultVT = Join[{relief}, Table[Table[{fitsVT[[i, j, 1]], (fitsVT[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx}], {i, Length[fitParams]}]];
				
                Return[<| "lmSet" -> lmSet, "fitParams" -> fitParams, "resultVT" -> resultVT |>]
]


THMethod[datasetWells_, timeNH_, h_]:= 
Module[{
                dataWellsTH,
                lmSet,
                i,
                j,
                fitsTH,
                fitsTHParams,
                dx,
                len,
                relief,
                positions,
                wellsSet,
                thSet,
                errors,
                interpolationErrors,
                resultTH                
},               
				dx = timeNH[[1, 2, 1]]; (*x step on section*)
	            len = timeNH[[1, -1, 1]]; (*length of section*)
                
                dataWellsTH = Table[Values[Normal[datasetWells[Select[#horizon == i&]][[All, {4, 5}]]]], {i, 2, Max[datasetWells[All, 3]]}]; (*h, t*)              
                relief = Values[Normal[datasetWells[Select[#horizon == 1&]][[All, {2, 4}]]]];
                lmSet = Table[LinearModelFit[dataWellsTH[[i]], h, h], {i, 1, Length[dataWellsTH]}];	
                fitsTHParams = Table[lmSet[[i]]["BestFitParameters"], {i, 1, Length[lmSet]}];
                fitsTH = Table[Table[{(j - 1) dx, (2timeNH[[i + 1, j, 2]]/fitsTHParams[[i, 2]] - fitsTHParams[[i, 1]]/fitsTHParams[[i, 2]])}, {j, len/dx}], {i, Length[fitsTHParams]}];
								
				(*find errors*)
				positions = DeleteDuplicates[Normal[datasetWells[All, "x"]]];
				wellsSet = dataWellsTH[[All, All, 1]];
				thSet = Table[Flatten[Cases[fitsTH[[i]], {#,__}]&/@positions, 1], {i, 1, Length[fitsTH]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (wellsSet - thSet)[[i, j]]},{j, 1, Length[positions]}],{i, 1, Length[wellsSet]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx}], {i, Length[errors]}];
				
				resultTH = Join[{relief}, Table[Table[{fitsTH[[i, j, 1]], (fitsTH[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx}], {i, Length[fitsTHParams]}]];
				
                Return[<|"resultTH" -> resultTH, "lmSet" -> lmSet, "fitsTHParams" -> fitsTHParams|>]
]


VaveMethod[datasetWells_, timeNH_]:= 
Module[{
                
                i,
                j,
                dx,
                len,
                tableTH,
                relief,
                positions,
                hFromWells,
                htSet,
                errors,
                interpolationErrors,
                result,
                vAve,
                hPredictTable,
                hPredict
},               
				dx = timeNH[[1, 2, 1]];
	            len = timeNH[[1, -1, 1]];
                
                tableTH = Table[Values[Normal[datasetWells[Select[#horizon == i&]][[All, {4, 5}]]]], {i, 2, Max[datasetWells[All, 3]]}]; (*pairs (time, depth)*)          
                relief = Values[Normal[datasetWells[Select[#horizon == 1&]][[All, {2, 4}]]]]; (*(x, depth) of first horizon*)
                
                vAve = Table[Mean[Table[-2tableTH[[i, j, 1]]/tableTH[[i, j, 2]], {j, Length[tableTH[[i]]]}]], {i, Length[tableTH]}]; (*define average velocity for each horizon for each position on section; depth should be withminus always otherwise...*)
                hPredictTable = Table[Table[{(j - 1) dx, -vAve[[i]] * timeNH[[i + 1, j, 2]]}, {j, len/dx}], {i, Length[vAve]}];
								
				(*find errors*)
				positions = DeleteDuplicates[Normal[datasetWells[All, "x"]]];
				hFromWells = tableTH[[All, All, 1]];
				hPredict = Table[Flatten[Cases[hPredictTable[[i]], {#,__}]&/@positions, 1], {i, 1, Length[hPredictTable]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (hFromWells - hPredict)[[i, j]]},{j, 1, Length[positions]}],{i, 1, Length[hFromWells]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx}], {i, Length[errors]}];
				
				result = Join[{relief}, Table[Table[{hPredictTable[[i, j, 1]], (hPredictTable[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx}], {i, Length[hPredictTable]}]];
				
                Return[<|"result" -> result, "vAve" -> vAve, "hPredictTable" -> hPredictTable|>]
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
