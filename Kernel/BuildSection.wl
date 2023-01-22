(* ::Package:: *)

(* ::Chapter:: *)
(*BuildSection*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[BuildDepthSection, BuildVelocitySection, BuildTimeSection, BuildWellSet, CheckVelocity, HT, VT]


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


HT::usage = 
"HT[datasetWells, timeNH, t]"


VT::usage = 
"VT[datasetWells, datasetVelModel, timeNH, t]"


(* ::Section::Closed:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


BuildDepthSection[listH_, alpha_, len_, dx_, anoMaxDisp_, hTapering_, type_] := 
Module[{
                horNH,
                horizons,
                anoStartFiltRad,
                anoARParam,
                daySurface,    
                rfAno,
                localAnomaly, 
                localAnomalyFiltered,
                listSumH,
                table,
                TotalH,
                max,
                i,
                numOfhor
},
				If[hTapering >= 1/2 ((len/dx)- anoStartFiltRad) listH[[1]]/TotalH , Return["hTapering is too big"]];
 
				anoStartFiltRad = 5;
				anoARParam = {0.1, 1};
                daySurface = Table[0, len/dx + 1];    
                localAnomaly = Table[0, len/dx + 1];
                TotalH = Total[listH];
                listSumH = Table[Sum[listH[[i]], {i, 1, i}], {i, 1, Length[listH]}];
                
                table =Join[{daySurface}, Table[-Sum[listH[[i]], {i, 1, i}], {i, Length[listH]}, {j, len/dx + 1}]];
                horizons = Table[(table[[i]][[j]] - (len - (j-1) dx) Tan[alpha]), {i, Length[listH] + 1}, {j, len/dx + 1}];
                rfAno = RandomFunction[ARProcess[{anoARParam[[1]]}, anoARParam[[2]]], {1, len+1}];
                localAnomaly = anoMaxDisp GaussianFilter[rfAno["SliceData", Range[1, len+1, dx]], anoStartFiltRad][[1]];
                listSumH = Join[{1}, Table[Sum[listH[[i]], {i, 1, i}], {i, 1, Length[listH]}]];
                
                i = 1;
                While[i < Length[horizons] + 1, 
                            horizons[[-type*i]] += localAnomaly;
                            
                               localAnomalyFiltered = anoMaxDisp GaussianFilter[rfAno["SliceData",Range[1, len + 1, dx]],(anoStartFiltRad + hTapering TotalH/listSumH[[-i]])][[1]];
                               max = Abs[Max[Table[type*(localAnomaly[[j]] - localAnomalyFiltered[[j]]), {j, 1, Length[localAnomaly]}]]];
                               localAnomaly = Table[localAnomalyFiltered[[j]] + type*max, {j, 1, Length[localAnomalyFiltered]}];
                               i++;
                ];
    
                horNH = Table[{(j-1) dx, Round[horizons[[i]][[j]]]}, {i, Length[listH] + 1}, {j, len/dx + 1}];
    
                Return[<|"horNH" -> horNH|>]
]


BuildVelocitySection[horNH_, listV_] := 
Module[{
                velModel,
                i, 
                j, 
                dh,
                dx,
                f,
                maxH,
                horNHforVelMod,
                startV,
                layerthickness,
                horNHpart,
                assocVelModel,
                datasetVelModel
                
},
                dx = horNH[[1, 2, 1]] - horNH[[1, 1, 1]];
                f[dh_] := startV + 100 Sqrt[dh] + RandomInteger[{-startV/100, startV/100}];(*velocity trend in layer, should be into Module params*)
                horNHpart = Part[horNH, All, All, -1];
                maxH = -Min[horNH[[-1]]];
                velModel = Table[{0, 0, 0, 0}, {i, 1, Length[horNHpart], 1}, {j, 1, Length[horNHpart[[1]]]}, {dh, 1, maxH}];
                horNHforVelMod = Join[horNHpart, {Table[-(maxH + 1), {i, 1, Length[horNHpart[[1]]]}]}];
               
                For[i = 1, i <= Length[horNHforVelMod] - 1, i++,
                    startV = listV[[i]];
                    For[j = 1, j <= Length[horNHpart[[1]]], j++,
                        layerthickness = Ceiling[Abs[horNHforVelMod[[i + 1]][[j]] - horNHforVelMod[[i]][[j]]]];
                        For[dh = 0, dh <= layerthickness, dh++, 
                            velModel[[i]][[j]][[dh + 1]] = {i, (j-1) dx, Round[horNHpart[[i, j]]] - dh, Round[f[dh]]};
                        ]
                    ]
                ];  
                
                velModel = Table[DeleteCases[velModel[[i, j]], {0, 0, 0, 0}], {i, 1, Length[horNHpart]}, {j, 1, Length[horNHpart[[1]]]}];
                assocVelModel = Map[<|"layer" -> #[[1]], "x" -> #[[2]], "h" -> #[[3]], "v" -> #[[4]]|>&, Flatten[velModel, 2]];
				datasetVelModel = Dataset[assocVelModel];
                
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
                dx = horNH[[1, 2, 1]] - horNH[[1, 1, 1]];
                timeNH = Table[{0, 0}, {i, 1, Length[horNH]}, {j, 1, Length[horNH[[1]]]}];
        
                For[i = 1, i <= Length[horNH], i++, 
                For[j = 1, j <= Length[horNH[[i]]], j++,
                    If[i == 1, timeNH[[i]][[j]] = {(j - 1) dx, 0},
                        dv = Table[Mean[Flatten[Part[velModel[[All]][[All, j]]][[All]][[k, All, 4]]]], {k, i}];
                        dh = Table[Abs[horNH[[k, j, 2]]],{k, i}];
                        If[Length[dh] == Length[dv], dt = Abs[N[2 dh/dv]], Return["mistake"]];
                        timeNH[[i]][[j]] = {(j - 1) dx, Total[dt]}
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
                lmHandT,
                wellsNums,
                assocWells,
                datasetWells
},
                
                dx = horNH[[1, 2, 1]] - horNH[[1, 1, 1]];
                len = horNH[[1, -1, 1]];
                
                If[wellsLocation == "max",
                    positions = Round[FindPeaks[horNH[[-1]][[All, 2]], 2]][[All, 1]];
                    If[positions == {}, Return["No peaks!"], If[Length[positions] > numOfWells, positions = Sort[RandomSample[positions, numOfWells]]]],
                    If[wellsLocation == "min", 
                        positions = Round[FindPeaks[-horNH[[-1]][[All, 2]], 2]][[All, 1]];
                        If[positions == {}, Return["No peaks!"], If[Length[positions] > numOfWells, positions = Sort[RandomSample[positions, numOfWells]]]],
                        If[wellsLocation == "regular",
                            positions = Range[Round[len/dx/numOfWells/2], len/dx, Round[len/dx/numOfWells]],
                            If[wellsLocation == "random",
                                positions = Sort[RandomSample[Range[1, len/dx, 1], numOfWells]], (*mb povtory !!!*)
                                Print["wellsLocation undefined"]]]
                ]
                ];
               
                wellsNums = Table[Table[{Range[1, Length[positions]][[j]], dx *(positions[[j]] - 1)}, {i, Length[horNH]}], {j, Length[positions]}];
                wells = Flatten[Table[Table[{wellsNums[[k, i, 1]], wellsNums[[k, i, 2]], i, N[Interpolation[horNH[[i]], dx *(positions[[k]] - 1)]], N[Interpolation[timeNH[[i]], dx *(positions[[k]]-1)]]}, {i, Length[horNH]}],{k, Length[positions]}], 1];
				assocWells = Map[<|"well" -> #[[1]], "x" -> #[[2]], "horizon" -> #[[3]], "depth" -> #[[4]], "time" -> #[[5]]|>&, wells];
				datasetWells = Dataset[assocWells];
				
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


HT[datasetWells_, timeNH_, t_]:= 
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
                fitsHT = Table[Table[{(j - 1) dx, (fitsHTParams[[i, 2]] * timeNH[[i + 1, j, 2]] + fitsHTParams[[i, 1]])}, {j, len/dx}], {i, Length[fitsHTParams]}];
								
				(*find errors*)
				positions = DeleteDuplicates[Normal[datasetWells[All, "x"]]];
				wellsSet = dataWellsHT[[All, All, 2]];
				htSet = Table[Flatten[Cases[fitsHT[[i]], {#,__}]&/@positions,1], {i, 1, Length[fitsHT]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (wellsSet - htSet)[[i, j]]},{j, 1, Length[positions]}],{i, 1, Length[wellsSet]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx}], {i, Length[errors]}];
				
				resultHT = Join[{relief}, Table[Table[{fitsHT[[i, j, 1]], (fitsHT[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx}], {i, Length[fitsHTParams]}]];
				
                Return[<|"dataWellsHT"->dataWellsHT, "resultHT" -> resultHT, "lmSet" -> lmSet, "fitsHTParams" -> fitsHTParams|>]
]


VT[datasetWells_, datasetVelModel_, timeNH_, t_]:= 
Module[{
                wellsXTH,
                lmSet,
                i,
                j,
                setVT,
                x,
                h,
                v,
                datasetVT,
                fitParams,
                dx,
                len,
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
				vtSet = Table[Flatten[Cases[fitsVT[[i]], {#,__}]&/@positions,1], {i, 1, Length[fitsVT]}][[All, All, 2]];
				errors = Table[Table[{positions[[j]], (wellsSet - vtSet)[[i, j]]},{j, 1, Length[positions]}],{i, 1, Length[wellsSet]}];
				interpolationErrors = Table[Table[{(j - 1)dx, Interpolation[errors[[i]], (j - 1)dx, Method -> "Spline"]}, {j, len/dx}], {i, Length[errors]}];
				
				resultVT = Join[{relief}, Table[Table[{fitsVT[[i, j, 1]], (fitsVT[[i, j, 2]] + interpolationErrors[[i, j, 2]])}, {j, len/dx}], {i, Length[fitParams]}]];
				
                Return[<| "lmSet" -> lmSet, "setVT" -> setVT , "fitParams" -> fitParams, "resultVT" -> resultVT |>]
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
