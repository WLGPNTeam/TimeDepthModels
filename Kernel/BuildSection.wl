(* ::Package:: *)

(* ::Chapter:: *)
(*BuildSection*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[BuildDepthSection, BuildVelocitySection, BuildTimeSection, BuildWellSet]


BuildDepthSection::usage = 
"BuildDepthSection[listH, alpha, len, dx, anoMaxDisp, hTapering]"


BuildVelocitySection::usage =
"BuildVelocitySection[horNHsorted, listV]"


BuildTimeSection::usage = 
"BuildTimeSection[horNHsorted, velModel]"


BuildWellSet::usage = 
"BuildWellSet[horNHsorted, timeNH, numOfWells, wellsLocation]"


(* ::Section::Closed:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


BuildDepthSection[listH_, alpha_, len_, dx_, anoMaxDisp_, hTapering_] := 
Module[{
                horNH,
                anoStartFiltRad,
                anoARParam,
                daySurface,    
                rfAno,
                localAnomaly, 
                localAnomalyFiltered,
                horNHsorted,
                listSumH,
                table,
                TotalH,
                max,
                i
 },
				If[hTapering >= 1/2 ((len/dx)-anoStartFiltRad) listH[[1]]/TotalH , Return["hTapering is too big"]];
 
				anoStartFiltRad = 5;
				anoARParam = {0.1, 1};
                daySurface = Table[0, len/dx + 1];    
                localAnomaly = Table[0, len/dx + 1];
                TotalH = Total[listH];
                 
                table =Join[{daySurface}, Table[-Sum[listH[[i]], {i, 1, i}], {i, Length[listH]}, {j, len/dx + 1}]];
                horNH = Table[(table[[i]][[j]] - (len - (j-1) dx) Tan[alpha]), {i, Length[listH] + 1}, {j, len/dx + 1}];
                rfAno = RandomFunction[ARProcess[{anoARParam[[1]]}, anoARParam[[2]]], {1, len+1}];
                localAnomaly = anoMaxDisp GaussianFilter[rfAno["SliceData", Range[1, len+1, dx]], anoStartFiltRad][[1]];
                listSumH = Join[{1}, Table[Sum[listH[[i]], {i, 1, i}], {i, 1, Length[listH]}]];
                
                i = 1;
                While[i < Length[horNH]+1, 
                            horNH[[-i]] += localAnomaly;
                            localAnomalyFiltered = anoMaxDisp GaussianFilter[rfAno["SliceData",Range[1, len+1, dx]],(anoStartFiltRad + hTapering TotalH/listSumH[[-i]])][[1]];
                            max = Max[Table[localAnomaly[[j]] - localAnomalyFiltered[[j]], {j, 1, Length[localAnomaly]}]];
                            localAnomaly = Table[localAnomalyFiltered[[j]] + max, {j, 1, Length[localAnomalyFiltered]}];
                            i++;
                ];
    
                horNHsorted = Table[{(j-1) dx, Round[horNH[[i]][[j]]]}, {i, Length[listH] + 1}, {j, len/dx+1}];
    
                Return[<|"horNH" -> horNH, "horNHsorted" -> horNHsorted|>]
]

BuildVelocitySection[horNHsorted_, listV_] := 
Module[{
                velModel,
                i, 
                j, 
                dh,
                f,
                maxH,
                horNHforVelMod,
                startV,
                layerthickness,
                horNH
                
},
                dx = horNHsorted[[1,2,1]] - horNHsorted[[1,1,1]];
                f[dh_] := 30 Sqrt[dh];(*velocity trend in layer, should be into Module params*)
                horNH = Part[horNHsorted, All, All,-1];
                maxH = -Min[horNHsorted[[-1]]];
                velModel = Table[{0, 0, 0}, {i, 1, Length[horNH], 1}, {j, 1, Length[horNH[[1]]]}, {dh, 1, maxH}];
                horNHforVelMod = Join[horNH, {Table[-(maxH + 1), {i, 1, Length[horNH[[1]]]}]}];
               
                For[i = 1, i <= Length[horNHforVelMod] - 1, i++,
                    startV = listV[[i]];
                    For[j = 1, j <= Length[horNH[[1]]], j++,
                        layerthickness = Ceiling[Abs[horNHforVelMod[[i + 1]][[j]] - horNHforVelMod[[i]][[j]]]];
                        For[dh = 0, dh <= layerthickness, dh++, 
                            velModel[[i]][[j]][[dh + 1]] = {j dx, Round[horNH[[i, j]]] - dh, Round[startV + f[dh]]};
                        ]
                    ]
                ];  
                
                velModel = Table[DeleteCases[velModel[[i, j]], {0, 0, 0}], {i, 1, Length[horNH]}, {j, 1, Length[horNH[[1]]]}];
                
                <|"velModel" -> velModel|>
]

BuildTimeSection[horNHsorted_, velModel_] := 
Module[{
                timeNH,
                dx,
                i,
                j,
                dv,
                dh,
				dt
},
                dx = horNHsorted[[1, 2, 1]] - horNHsorted[[1, 1, 1]];
                timeNH = Table[{0, 0}, {i, 1, Length[horNHsorted]}, {j, 1, Length[horNHsorted[[1]]]}];
        
                For[i = 1, i <= Length[horNHsorted], i++, 
                For[j = 1, j <= Length[horNHsorted[[i]]], j++,
                    If[i == 1, timeNH[[i]][[j]] = {(j - 1) dx, 0},
                        dv = Table[Mean[Flatten[Part[velModel[[All]][[All, j]]][[All]][[k, All, 3]]]],{k, i}];
                        dh = Table[Abs[horNHsorted[[k, j, 2]]],{k,i}];
                        If[Length[dh]==Length[dv],dt=N[-2dh/dv],Return["mistake"]];
                        timeNH[[i]][[j]] = {(j - 1) dx, Total[dt]}
                        ]
                    ]
                ];
                
                Return[<|"timeNH" -> timeNH|>]
]

BuildWellSet[horNHsorted_ ,timeNH_, numOfWells_, wellsLocation_] := 
Module[{
                wells,
                positions,
                dx,
                len,
                i,
                k,
                lm,
                wellsNums,
                assocWells,
                datasetWells
},
                
                dx = horNHsorted[[1, 2, 1]] - horNHsorted[[1, 1, 1]];
                len = horNHsorted[[1, -1, 1]];
                
                If[wellsLocation == "max",
                    positions = Round[FindPeaks[horNHsorted[[-1]][[All, 2]], 2]];
                    If[positions == {}, Return["No peaks!"], If[Length[positions] > numOfwells, positions = Sort[RandomSample[positions, numOfwells]]]],
                    If[wellsLocation == "min", 
                        positions = Round[FindPeaks[-horNHsorted[[-1]][[All, 2]], 2]];
             If[positions == {}, Return["No peaks!"], If[Length[positions] > numOfwells, positions = Sort[RandomSample[positions, numOfwells]]]],
             If[wellsLocation == "regular",
                            positions = Table[{Range[2, len/dx, Round[len/dx/numOfWells]][[i]], 0},{i, numOfWells}],
                            If[wellsLocation == "random",
                                positions = Table[{RandomSample[Range[2, len/dx, 3], numOfWells][[i]], 0},{i, numOfWells}],
                                Print["wellsLocation undefined"]
                            ]
                        ]
                    ]
                ];
               
                Print[positions];
                wellsNums = Table[Table[{Range[1, Length[positions]][[j]], dx (positions[[j, 1]]-1)},{i,Length[horNHsorted]}],{j, Length[positions]}];
                Table[lm[i_, x_]:={Interpolation[horNHsorted[[i]], x], Interpolation[timeNH[[i]], x]},{i, Length[horNHsorted]}];
                wells = Flatten[Table[Table[{wellsNums[[k, i, 1]], wellsNums[[k, i, 2]], StringJoin["Hor ", ToString[i]], N[lm[i, k][[1]]], N[lm[i, k][[2]]]},{i, Length[horNHsorted]}],{k, Length[positions]}],1];
				assocWells = Map[<|"well" -> #[[1]], "x" -> #[[2]], "hor" -> #[[3]], "depth" -> #[[4]], "time" -> #[[5]]|>&, wells];
				datasetWells = Dataset[assocWells];
				
                Return[<|"wells" -> wells, "datasetWells" -> datasetWells|>]
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
