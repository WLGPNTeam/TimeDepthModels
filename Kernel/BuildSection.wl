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
                anoStartFiltRad = 5,
                anoARParam = {0.1, 1},
                daySurface = Table[0, len/dx + 1],    
                rfAno,
                localAnomaly = Table[0, len/dx + 1], 
                localAnomalyFiltered,
                horNHsorted,
                i = 1,
                max,
                listSumH,
                table,
                TotalH = Total[listH]
 },
                table =Join[{daySurface}, Table[-Sum[listH[[i]], {i, 1, i}], {i, Length[listH]}, {j, len/dx + 1}]];
                horNH = Table[(table[[i]][[j]] - (len - (j-1) dx) Tan[alpha]), {i, Length[listH] + 1}, {j, len/dx + 1}];
                rfAno = RandomFunction[ARProcess[{anoARParam[[1]]}, anoARParam[[2]]], {1, len+1}];
                localAnomaly = anoMaxDisp GaussianFilter[rfAno["SliceData", Range[1, len+1, dx]], anoStartFiltRad][[1]];
                listSumH = Join[{1}, Table[Sum[listH[[i]], {i, 1, i}], {i, 1, Length[listH]}]];
    
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
                dx = horNHsorted[[1,2,1]] - horNHsorted[[1,1,1]],
                horNH
                
},
                f[dh_] := 30 Sqrt[dh];(*velocity trend in layer, should be into Module params*)
                horNH = Part[horNHsorted,All,All,-1];
                maxH = -Min[horNHsorted[[-1]]];
                velModel = Table[{0, 0, 0}, {i, 1, Length[horNH], 1}, {j, 1, Length[horNH[[1]]]}, {dh, 1, maxH}];
                horNHforVelMod = Join[horNH, {Table[-(maxH + 10), {i, 1, Length[horNH[[1]]]}]}];
               
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
                v,
                h
},
                dx = horNHsorted[[1, 2, 1]] - horNHsorted[[1, 1, 1]];
                timeNH = Table[{0, 0}, {i, 1, Length[horNHsorted]}, {j, 1, Length[horNHsorted[[1]]]}];
        
                For[i = 1, i <= Length[horNHsorted], i++, 
                For[j = 1, j <= Length[horNHsorted[[i]]], j++,
                    If[i == 1, timeNH[[i]][[j]] = {(j - 1) dx, 0},
                        v = Mean[Flatten[Part[velModel[[All]][[All, j]]][[All]][[1 ;; i, All, 3]]]];
                        h = Abs[horNHsorted[[i, j, 2]]];
                        timeNH[[i]][[j]] = {(j - 1) dx, -N[2 h/v]}
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
                lm
},
                dx = horNHsorted[[1, 2, 1]] - horNHsorted[[1, 1, 1]];
                len = horNHsorted[[1, -1, 1]];
                
                If[wellsLocation == "max",
                    positions = FindPeaks[horNHsorted[[-1]][[All, 2]], 2],
                    If[wellsLocation == "min", 
                        positions = FindPeaks[-horNHsorted[[-1]][[All, 2]],2],
                        If[wellsLocation == "regular",
                            positions = Table[{Range[2, len/dx, Round[len/dx/numOfWells]][[i]], 0},{i, numOfWells}],
                            If[wellsLocation == "random",
                                positions = Table[{RandomSample[Range[2, len/dx, 3], numOfWells][[i]], 0},{i, numOfWells}],
                                Print["wellsLocation undefined"]
                            ]
                        ]
                    ]
                ];
                
                Table[lm[i_, x_]:={Interpolation[horNHsorted[[i]], x], Interpolation[timeNH[[i]], x]},{i, Length[horNHsorted]}];
                wells = Table[{k, dx (positions[[k, 1]]-1), Table[{StringJoin["Hor ", ToString[i]], N[lm[i, k][[1]]], N[lm[i, k][[2]]]},{i, Length[horNHsorted]}]},{k, Length[positions]}];

                Return[<|"wells" -> wells|>]
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
