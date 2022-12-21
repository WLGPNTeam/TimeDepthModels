(* ::Package:: *)

(* ::Chapter:: *)
(*BuildSection*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[BuildDepthSection, BuildVelocitySection]


BuildDepthSection::usage = 
"BuildDepthSection[listH, alpha, len, dx, anoMaxDisp, hTapering]"

BuildVelocitySection::usage =
"BuildVelocitySection[horNH, listV]"


(* ::Section:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)

BuildDepthSection[listH_, alpha_, len_, dx_, anoMaxDisp_, hTapering_] := 
Module[{
                horNH,
                anoStartFiltRad = 5,
                surfStartFiltRad = 5,
                surfMaxDisp = 1/2 listH[[1]],
                surfDefault = 5,
                surfARParam = {0.5, 10},
                anoARParam = {0.1, 1},
                daySurface = Table[5, len/dx],
                rfDaySurface,
                rfAno,
                localAnomaly = Table[0, len/dx], 
                localAnomalyFiltered,
                horNHsorted,
                i = 1,
                max,
                listSumH,
                table
 },
                rfDaySurface = RandomFunction[ARProcess[{anoARParam[[1]]}, anoARParam[[2]]], {1, len}];
                daySurface = surfMaxDisp*GaussianFilter[rfDaySurface["SliceData", Range[1, len, dx]], surfStartFiltRad][[1]];
                table =Join[{daySurface},Table[-Sum[listH[[i]], {i, 1, i}], {i, Length[listH]}, {j, len/dx}]];
                horNH = Table[(table[[i]][[j]] - (len - j*dx)*Tan[alpha]), {i, Length[listH] + 1}, {j, len/dx}];
                rfAno = RandomFunction[ARProcess[{surfARParam[[1]]}, surfARParam[[2]]], {1, len}];
                localAnomaly = anoMaxDisp*GaussianFilter[rfAno["SliceData", Range[1, len, dx]], anoStartFiltRad][[1]];
                listSumH = Join[{1}, Table[Sum[listH[[i]], {i, 1, i}], {i, 1, Length[listH]}]];
    
                ];
                            i++;
    While[i < Length[horNH] + 1, 
                            max = Max[Table[localAnomaly[[j]] - localAnomalyFiltered[[j]], {j, 1, Length[localAnomaly]}]];
    
                horNHsorted = Table[{dx*j, horNH[[i]][[j]]}, {i, Length[listH] + 1}, {j, len/dx}];
    
                Return[<|"horNH" -> horNH, "horNHsorted" -> horNHsorted|>]
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


BuildVelocitySection[horNH_, listV_] := 
Module[{
                velModel,
                i, j, dh,
                f,
                sortedH,
                maxH,
                horNHforVelMod,
                startV,
                layerthickness
},
                f[dh_] := 50*Sqrt[dh];(*в параметры перенести*)
                sortedH = Table[{dx*j, horNH[[i]][[j]]}, {i, Length[listH] + 1}, {j, Length[horNH[[1]]]}];
                maxH = -Min[sortedH[[-1]]];
                velModel = Table[{0, 0, 0}, {i, 1, Length[horNH], 1}, {j, 1, Length[horNH[[1]]]}, {dh, 1, maxH}];
                horNHforVelMod = Join[horNH, {Table[-(maxH + 10), {i, 1, Length[horNH[[1]]]}]}];
                For[i = 1, i <= Length[horNHforVelMod] - 1, i++,
                    startV = listV[[i]];
                    For[j = 1, j <= Length[horNH[[1]]], j++,
                        layerthickness = Ceiling[Abs[horNHforVelMod[[i + 1]][[j]] - horNHforVelMod[[i]][[j]]]];
                        For[dh = 0, dh <= layerthickness, dh++, 
                            velModel[[i]][[j]][[dh + 1]] = {j*dx, Round[horNH[[i, j]]] - dh, Round[Evaluate[startV + f[dh]]]};
                        ]
                    ]
                ];  
                
                velModel = Table[DeleteCases[velModel[[i, j]], {0, 0, 0}], {i, 1, Length[horNH]}, {j, 1, Length[horNH[[1]]]}];
                
                <|"velModel" -> velModel|>
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
