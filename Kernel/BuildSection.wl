(* ::Package:: *)

(* ::Chapter:: *)

(*BuildSection*)

(* ::Section:: *)

(*Package*)

BeginPackage["WLGPNTeam`TimeDepthModels`"]

(* ::Section:: *)

(*Names*)

ClearAll[BuildDepthSection]

BuildDepthSection::usage = "BuildDepthSection[numOfLayers, listH, alpha, len, dx, anoMaxDisp, hTapering]"

(* ::Section:: *)

(*Private context*)

Begin["`Private`"]

(* ::Section:: *)

(*Implementation*)

BuildDepthSection[numOfLayers_, listH_, alpha_, len_, dx_, anoMaxDisp_,
     hTapering_] :=
    Module[{horizonNH, anoStartFiltRad = 5, surfStartFiltRad = 5, surfMaxDisp
         = 1/2 listH[[1]], surfDefault = 5, surfARParam = {0.5, 10}, anoARParam
         = {0.1, 1}, daySurface = Table[5, len / dx], rfDaySurface, rfAno, localAnomaly
         = Table[0, len / dx], localAnomalyFiltered, horNHsorted, i = 1, max,
         listSumH, table},
        rfDaySurface = RandomFunction[ARProcess[{anoARParam[[1]]}, anoARParam
            [[2]]], {1, len}];
        daySurface = surfMaxDisp * GaussianFilter[rfDaySurface["SliceData",
             Range[1, len, dx]], surfStartFiltRad][[1]];
        table = Join[{daySurface}, Table[-Sum[listH[[i]], {i, 1, i}],
             {i, numOfLayers}, {j, len / dx}]];
        horizonNH = Table[(table[[i]][[j]] - (len - j * dx) * Tan[alpha
            ]), {i, numOfLayers + 1}, {j, len / dx}];
        rfAno = RandomFunction[ARProcess[{surfARParam[[1]]}, surfARParam
            [[2]]], {1, len}];
        localAnomaly = anoMaxDisp * GaussianFilter[rfAno["SliceData",
             Range[1, len, dx]], anoStartFiltRad][[1]];
        listSumH = Join[{1}, Table[Sum[listH[[i]], {i, 1, i}], {i, 1,
             numOfLayers}]];
        While[
            i < Length[horizonNH] + 1
            ,
            horizonNH[[-i]] += localAnomaly;
            localAnomalyFiltered = anoMaxDisp * GaussianFilter[rfAno[
                "SliceData", Range[1, len, dx]], (anoStartFiltRad + hTapering TotalH 
                / listSumH[[-i]])][[1]];
            max = Max[Table[localAnomaly[[j]] - localAnomalyFiltered[[
                j]], {j, 1, Length[localAnomaly]}]];
            localAnomaly = Table[localAnomalyFiltered[[j]] + max, {j,
                 1, Length[localAnomalyFiltered]}];
            i++;
        ];
        horNHsorted = Table[{dx * j, horizonNH[[i]][[j]]}, {i, numOfLayers
             + 1}, {j, len / dx}];
        <|"horizonNH" -> horizonNH, "horNHsorted" -> horNHsorted|>
    ]

(* ::Section:: *)

(*End private*)

End[]

(* ::Section:: *)

(*End package*)

EndPackage[]