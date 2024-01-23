(* ::Package:: *)

(* ::Chapter:: *)
(*BuildSection*)


(* ::Section:: *)
(*Package*)


BeginPackage["WLGPNTeam`TimeDepthModels`"]


(* ::Section:: *)
(*Names*)


ClearAll[DepthSection, VelocitySection, TimeSection, WellSet]


DepthSection::usage = 
"DepthSection[listH, slopes, len, dx, type, hDispersion, hTapering, gaussRadius]"


(*DepthSection has another variant for building combined section*)


VelocitySection::usage =
"VelocitySection[horizons, listV, trends, anomalies]"


TimeSection::usage = 
"TimeSection[horizons, model]"


WellSet::usage = 
"WellSet[horizons, timeSection, count, type]"


(* ::Section::Closed:: *)
(*Private context*)


Begin["`Private`"]


(* ::Section:: *)
(*Implementation*)


DepthSection[listH_, slopes_, len_, dx_, type_, hDispersion_, foldingParameters_] := 
Module[{
                i,
                j,
                horizons,
                hTapering,
                radius,
                ARparameters,
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
		       
		       If[hDispersion != 0,
		        hTapering = foldingParameters[[1]];
				radius = foldingParameters[[2]]; 
				If[hTapering >= 1/2 ((len/dx) - radius) listH[[1]]/totalH , Return["hTapering is too big"]]
				]; 
 
				(*GaussianFilter radius for layer with the highest dispersion*)
				ARparameters = {0.1, 1}; (*parametres for ARProcess*) (*hmmmm....*)
                
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
						(*else*)
						For[j = 1, j <= len/dx + 1, j++, table[[i, j]] = -N[Tan[slopes[[i]]]*(j - 1) * dx]];
						dh1 = Abs[Max[Table[table[[i, j]] - table[[i - 1, j]], {j, len/dx + 1}]]];
						dh2 = sums[[i - 1]];(*If[dh1 >= dh2, max = dh1, max = dh2];*)
					    table[[i, All]] = Function[x, x - Max[dh1, dh2] * 1.1]/@table[[i, All]]
                    ] 
                  ]
                ];
                (*make the highest anomaly*)
                If[hDispersion != 0,
                  anomaly = Table[0, len/dx + 1]; (*init anomaly*);
                  data = RandomFunction[ARProcess[{ARparameters[[1]]}, ARparameters[[2]]], {1, len + 1}]; 
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
                If[MatchQ[table[[1]], 
                Table[0., {j, len/dx + 1}]], 
                horizons = Table[{(j - 1) dx, Round[table[[i, j]]]}, {i, Length[listH] + 1}, {j, len/dx + 1}],
                max = Max[table[[1]]]; (*need to shift horizons under {0,0...0}*)
                surface = Table[{(j - 1) dx, 0}, {j, len/dx + 1}];
                horizons = Join[{surface}, Table[{(j - 1) dx, (Round[table[[i, j]]] - 1.1 * max)}, {i, Length[listH] + 1}, {j, len/dx + 1}]]                  
                  ];
                                    
                Return[<|"horizons" -> horizons, "table" -> table|>]
]


VelocitySection[horizons_, listV_, trendsQ_, anomaliesQ_] := 
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
                
                f[vAve_, trendQ_, anomalyQ_, dh_] := vAve + trendQ * 10*Sqrt[dh] +  RandomReal[{-anomalyQ*vAve, anomalyQ*vAve}];(*vAve - velocity from listV for each layer, trendsQ - depth trend exist or not, anomaliesQ - anomalies exist or not*)
                                              
                model = Table[{0, 0, 0, 0}, {i, Length[horizons]}, {j, Length[horizons[[i]]]}, {dh, maxH}]; 
                (*init model [horizons count x pickets count x layer thickness]; use maxH here beccause thickness is different at each (i,j) position; zeros will be deleted after)*)
                horizonsJoined = Join[Part[horizons, All, All, -1], {Table[-(maxH + 1), {j, Length[horizons[[1]]]}]}]; (*modificated horizons*)
               
                For[i = 1, i <= Length[horizons], i++, (*go into layer*)
                    vAve = listV[[i]];(*init vAve for the layer*)
                    trendQ = trendsQ[[i]];(*trendQ = 1 or 0*)
                    anomalyQ = anomaliesQ[[i]]; (*anomalyQ = 1 or 0*)
                    For[j = 1, j <= Length[horizons[[i]]], j++, (*go to j position on section*)
                        layerThickness = Ceiling[Abs[horizonsJoined[[i + 1]][[j]] - horizonsJoined[[i]][[j]]]]; (*define i layer thickness on j position*)
                        For[dh = 0, dh <= layerThickness, dh++, (*go vertical and evaluate velocity on each dh step*)
                            model[[i]][[j]][[dh + 1]] = {i, (j - 1) dx, Round[horizonsJoined[[i, j]]] - dh, N[f[vAve, trendQ, anomalyQ, dh]]}; (* = {layer, x coordinate, h cootdinate, velocity}*)
                        ]
                    ]
                ];  
                
                model = Table[DeleteCases[model[[i, j]], {0, 0, 0, 0}], {i, Length[horizons]}, {j, Length[horizons[[i]]]}]; (*delete (0, 0, 0, 0) cases*)
                dataset = Dataset[Map[<|"layer" -> #[[1]], "x" -> #[[2]], "z" -> -#[[3]], "v" -> #[[4]]|>&, Flatten[model, 2]]]; (*make dataset*)
                
                <|"model" -> model, "dataset" -> dataset|>
]


TimeSection[horizons_, model_] := 
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
                time = Table[{0, 0}, {i, 1, Length[horizons]}, {j, 1, Length[horizons[[1]]]}]; (*init time [num of horizons x num of pickets] *)
        
       
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


WellSet[horizons_, timeSection_, count_, type_] := 
Module[{
                i,
                k,
                dx,
                len,
                table,
                positions,                           
                ds
},
                
                dx = horizons[[1, 2, 1]] - horizons[[1, 1, 1]]; (*section step*)
                len = horizons[[1, -1, 1]]; (*length of section*)
                
                If[type == "max",
                    positions = Last @ SortBy[Length] @ Table[Round[FindPeaks[horizons[[i]][[All, 2]], 2]][[All, 1]], {i, Length[horizons]}];(*find max values of the horizon with the highet dispersion*)
                    If[positions == {}, Return["No peaks!"], If[Length[positions] > count, positions = Sort[RandomSample[positions, count]]]], (*choose part of found values determined by wellCount*)
                    If[type == "min", 
                        positions = Last @ SortBy[Length] @ Table[Round[FindPeaks[-horizons[[i]][[All, 2]], 2]][[All, 1]], {i, Length[horizons]}]; (*find min values of the horizon with the highet dispersion*)
                        If[positions == {}, Return["No peaks!"], If[Length[positions] > count, positions = Sort[RandomSample[positions, count]]]], (*choose part of found values determined by wellCount*)
                        If[type == "regular",
                            positions = Range[Round[len/dx/count/2], len/dx, Round[len/dx/count]], (*wells on the regular net*)
                            If[type == "random",
                                positions = Sort[RandomSample[Range[1, len/dx, 1], count]], (*wells positions are random*)
                                Print["wellType undefined"]]] (*or else*)
                ]
                ];
                               
                table = Flatten[Table[Table[{k, dx (positions[[k]] - 1), i, N[Interpolation[horizons[[i]], dx *(positions[[k]] - 1)]], N[Interpolation[timeSection[[i]], dx * (positions[[k]] - 1)]]}, {i, Length[horizons]}], {k, Length[positions]}], 1]; (*table (wellNum, position on section, horizon, depth, time)*)
				ds = Dataset[Map[<|"well" -> #[[1]], "x" -> #[[2]], "horizon" -> #[[3]], "depth" -> #[[4]], "time" -> #[[5]]|>&, table]]; (*make dataset of wells*)
				
                Return[<|"table" -> table, "ds" -> ds, "positions" -> dx positions|>]
]


DepthSection[horizonsTop_, horizonsBot_]:=
Module[{
                trend,
                horizonsT,
                horizonsResult,
                length,
                dx,
                t,
                i,
                j,
                tab,
                surface
},
                length = Length[horizonsBot[[1]]];
                horizonsT = horizonsTop;
                horizonsResult = horizonsBot;
                (*added a trend from the last horizon of the first set to the 2nd horizon of the second set*)
                trend = # - Mean[horizonsT[[-1, All, 2]]]&/@horizonsT[[-1, All, 2]]; 
                
                horizonsResult[[2, All, 2]] += trend;   (*1st is zero*)

                (*check if the 2nd and 3rd horizons of the second set intersect, add the delta*)
                t = Table[horizonsResult[[3, j, 2]] - horizonsResult[[2, j, 2]], {j, length}]; 
                If[Not[AllTrue[t, # < 0&]], 
                    For[j = 1, j <= Length[horizonsResult[[2]]], j++, horizonsResult[[2, j, 2]] += 1.1 Max[t]]
                    ];

                (*remove the uppermost horizon*)
                horizonsResult = Delete[horizonsResult, 1];
                horizonsT = Delete[horizonsT, 1];

                (*creating a table of dh values for the remaining horizons above*)
                tab = Table[Table[{horizonsT[[-i - 1, j, 1]], -horizonsT[[-i, j, 2]] + horizonsT[[-i - 1, j, 2]]}, {j, length}], {i, Length[horizonsT] - 1}];

                (*forming the overlying horizons
                check that each new horizon does not intersect with the underlying one
                move it up if it intersects by the value of the intersection x 1.1*)
                For[i = 1, i <= Length[tab], i++, z = Table[{tab[[i, j, 1]], horizonsResult[[1, j, 2]] + tab[[i, j, 2]]}, {j, length}];
                    horizonsResult = Join[{z}, horizonsResult];
                    t = Table[horizonsResult[[2, j, 2]] - horizonsResult[[1, j, 2]], {j, length}];
                    If[Not[AllTrue[t, # < 0&]], For[j = 1, j <= Length[horizonsResult[[1]]], j++, horizonsResult[[1, j, 2]] += 1.1 Max[t]]]
                ];

                (*all horizons get down zero*)
                t = horizonsResult[[1, All, 2]];
                If[Not[AllTrue[t, # < 0&]], 
                    For[i = 1,i <= Length[horizonsResult],i++,
                        For[j = 1, j <= Length[horizonsResult[[i]]], j++, horizonsResult[[i, j, 2]] -= 1.1 Max[t]
                        ]
                    ]
                ];

                (*add surface (zeros)*)
                dx = horizonsBot[[1, 2, 1]]- horizonsBot[[1, 1, 1]];
                surface = Table[{(j - 1) dx, 0}, {j, length + 1}];
                horizonsResult = Join [{surface}, horizonsResult]; 
                
                Return[<|"horizons" -> horizonsResult|>]
]


(* ::Section:: *)
(*End private*)


End[]


(* ::Section:: *)
(*End package*)


EndPackage[]
