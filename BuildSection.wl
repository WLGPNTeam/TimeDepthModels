BuildDepthSection[numOfLayers_, listH_, alpha_, len_, dx_, anoMaxDisp_, hTapering_] := Module[ (*количество слоев, мощности, угол наклона, длина профиля, 
  шаг, степень сглаживания*)
  {
   horizonNH,
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
  
  table = Join[{daySurface}, Table[-Sum[listH[[i]], {i, 1, i}], {i, numOfLayers}, {j, len/dx}]];
  horizonNH = Table[(table[[i]][[j]] - (len - j*dx)*Tan[alpha]), {i, numOfLayers + 1}, {j, len/dx}];
    
  rfAno = RandomFunction[ARProcess[{surfARParam[[1]]}, surfARParam[[2]]], {1, len}];
  localAnomaly = anoMaxDisp*GaussianFilter[rfAno["SliceData", Range[1, len, dx]], anoStartFiltRad][[1]];
  listSumH = Join[{1}, Table[Sum[listH[[i]], {i, 1, i}], {i, 1, numOfLayers}]];
  While[i < Length[horizonNH] + 1,
   
   horizonNH[[-i]] += localAnomaly;
   localAnomalyFiltered = anoMaxDisp*GaussianFilter[rfAno["SliceData", Range[1, len, dx]], (anoStartFiltRad + hTapering TotalH/listSumH[[-i]])][[1]];
   max = Max[Table[localAnomaly[[j]] - localAnomalyFiltered[[j]], {j, 1, Length[localAnomaly]}]];
   localAnomaly = Table[localAnomalyFiltered[[j]] + max, {j, 1, Length[localAnomalyFiltered]}];
   i++;
   ];
  
  horNHsorted = Table[{dx*j, horizonNH[[i]][[j]]}, {i, numOfLayers + 1}, {j, len/dx}];
  
  <|"horizonNH" -> horizonNH, "horNHsorted" -> horNHsorted|>
  	]

BuildVelocitySection[horNH_, len_, dx_, listV_] := Module[ (*убрать len и dx  из параметров и ненужные переменные внутри*)
  {
   velModel,
   i, j,
   sortedH,
   maxH,
   dopH,
   horNHforVelMod,
   startV0,
   layerthickness,
   modelForPlot
   },
  f[dh_] := 50*Sqrt[dh];(* в параметры перенести*)
  sortedH = Table[{dx*j, horNH[[i]][[j]]}, {i, numOfLayers + 1}, {j, len/dx}];
  maxH = -Min[sortedH[[-1]]];
  velModel = Table[{0, 0, 0}, {i, 1, Length[horNH], 1}, {j, 1, len/dx, 1}, {dh, 1, maxH}];
  dopH = Table[-(maxH + 10), {i, 1, len/dx, 1}];
  horNHforVelMod = Join[horNH, {dopH}];
  For[i = 1, i <= Length[horNHforVelMod] - 1, i++,
   startV0 = V0[[i]];
   
   For[j = 1, j <= len/dx, j++, 
    layerthickness = Ceiling[Abs[horNHforVelMod[[i + 1]][[j]] - horNHforVelMod[[i]][[j]]]]; 
    dh = Ceiling[layerthickness/10];
    For[dh = 0, dh <= layerthickness, dh++, 
     velModel[[i]][[j]][[dh + 1]] = {j*dx, Round[horNH[[i, j]]] - dh, Round[Evaluate[startV0 + f[dh]]]};
     ]
    ]
   ];
  
  modelForPlot = Table[DeleteCases[velModel[[i, j]], {0, 0, 0}], {i, 1, 
     Length[horNH], 1}, {j, 1, len/dx, 1}];
  <|"velModel" -> velModel, "modelForPlot" -> modelForPlot|>
  ]