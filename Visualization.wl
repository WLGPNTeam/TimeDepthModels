  PlotSection[horNHsorted_] := Module[
  {},
  ListLinePlot[horNHsorted,
                GridLinesStyle -> Directive[Thick, Blue],
   	            FrameStyle -> Directive[Black, 18],
               	Filling -> Bottom, Frame -> True, ImageSize -> 800,
   	            PlotLabels -> Map["Hor " <> ToString[#] &, (Range[numOfLayers + 1] - 1)],
               	PlotLabel -> "Depth Model", 
                LabelStyle -> Directive[18, Bold, Gray]]
  ]

PlotVelocity[modelForPlot_] := Module[
  {},
  ListContourPlot[Flatten[modelForPlot, 2], PlotTheme -> "Detailed"]
  ]