  PlotSection[horNHsorted_] := Module[
<<<<<<< f8072b740f5b74be51de6c321943471169735a53
  {},
=======
  {
      
    },
>>>>>>> update
  ListLinePlot[horNHsorted,
                GridLinesStyle -> Directive[Thick, Blue],
   	            FrameStyle -> Directive[Black, 18],
               	Filling -> Bottom, Frame -> True, ImageSize -> 800,
   	            PlotLabels -> Map["Hor " <> ToString[#] &, (Range[numOfLayers + 1] - 1)],
               	PlotLabel -> "Depth Model", 
                LabelStyle -> Directive[18, Bold, Gray]]
  ]

PlotVelocity[velModel_] := Module[
  {

     },
  ListContourPlot[Flatten[velModel, 2], PlotTheme -> "Detailed"]
  ]