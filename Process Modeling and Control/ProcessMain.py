from Process import ProcessControlandModeling

Sys1 = ProcessControlandModeling(Numerator = 1, Denominator = [2, 3, 3])

Sys2 = ProcessControlandModeling(Numerator = 1, Denominator = [2, 3, 3], TimeDelay=5, N=1)

Sys3 = ProcessControlandModeling(Numerator = [1, 1], Denominator = [2, 3, 3])

Sys4 = ProcessControlandModeling(Numerator = [1, 1], Denominator = [2, 3, 3])

Sys5 = ProcessControlandModeling(Numerator = 1, Denominator = [2, 3, 3], TimeDelay=5, N=1)


NumPoints = 150
TFinal    = 15         

Output1, Time = Sys1.InputFunction(3, NumPoints, TFinal, 'Step')

Output2, Time = Sys2.InputFunction(3, NumPoints, TFinal, 'Step')

Output3, Time = Sys3.InputFunction(3, NumPoints, TFinal, 'Step')

Output5, Time = Sys5.InputFunction(3, NumPoints, TFinal, 'Step')


import numpy as np

U = np.zeros_like(Time)
U[0:10] = 1

Output4, Time = Sys4.InputFunction(3, NumPoints, TFinal, 'Square', T= Time, U = U)

# Sys1.PlotResults(Output2 = Output2, Input2Type = 'Step',  Input2Magnitude = 3, Sys2TimeDelay=0, \
#                 Output3 = Output3, Input3Type = 'Step',   Input3Magnitude = 3, Sys3TimeDelay=5,  \
#                 Output4 = Output4, Input4Type = 'Square', Input4Magnitude = 3, Sys4TimeDelay=0, \
#                 YUnit = 'Concentration')


from Process import FitFOPDT


FitFOPDT(Output5, 3)
