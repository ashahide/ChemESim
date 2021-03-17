from Process import TransferFunction

Num = 3

Den = [20, 3, 3]

Sys1 = TransferFunction(Numerator = Num, Denominator = Den)

Sys2 = TransferFunction(Numerator = Num, Denominator = Den, TimeDelay=5, N=1)

Sys3 = TransferFunction(Numerator = Num, Denominator = Den)

Sys4 = TransferFunction(Numerator = Num, Denominator = Den, TimeDelay=5, N=1)

Sys5 = TransferFunction(Numerator = Num,  Denominator = Den)

Sys6 = TransferFunction(Numerator = Num, Denominator = Den, TimeDelay=5, N=1)



NumPoints = 200
TFinal    = 50        

Output1, Time1 = Sys1.InputFunction(Num, NumPoints, TFinal, 'Step')

Output2, Time = Sys2.InputFunction(Num, NumPoints, TFinal, 'Step')

Output3, Time = Sys3.InputFunction(Num, NumPoints, TFinal, 'Impulse')

Output4, Time = Sys4.InputFunction(Num, NumPoints, TFinal, 'Impulse')

import numpy as np

U = np.zeros_like(Time)
U[0:10] = 1

Output5, Time = Sys5.InputFunction(Num, NumPoints, TFinal, 'Square', T= Time1, U = U)

Output6, Time = Sys6.InputFunction(Num, NumPoints, TFinal, 'Square', T= Time1, U = U)

from Process import CompareResults

CompareResults(Sys1, Sys2, Sys3 = Sys3, Sys4 = Sys4, Sys5 = Sys5, Sys6 = Sys6, YUnit = 'Concentration')




# Sys1.PlotResults(Output2 = Output2, Input2Type = 'Step',     Input2Magnitude = 3,   Sys2TimeDelay=5, \
#                 Output3 = Output3,  Input3Type = 'Impulse',  Input3Magnitude = 3,   Sys3TimeDelay=0, \
#                 Output4 = Output4,  Input4Type = 'Impulse',  Input4Magnitude = 3,   Input4TimeDefined = U,   Sys4TimeDelay=5, \
#                 Output5 = Output5,  Input5Type = 'Square',   Input5Magnitude = 3,   Input5TimeDefined = U, Sys5TimeDelay=0, \
#                 Output6 = Output6,  Input6Type = 'Square',   Input6Magnitude = 3,   Input6TimeDefined = U, Sys6TimeDelay=5, \
#                 YUnit = 'Concentration')

