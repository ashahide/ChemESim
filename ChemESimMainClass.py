from ProcessModelingAndControlClass import TransferFunction

Num = 3

Den = [20, 3, 3]

Sys1 = TransferFunction(Numerator = Num, Denominator = Den)

Sys2 = TransferFunction(Numerator = Num, Denominator = Den, TimeDelay=5, N=1)

Sys3 = TransferFunction(Numerator = Num, Denominator = Den)

Sys4 = TransferFunction(Numerator = Num, Denominator = Den, TimeDelay=5, N=1)

Sys5 = TransferFunction(Numerator = Num,  Denominator = Den)

Sys6 = TransferFunction(Numerator = Num, Denominator = Den, TimeDelay=5, N=1)

import numpy as np
NumPoints = 200
TFinal    = 50        

T = np.linspace(0, TFinal, NumPoints)

U = np.zeros_like(T)

U[0:50] = 1

Output1, Time1 = Sys1.InputFunction(Num, NumPoints, TFinal, 'Step')

Output2, Time = Sys2.InputFunction(Num, NumPoints, TFinal, 'Step')

Output3, Time = Sys3.InputFunction(Num, NumPoints, TFinal, 'Impulse')

Output4, Time = Sys4.InputFunction(Num, NumPoints, TFinal, 'Impulse')

Output5, Time = Sys5.InputFunction(Num, NumPoints, TFinal, 'Square', T= Time1, U = U)

Output6, Time = Sys6.InputFunction(Num, NumPoints, TFinal, 'Square', T= Time1, U = U)

Sys1.PlotResponse()

from ProcessModelingAndControlClass import CompareResults

CompareResults(Sys1, Sys2, Sys3, Sys4, Sys5, Sys6, YUnit = 'Output')





    
