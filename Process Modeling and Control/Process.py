class TransferFunction:

    """
    Links to cite:

    [1] http://techteach.no/python_control/python_control.pdf

    """

    def __init__(self, Numerator, Denominator, TimeDelay = None, N = None):

        self.Numerator   = Numerator
        self.Denominator = Denominator
        self.TimeDelay   = TimeDelay
        self.N           = N #for pade approximation [1]

        """
        Using method for time delay creation from [1]

        """
        import control 

        if self.TimeDelay is not None and self.N is not None:

            PadeNumerator, PadeDenominator = control.pade(self.TimeDelay, self.N)

            G_Pade = control.tf(PadeNumerator, PadeDenominator)

            G_NoDelay = control.tf(self.Numerator, self.Denominator)

            G = G_Pade*G_NoDelay
            

        else:
            G = control.tf(self.Numerator, self.Denominator)

        self.G = G

    def InputFunction(self, Magnitude, NumPoints, TFinal, Type, T = None, U = None):

        import control.matlab
        import control
        import numpy as np 
        

        t = np.linspace(0, TFinal, NumPoints)

        if Type == 'Step':
            Output, Time = control.matlab.step(self.G, t)

        elif Type == 'Impulse':
            Output, Time = control.matlab.impulse(self.G, t)

        elif Type == 'Square':
            T, Output, Time  = control.forced_response(self.G, T=T, U=U)

        self.Magnitude = Magnitude
        self.InputFunction      = Type
        self.Output            = Magnitude*Output
        self.Time              = Time
        self.T  = T
        self.U  = U


        return Magnitude*Output, Time

    def PolesAndZeros(self):
        import control.matlab

        Poles = control.matlab.pole(self.G)
        Zeros = control.matlab.zero(self.G)

        print("Poles: {}".format(Poles))
        print("Zeros: {}".format(Zeros))
        
        return Poles, Zeros

    def RootLocus(self):

        import control.matlab
        import matplotlib.pyplot as plt 

        control.matlab.rlocus(self.G)

        plt.show()

        return

    def BodePlot(self, Units = None):
        import control

        if Units == 'dB':
            BodeMagnitude, BodePhase, BodeOmega = control.bode_plot(self.G, dB = True)
        elif Units == 'Hz':
            BodeMagnitude, BodePhase, BodeOmega = control.bode_plot(self.G, Hz = True)
        elif Units == 'degree':
            BodeMagnitude, BodePhase, BodeOmega = control.bode_plot(self.G, deg = True)
        else:
            BodeMagnitude, BodePhase, BodeOmega = control.bode_plot(self.G)

        return BodeMagnitude, BodePhase, BodeOmega

    def Damping(self):
        import control.matlab

        Frequencies, Damping, Poles = control.matlab.damp(self.G)

        return Frequencies, Damping, Poles


def CompareResults(Sys1, Sys2, Sys3 = None, Sys4 = None, Sys5 = None, Sys6 = None, YUnit = None):

    def FindMaxNonZero(a, T):
        return T[int(max([i for i, e in enumerate(a) if e != 0]))]

    if Sys1.TimeDelay == None:
        Sys1.TimeDelay = 0

    import matplotlib.pyplot as plt
    plt.style.use('seaborn-darkgrid')

    LineWidth = 3

    """
    Figure 1: Plot individual systems wit inputs
    """
    plt.figure(figsize=(14,12))

    plt.tight_layout()

    h = 3

    plt.subplot(3,2,1)

    if Sys1.InputFunction == 'Square':
        plt.hlines(Sys1.Magnitude, 0, FindMaxNonZero(Sys1.U, Sys1.T), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys1.InputFunction, Sys1.Magnitude))

        plt.vlines(0, 0, Sys1.Magnitude)
        plt.vlines(FindMaxNonZero(Sys1.U, Sys1.T), 0, Sys1.Magnitude,  linewidth = LineWidth)
    
        color = 'g'
    
    elif Sys1.InputFunction == 'Impulse':
        plt.vlines(0, 0, Sys1.Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Sys1.Magnitude))

        color = 'r'

    else:
        plt.hlines(Sys1.Magnitude, 0, Sys1.Time[-1], linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys1.InputFunction, Sys1.Magnitude))
        try:
            plt.vlines(Sys1.T[0], 0, Sys1.Magnitude,  linewidth = LineWidth)
        except:
            plt.vlines(0, 0, Sys1.Magnitude,  linewidth = LineWidth)

        color = 'b'

    
    plt.plot(Sys1.Time, Sys1.Output, color = color, linewidth = LineWidth, label = 'TF 1 (Time Delay = {})'.format(Sys1.TimeDelay))

    plt.xlabel('Time')
    plt.ylabel(YUnit + ' (Deviation)')
    plt.legend(loc = 'best')

    
    plt.subplot(3,2,2)

    if Sys2.InputFunction == 'Square':
        plt.hlines(Sys2.Magnitude, 0, FindMaxNonZero(Sys2.U, Sys2.T), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys2.InputFunction, Sys2.Magnitude))

        plt.vlines(0, 0, Sys2.Magnitude)
        plt.vlines(FindMaxNonZero(Sys2.U, Sys2.T), 0, Sys2.Magnitude, linewidth = LineWidth)

        color = 'g'

    elif Sys2.InputFunction == 'Impulse':
        plt.vlines(0, 0, Sys2.Magnitude, linewidth = LineWidth)

        color = 'r'

    else:
        plt.hlines(Sys2.Magnitude, 0, Sys1.Time[-1], linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys2.InputFunction, Sys2.Magnitude))

        try:
            plt.vlines(Sys2.T[0], 0, Sys2.Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Sys2.Magnitude))
        except:
            plt.vlines(0, 0, Sys2.Magnitude, linewidth = LineWidth)     

        color = 'b'          

    
    plt.plot(Sys1.Time, Sys2.Output, color = color, linewidth = LineWidth,  label = 'TF 2 (Time Delay = {})'.format(Sys2.TimeDelay))

    plt.xlabel('Time')
    plt.ylabel(YUnit + ' (Deviation)')
    plt.legend(loc = 'best')


    if Sys3 is not None:
        plt.subplot(3,2,h)
        
        if Sys3.InputFunction == 'Square':
            plt.hlines(Sys3.Magnitude, 0, FindMaxNonZero(Sys3.T), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys3.InputFunction, Sys3.Magnitude))

            plt.vlines(0, 0, Sys3.Magnitude, linewidth = LineWidth)
            plt.vlines(FindMaxNonZero(Sys3.T), 0, Sys3.Magnitude, linewidth = LineWidth)

            color = 'g'

        elif Sys3.InputFunction == 'Impulse':
            plt.vlines(0, 0, Sys3.Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Sys3.Magnitude))

            color = 'r'
        
        else:
            try:
                plt.vlines(Sys3.T[0], 0, Sys3.Magnitude ,linewidth = LineWidth)
            except:
                plt.vlines(0, 0, Sys3.Magnitude, linewidth = LineWidth)

            color = 'b'

        plt.plot(Sys1.Time, Sys3.Output, color = color, linewidth = LineWidth,  label = 'TF 3 (Time Delay = {})'.format(Sys3.TimeDelay))

        plt.xlabel('Time')
        plt.ylabel(YUnit + ' (Deviation)')
        plt.legend(loc = 'best')

        h += 1

    if Sys4 is not None:
        plt.subplot(3,2,h)

        if Sys4.InputFunction == 'Square':
            plt.hlines(Sys4.Magnitude, 0, FindMaxNonZero(Sys4.T), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys4.InputFunction, Sys4.Magnitude))

            plt.vlines(0, 0, Sys4.Magnitude, linewidth = LineWidth)
            plt.vlines(FindMaxNonZero(Sys4.T), 0, Sys4.Magnitude, linewidth = LineWidth)

            color = 'g'

        elif Sys4.InputFunction == 'Impulse':
            plt.vlines(0, 0, Sys4.Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Sys5.Magnitude))

            color = 'r'
        
        else:
            plt.hlines(Sys4.Magnitude, 0, Sys1.Time[-1], linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys4.InputFunction, Sys4.Magnitude))

            try:
                plt.vlines(Sys4.T[0], 0, Sys4.Magnitude, linewidth = LineWidth)
            except:
                plt.vlines(0, 0, Sys4.Magnitude, linewidth = LineWidth)

            color = 'b'

        plt.plot(Sys1.Time, Sys4.Output, color = color, linewidth = LineWidth,  label = 'TF 4 (Time Delay = {})'.format(Sys4.TimeDelay))

        plt.xlabel('Time')
        plt.ylabel(YUnit + ' (Deviation)')
        plt.legend(loc = 'best')

        h += 1

    if Sys5 is not None:
        plt.subplot(3,2,h)

        if Sys5.InputFunction == 'Square':
            plt.hlines(Sys5.Magnitude, 0, FindMaxNonZero(Sys5.U, Sys5.T), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys5.InputFunction, Sys5.Magnitude))

            plt.vlines(0, 0, Sys5.Magnitude, linewidth = LineWidth)
            plt.vlines(FindMaxNonZero(Sys5.U, Sys5.T), 0, Sys5.Magnitude, linewidth = LineWidth)

            color = 'g'

        elif Sys5.InputFunction == 'Impulse':
            plt.vlines(0, 0, Sys5.Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Sys5.Magnitude))

            color = 'r'

        else:
            plt.hlines(Sys5.Magnitude, 0, Sys1.Time[-1], linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys5.InputFunction, Sys5.Magnitude))

            try:
                plt.vlines(Sys5.T[0], 0, Sys5.Magnitude, linewidth = LineWidth)
            except:
                plt.vlines(0, 0, Sys5.Magnitude, linewidth = LineWidth)     

            color = 'b'

        plt.plot(Sys1.Time, Sys5.Output, color = color, linewidth = LineWidth,  label = 'TF 5 (Time Delay = {})'.format(Sys5.TimeDelay))

        plt.xlabel('Time')
        plt.ylabel(YUnit + ' (Deviation)')
        plt.legend(loc = 'best')

        h += 1

    if Sys6 is not None:
        plt.subplot(3,2,h)

        if Sys6.InputFunction == 'Square':
            plt.hlines(Sys6.Magnitude, 0, FindMaxNonZero(Sys6.U, Sys6.T), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys6.InputFunction, Sys6.Magnitude))

            plt.vlines(0, 0, Sys6.Magnitude, linewidth = LineWidth)
            plt.vlines(FindMaxNonZero(Sys6.U, Sys6.T), 0, Sys6.Magnitude, linewidth = LineWidth)
            
            color = 'g'
        
        elif Sys6.InputFunction == 'Impulse':
            plt.vlines(0, 0, Sys6.Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Sys6.Magnitude))

            color = 'r'

        else:
            plt.hlines(Sys6.Magnitude, 0, Sys1.Time[-1], linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Sys6.InputFunction, Sys6.Magnitude))

            try:
                plt.vlines(Sys6.T[0], 0, Sys6.Magnitude, linewidth = LineWidth)
            except:
                plt.vlines(0, 0, Sys6.Magnitude, linewidth = LineWidth)       

            color = 'b'

        plt.plot(Sys1.Time, Sys6.Output, color = color, linewidth = LineWidth,  label = 'TF 6 (Time Delay = {})'.format(Sys6.TimeDelay))

        plt.xlabel('Time')
        plt.ylabel(YUnit + ' (Deviation)')
        plt.legend(loc = 'best')

        h += 1

    """
    Figure 2: Plot all systems at once
    """

    def SelectColor(Type):
        if Type == 'Square':
            color = 'g'

        elif Type == 'Impulse':
            color = 'r'

        elif Type == 'Step':
            color = 'b'

        return color

    StyleList = ['-', '--', '-o', '-x', '-+', '-*']
    ColorList = []

    plt.figure(figsize=(14,12))

    ColorList.append(SelectColor(Sys1.InputFunction))

    style = StyleList[ColorList.count(SelectColor(Sys1.InputFunction)) - 1]
    
    plt.plot(Sys1.Time, Sys1.Output, SelectColor(Sys1.InputFunction) + style, linewidth = LineWidth, \
        label =  'TF 1 ({} Input, Magnitude = {}, Time Delay = {})'.format(Sys1.InputFunction, Sys1.Magnitude, Sys1.TimeDelay))

    

    if Sys2 is not None:
        ColorList.append(SelectColor(Sys2.InputFunction))

        style = StyleList[ColorList.count(SelectColor(Sys2.InputFunction)) - 1]

        plt.plot(Sys1.Time, Sys2.Output, SelectColor(Sys2.InputFunction) + style, linewidth = LineWidth, \
            label = 'TF 2 ({} Input, Magnitude = {}, Time Delay = {})'.format(Sys2.InputFunction, Sys2.Magnitude, Sys2.TimeDelay))

    if Sys3 is not None:
        ColorList.append(SelectColor(Sys3.InputFunction))

        style = StyleList[ColorList.count(SelectColor(Sys3.InputFunction)) - 1]

        plt.plot(Sys1.Time, Sys3.Output,  SelectColor(Sys3.InputFunction) + style, linewidth = LineWidth, \
            label = 'TF 3 ({} Input, Magnitude = {}, Time Delay = {})'.format(Sys3.InputFunction, Sys3.Magnitude, Sys3.TimeDelay))

    if Sys4 is not None:
        ColorList.append(SelectColor(Sys4.InputFunction))

        style = StyleList[ColorList.count(SelectColor(Sys4.InputFunction)) - 1]

        plt.plot(Sys1.Time, Sys4.Output, SelectColor(Sys4.InputFunction) + style, linewidth = LineWidth, \
            label = 'TF 4 ({} Input, Magnitude = {}, Time Delay = {})'.format(Sys4.InputFunction, Sys4.Magnitude, Sys4.TimeDelay))

    if Sys5 is not None:
        ColorList.append(SelectColor(Sys5.InputFunction))

        style = StyleList[ColorList.count(SelectColor(Sys5.InputFunction)) - 1]

        plt.plot(Sys1.Time, Sys5.Output, SelectColor(Sys5.InputFunction) + style, linewidth = LineWidth, \
            label = 'TF 5 ({} Input, Magnitude = {}, Time Delay = {})'.format(Sys5.InputFunction, Sys5.Magnitude, Sys5.TimeDelay))

    if Sys6 is not None:
        ColorList.append(SelectColor(Sys6.InputFunction))

        style = StyleList[ColorList.count(SelectColor(Sys6.InputFunction)) - 1]

        plt.plot(Sys1.Time, Sys6.Output,  SelectColor(Sys6.InputFunction) + style, linewidth = LineWidth, \
            label = 'TF 6 ({} Input, Magnitude = {}, Time Delay = {})'.format(Sys6.InputFunction, Sys6.Magnitude, Sys6.TimeDelay))

    plt.legend()
    plt.xlabel('Time', fontsize = 14)

    if YUnit is not None:
        plt.ylabel(YUnit +  ' as Deviation Variable', fontsize = 14)

    plt.show()

    return


def FitFOPDT(Data, M):
    def Model(x):
        import numpy as np 
        K     = x[0] 
        tau   = x[1]
        theta = x[2]

        return [K*M*(1 - np.exp(-max(i - theta, 0)/tau)) for i in t]

    def ObjectiveFunction(x):
        return sum((Model(x) - Data))**2

    import numpy as np
    t = np.linspace(0, len(Data), len(Data))

    from scipy.optimize import minimize
    Guess = [1, 1, 1]
    solution = minimize(ObjectiveFunction, Guess, bounds = ((0, 10), (0, 10), (0, 10)))

    EstimatedParameters = solution.x

    print('K = {}'.format(EstimatedParameters[0]))
    print('tau = {}'.format(EstimatedParameters[1]))
    print('theta = {}'.format(EstimatedParameters[2]))

    OptimalModelOutput = Model(EstimatedParameters)

    import matplotlib.pyplot as plt

    plt.plot(Data, label = 'Process Data')
    plt.plot(OptimalModelOutput, label = 'Optimal Model Output')

    plt.legend()

    plt.show()

    return
  





