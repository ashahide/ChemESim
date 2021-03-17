class ProcessControlandModeling:

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

        self.ResponseMagnitude = Magnitude
        self.ResponseType      = Type
        self.Output            = Magnitude*Output
        self.Time              = Time
        self.InputTimeDefined  = T


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


    def PlotResults(self, Output2 = None, Input2Type = None, Input2Magnitude = None, Input2TimeDefined = None, Sys2TimeDelay = None, \
                          Output3 = None, Input3Type = None, Input3Magnitude = None, Input3TimeDefined = None, Sys3TimeDelay = None, \
                          Output4 = None, Input4Type = None, Input4Magnitude = None, Input4TimeDefined = None, Sys4TimeDelay = None, \
                          Output5 = None, Input5Type = None, Input5Magnitude = None, Input5TimeDefined = None, Sys5TimeDelay = None, \
                          Output6 = None, Input6Type = None, Input6Magnitude = None, Input6TimeDefined = None, Sys6TimeDelay = None, \
                          YUnit = None):

        def FindMaxNonZero(a):
            return a[int(max([i for i, e in enumerate(a) if e != 0]))]

        if self.TimeDelay == None:
            self.TimeDelay = 0
  
        import matplotlib.pyplot as plt
        plt.style.use('seaborn-darkgrid')

        LineWidth = 3

        """
        Figure 1: Plot individual systems wit inputs
        """
        plt.figure(figsize=(14,12))

        plt.tight_layout()

        h = 2

        plt.subplot(3,2,1)

        if self.InputFunction == 'Square':
            plt.hlines(self.ResponseMagnitude, 0, FindMaxNonZero(self.InputTimeDefined), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(self.ResponseType, self.ResponseMagnitude))

            plt.vlines(0, 0, self.ResponseMagnitude)
            plt.vlines(FindMaxNonZero(self.InputTimeDefined), 0, self.ResponseMagnitude,  linewidth = LineWidth)
        
            color = 'g'
        
        elif self.InputFunction == 'Impulse':
            plt.vlines(0, 0, self.ResponseMagnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(self.ResponseMagnitude))

            color = 'r'

        else:
            plt.hlines(self.ResponseMagnitude, 0, self.Time[-1], linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(self.ResponseType, self.ResponseMagnitude))
            try:
                plt.vlines(self.InputTimeDefined[0], 0, self.ResponseMagnitude,  linewidth = LineWidth)
            except:
                plt.vlines(0, 0, self.ResponseMagnitude,  linewidth = LineWidth)

            color = 'b'

        
        plt.plot(self.Time, self.Output, color = color, linewidth = LineWidth, label = 'TF 1 (Time Delay = {})'.format(self.TimeDelay))

        plt.xlabel('Time')
        plt.ylabel(YUnit + ' (Deviation)')
        plt.legend(loc = 'best')

        if Output2 is not None:
            plt.subplot(3,2,h)

            if Input2Type == 'Square':
                plt.hlines(Input2Magnitude, 0, FindMaxNonZero(Input2TimeDefined), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Input2Type, Input2Magnitude))

                plt.vlines(0, 0, Input2Magnitude)
                plt.vlines(FindMaxNonZero(Input2TimeDefined), 0, Input2Magnitude, linewidth = LineWidth)

                color = 'g'

            elif Input2Type == 'Impulse':
                plt.vlines(0, 0, Input2Magnitude, linewidth = LineWidth)

                color = 'r'

            else:
                plt.hlines(Input2Magnitude, 0, self.Time[-1], linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Input2Type, Input2Magnitude))

                try:
                    plt.vlines(Input2TimeDefined[0], 0, Input2Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Input2Magnitude))
                except:
                    plt.vlines(0, 0, Input2Magnitude, linewidth = LineWidth)     

                color = 'b'          

            
            plt.plot(self.Time, Output2, color = color, linewidth = LineWidth,  label = 'TF 2 (Time Delay = {})'.format(Sys2TimeDelay))

            plt.xlabel('Time')
            plt.ylabel(YUnit + ' (Deviation)')
            plt.legend(loc = 'best')

            h += 1

        if Output3 is not None:
            plt.subplot(3,2,h)
            
            if Input3Type == 'Square':
                plt.hlines(Input3Magnitude, 0, FindMaxNonZero(Input3TimeDefined), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Input3Type, Input3Magnitude))

                plt.vlines(0, 0, Input3Magnitude, linewidth = LineWidth)
                plt.vlines(FindMaxNonZero(Input3TimeDefined), 0, Input3Magnitude, linewidth = LineWidth)

                color = 'g'

            elif Input3Type == 'Impulse':
                plt.vlines(0, 0, Input3Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Input3Magnitude))

                color = 'r'
            
            else:
                try:
                    plt.vlines(Input3TimeDefined[0], 0, Input3Magnitude ,linewidth = LineWidth)
                except:
                    plt.vlines(0, 0, Input3Magnitude, linewidth = LineWidth)

                color = 'b'

            plt.plot(self.Time, Output3, color = color, linewidth = LineWidth,  label = 'TF 3 (Time Delay = {})'.format(Sys3TimeDelay))

            plt.xlabel('Time')
            plt.ylabel(YUnit + ' (Deviation)')
            plt.legend(loc = 'best')

            h += 1

        if Output4 is not None:
            plt.subplot(3,2,h)

            if Input4Type == 'Square':
                plt.hlines(Input4Magnitude, 0, FindMaxNonZero(Input4TimeDefined), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Input4Type, Input4Magnitude))

                plt.vlines(0, 0, Input4Magnitude, linewidth = LineWidth)
                plt.vlines(FindMaxNonZero(Input4TimeDefined), 0, Input4Magnitude, linewidth = LineWidth)

                color = 'g'

            elif Input4Type == 'Impulse':
                plt.vlines(0, 0, Input4Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Input5Magnitude))

                color = 'r'
            
            else:
                plt.hlines(Input4Magnitude, 0, self.Time[-1], linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Input4Type, Input4Magnitude))

                try:
                    plt.vlines(Input4TimeDefined[0], 0, Input4Magnitude, linewidth = LineWidth)
                except:
                    plt.vlines(0, 0, Input4Magnitude, linewidth = LineWidth)

                color = 'b'

            plt.plot(self.Time, Output4, color = color, linewidth = LineWidth,  label = 'TF 4 (Time Delay = {})'.format(Sys4TimeDelay))

            plt.xlabel('Time')
            plt.ylabel(YUnit + ' (Deviation)')
            plt.legend(loc = 'best')

            h += 1

        if Output5 is not None:
            plt.subplot(3,2,h)

            if Input5Type == 'Square':
                plt.hlines(Input5Magnitude, 0, FindMaxNonZero(Input5TimeDefined), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Input5Type, Input5Magnitude))

                plt.vlines(0, 0, Input5Magnitude, linewidth = LineWidth)
                plt.vlines(FindMaxNonZero(Input5TimeDefined), 0, Input5Magnitude, linewidth = LineWidth)

                color = 'g'

            elif Input5Type == 'Impulse':
                plt.vlines(0, 0, Input5Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Input5Magnitude))

                color = 'r'

            else:
                plt.hlines(Input5Magnitude, 0, self.Time[-1], linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Input5Type, Input5Magnitude))

                try:
                    plt.vlines(Input5TimeDefined[0], 0, Input5Magnitude, linewidth = LineWidth)
                except:
                    plt.vlines(0, 0, Input5Magnitude, linewidth = LineWidth)     

                color = 'b'

            plt.plot(self.Time, Output5, color = color, linewidth = LineWidth,  label = 'TF 5 (Time Delay = {})'.format(Sys5TimeDelay))

            plt.xlabel('Time')
            plt.ylabel(YUnit + ' (Deviation)')
            plt.legend(loc = 'best')

            h += 1

        if Output6 is not None:
            plt.subplot(3,2,h)

            if Input6Type == 'Square':
                plt.hlines(Input6Magnitude, 0, FindMaxNonZero(Input6TimeDefined), linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Input6Type, Input6Magnitude))

                plt.vlines(0, 0, Input6Magnitude, linewidth = LineWidth)
                plt.vlines(FindMaxNonZero(Input6TimeDefined), 0, Input6Magnitude, linewidth = LineWidth)
                
                color = 'g'
            
            elif Input6Type == 'Impulse':
                plt.vlines(0, 0, Input6Magnitude, linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Input6Magnitude))

                color = 'r'

            else:
                plt.hlines(Input6Magnitude, 0, self.Time[-1], linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Input6Type, Input6Magnitude))

                try:
                    plt.vlines(Input6TimeDefined[0], 0, Input6Magnitude, linewidth = LineWidth)
                except:
                    plt.vlines(0, 0, Input6Magnitude, linewidth = LineWidth)       

                color = 'b'

            plt.plot(self.Time, Output6, color = color, linewidth = LineWidth,  label = 'TF 6 (Time Delay = {})'.format(Sys6TimeDelay))

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

        ColorList.append(SelectColor(self.ResponseType))

        style = StyleList[ColorList.count(SelectColor(self.ResponseType)) - 1]
        
        plt.plot(self.Time, self.Output, SelectColor(self.ResponseType) + style, linewidth = LineWidth, label =  'TF 1 ({} Input, Magnitude = {}, Time Delay = {})'.format(self.ResponseType, self.ResponseMagnitude, self.TimeDelay))

        

        if Output2 is not None:
            ColorList.append(SelectColor(Input2Type))

            style = StyleList[ColorList.count(SelectColor(Input2Type)) - 1]

            plt.plot(self.Time, Output2, SelectColor(Input2Type) + style, linewidth = LineWidth, label = 'TF 2 ({} Input, Magnitude = {}, Time Delay = {})'.format(Input2Type, Input2Magnitude, Sys2TimeDelay))

        if Output3 is not None:
            ColorList.append(SelectColor(Input3Type))

            style = StyleList[ColorList.count(SelectColor(Input3Type)) - 1]

            plt.plot(self.Time, Output3,  SelectColor(Input3Type) + style, linewidth = LineWidth, label = 'TF 3 ({} Input, Magnitude = {}, Time Delay = {})'.format(Input3Type, Input3Magnitude, Sys3TimeDelay))

        if Output4 is not None:
            ColorList.append(SelectColor(Input4Type))

            style = StyleList[ColorList.count(SelectColor(Input4Type)) - 1]

            plt.plot(self.Time, Output4, SelectColor(Input4Type) + style, linewidth = LineWidth, label = 'TF 4 ({} Input, Magnitude = {}, Time Delay = {})'.format(Input4Type, Input4Magnitude, Sys4TimeDelay))

        if Output5 is not None:
            ColorList.append(SelectColor(Input5Type))

            style = StyleList[ColorList.count(SelectColor(Input5Type)) - 1]

            plt.plot(self.Time, Output5, SelectColor(Input5Type) + style, linewidth = LineWidth, label = 'TF 5 ({} Input, Magnitude = {}, Time Delay = {})'.format(Input5Type, Input5Magnitude, Sys5TimeDelay))

        if Output6 is not None:
            ColorList.append(SelectColor(Input6Type))

            style = StyleList[ColorList.count(SelectColor(Input6Type)) - 1]

            plt.plot(self.Time, Output6,  SelectColor(Input6Type) + style, linewidth = LineWidth, label = 'TF 6 ({} Input, Magnitude = {}, Time Delay = {})'.format(Input6Type, Input6Magnitude, Sys6TimeDelay))

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
  





