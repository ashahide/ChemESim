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


    def PlotResults(self, Output2 = None, Input2Type = None, Input2Magnitude = None, Sys2TimeDelay = None, \
                          Output3 = None, Input3Type = None, Input3Magnitude = None, Sys3TimeDelay = None, \
                          Output4 = None, Input4Type = None, Input4Magnitude = None, Sys4TimeDelay = None, \
                          YUnit = None):

        if self.TimeDelay == None:
            self.TimeDelay = 0

        import matplotlib.pyplot as plt

        plt.figure(figsize=(14,12))
        plt.plot(self.Time, self.Output, label = 'TF 1 ({} Input, Magnitude = {}, Time Delay = {})'.format(self.ResponseType, self.ResponseMagnitude, self.TimeDelay))

        if Output2 is not None:
            plt.plot(self.Time, Output2, label = 'TF 2 ({} Input, Magnitude = {}, Time Delay = {})'.format(Input2Type, Input2Magnitude, Sys2TimeDelay))

        if Output3 is not None:
            plt.plot(self.Time, Output3, label = 'TF 3 ({} Input, Magnitude = {}, Time Delay = {})'.format(Input3Type, Input3Magnitude, Sys3TimeDelay))

        if Output4 is not None:
            plt.plot(self.Time, Output4, label = 'TF 4 ({} Input, Magnitude = {}, Time Delay = {})'.format(Input4Type, Input4Magnitude, Sys4TimeDelay))

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

        return [K*M*(1 - np.exp(-np.max((i - theta), 0)/tau)) for i in t]

    def ObjectiveFunction(x):
        return sum((Model(x) - Data))**2

    import numpy as np
    t = np.linspace(0, len(Data), len(Data))

    from scipy.optimize import minimize
    Guess = [1, 2, 5]
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
  





