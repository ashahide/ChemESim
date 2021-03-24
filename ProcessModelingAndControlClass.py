class TransferFunction:

    """
    Links to cite:
    [1] Python control toolbox from CalTech
    [2] http://techteach.no/python_control/python_control.pdf
    """

    def __init__(self, Numerator, Denominator, TimeDelay = None, N = None):

        self.Numerator   = Numerator
        self.Denominator = Denominator
        self.TimeDelay   = TimeDelay
        self.N           = N #for pade approximation [2]

        """
        Using method for time delay creation from [2]
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
            self.Time              = Time

        elif Type == 'Impulse':
            Output, Time = control.matlab.impulse(self.G, t)
            self.Time              = Time

        elif Type == 'Square':
            T, Output  = control.forced_response(self.G, T=T, U=U)
            self.Time              = T

        self.Magnitude = Magnitude
        self.Input      = Type
        self.Output            = Magnitude*Output
  
        self.T  = T
        self.U  = U


        return Magnitude*Output, self.Time

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

    def PlotResponse(self):
        def FindMaxNonZero(a, T):
            return T[int(max([i for i, e in enumerate(a) if e != 0]))]

        import matplotlib.pyplot as plt 
        plt.style.use('seaborn-darkgrid')

        if self.Input == 'Square':
            color = 'g'
        
        elif self.Input == 'Impulse':
            color = 'r'

        else:
            color = 'b'

        LineWidth = 3

        fig, ax = plt.subplots()

        ax.plot(self.Output, color = color, linewidth = LineWidth, label = "TF")
        ax.set_ylabel("Output", color = color, fontsize = 14)
        ax.tick_params(axis = 'y', labelcolor = color, labelsize = 12)
        ax.set_xlabel('Time', fontsize = 14)

        ax2 = ax.twinx()
        ax2.set_ylabel("Input", color = 'k', fontsize = 14)
        ax2.tick_params(axis = 'y', labelcolor = 'k', labelsize = 12)


        if self.Input == 'Square':
            ax2.hlines(self.Magnitude, 0, FindMaxNonZero(self.U, self.T), color = 'k', linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(self.Input, self.Magnitude))

            ax2.vlines(0, 0, self.Magnitude, color = 'k',  linewidth = LineWidth)
            ax2.vlines(FindMaxNonZero(self.U, self.T), 0, self.Magnitude, color = 'k',  linewidth = LineWidth)
        
        elif self.Input == 'Impulse':
            ax2.vlines(0, 0, self.Magnitude, linewidth = LineWidth, color = 'k', label = 'Impulse Input, Magnitude = {}'.format(self.Magnitude))

        else:
            ax2.hlines(self.Magnitude, 0, self.Time[-1], color = 'k', linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(self.Input, self.Magnitude))
            try:
                ax2.vlines(self.T[0], 0, self.Magnitude, color = 'k',  linewidth = LineWidth)
            except:
                ax2.vlines(0, 0, self.Magnitude, color = 'k',  linewidth = LineWidth)

        plt.show()

def CompareResults(*Systems, YUnit = None):

    def FindMaxNonZero(a, T):
        return T[int(max([i for i, e in enumerate(a) if e != 0]))]

    import matplotlib.pyplot as plt
    plt.style.use('seaborn-darkgrid')

    LineWidth = 3

    """
    Figure 1: Plot individual systems wit inputs
    """
    plt.figure(figsize=(14,12))

    plt.tight_layout()

    for i in range(len(Systems)):

        plt.subplot(3,2,i+1)

        if Systems[i].Input == 'Square':
            plt.hlines(Systems[i].Magnitude, 0, FindMaxNonZero(Systems[i].U, Systems[i].T), color = 'k', linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Systems[i].Input, Systems[i].Magnitude))
            plt.vlines(0, 0, Systems[i].Magnitude, linewidth = LineWidth, color = 'k')
            plt.vlines(FindMaxNonZero(Systems[i].U, Systems[i].T), 0, Systems[i].Magnitude, color = 'k',  linewidth = LineWidth)

            color = 'g'


        elif Systems[i].Input == 'Impulse':
            plt.vlines(0, 0, Systems[i].Magnitude, color = 'k', linewidth = LineWidth, label = 'Impulse Input, Magnitude = {}'.format(Systems[i].Magnitude))

            color = 'r'

        else:
            plt.hlines(Systems[i].Magnitude, 0, Systems[i].Time[-1], color = 'k', linewidth = LineWidth, label = '{} Input, Magnitude = {}'.format(Systems[i].Input, Systems[i].Magnitude))
            try:
                plt.vlines(Systems[i].T[0], 0, Systems[i].Magnitude, color = 'k',  linewidth = LineWidth)
            except:
                plt.vlines(0, 0, Systems[i].Magnitude, color = 'k',  linewidth = LineWidth)

            color = 'b'

        plt.plot(Systems[i].Time, Systems[i].Output, color = color, linewidth = LineWidth, label = 'TF {} (Time Delay = {})'.format(i+1, Systems[i].TimeDelay))

        plt.xlabel('Time')
        plt.ylabel(YUnit + ' (Deviation)')
        plt.legend(loc = 'best')

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

    for i in range(len(Systems)):
        ColorList.append(SelectColor(Systems[i].Input))

        style = StyleList[ColorList.count(SelectColor(Systems[i].Input)) - 1]
    
        plt.plot(Systems[i].Time, Systems[i].Output, SelectColor(Systems[i].Input) + style, linewidth = LineWidth, \
            label =  f'TF {i+1} ({Systems[i].Input} Input, Magnitude = {Systems[i].Magnitude}, Time Delay = {Systems[i].TimeDelay})')

    
        plt.legend()
        plt.xlabel('Time', fontsize = 14)

        if YUnit is not None:
            plt.ylabel(YUnit +  ' as Deviation Variable', fontsize = 14)

    plt.show()

