class ReactionKinetics:
    """
    References: 
    [1] https://stackoverflow.com/questions/49896768/how-can-i-use-python-to-create-a-stoichiometric-matrix user: WolfgangK

    """

    def __init__(self, ReactionPhase):
        self.ReactionPhase = ReactionPhase

    def DefineReactions(self, ChemicalReactions):
        """
        Code so far copied directly from [1] (as of 03/23/2021)

        ChemicalReactions must be a single list in the form of ['R1: ', 'R2: ']

        Ex. ['R1: A + 2B + C <=> D', 'R2: A + B <=> C']

        """
        import pandas as pd
        import re  # regular expression
        
        def coeff_comp(s):
            # Separate stoichiometric coefficient and compound
            result = re.search('(?P<coeff>\d*)(?P<comp>.*)', s)
            coeff = result.group('coeff')
            comp = result.group('comp')
            if not coeff:
                coeff = 0                          # coefficient=0 if it is missing
            return comp, int(coeff)

        reactions_dict = {}                          # results dictionary

        for equation in ChemicalReactions:
            compounds = {}                           # dict -> compound: coeff 
            eq = equation.replace(' ', '')  
            r_id, reaction = eq.split(':')           # separate id from chem reaction
            lhs, rhs = reaction.split('<=>')         # split left and right hand side

            reagents = lhs.split('+')                # get list of reagents
            products = rhs.split('+')                # get list of products
            for reagent in reagents:
                comp, coeff = coeff_comp(reagent)
                compounds[comp] = - coeff            # negative on lhs
            for product in products:
                comp, coeff = coeff_comp(product)
                compounds[comp] = coeff              # positive on rhs
            reactions_dict[r_id] = compounds         

        self.ChemicalReactions    = ChemicalReactions
        self.StoichiometricMatrix = pd.DataFrame(reactions_dict).fillna(value=0).astype(int)
    
    def DefineReactorSystem(self, ReactorType, Inlets, Outlets):
        """
        Define system of reactors (single reactors or linked)

        ReactorType = List containing Batch, CSTR, PFR, or PBR
                        Ex. ReactorType = ['CSTR', 'PFR]

        Inlets  = List containing a separate list for each reactor. Each species in each list represents molar inlets for each
                that species to the reactor
                Ex. Inlets = [['A', 'B'], ['All']]

        Outlets = List containing a separate list for each reactor. Each species in each list represents molar outlets for each
                that species to the reactor. 'All' means every species should have some outlet (even if it's totally converted)
                Ex. Outlets = [['All'], ['All']]
        """
        import pandas as pd
        import numpy as np

        def FormatInletOutletMatrices(FlowList, Reactors, Type):
            FlowMatrix = pd.DataFrame(index = self.StoichiometricMatrix.index, columns = Reactors)

            for i in range(len(FlowList)):
                for j in range(len(FlowList[i])):
                    if FlowMatrix.index[j] == FlowList[i][j]:
                        FlowMatrix[FlowMatrix.columns[i]][FlowMatrix.index[j]] = 1
                    elif FlowList[i][j] == 'All':
                        FlowMatrix[FlowMatrix.columns[i]][:] = 1
                    else:
                        FlowMatrix[FlowMatrix.columns[i]][FlowMatrix.index[j]] = 0
    
            return FlowMatrix.fillna(value = 0)

        NumReactors = len(ReactorType)
        Batch = 1
        CSTR  = 1
        PFR   = 1
        PBR   = 1

        Reactors = []

        for i in range(NumReactors):
            if ReactorType[i] == 'Batch':
                Reactors.append(f"Batch {Batch}")
                Batch += 1

            elif ReactorType[i] == 'CSTR':
                Reactors.append(f"CSTR {CSTR}")
                CSTR += 1

            elif ReactorType[i] == 'PFR':
                Reactors.append(f"PFR {PFR}")
                PFR += 1

            elif ReactorType[i] == 'PBR':
                Reactors.append(f"PBR {PBR}")
                PBR += 1

        InletMatrix  = FormatInletOutletMatrices(Inlets, Reactors, Type = 'Inlet')
        OutletMatrix = FormatInletOutletMatrices(Outlets, Reactors, Type = 'Outlet')

        self.NumReactors         = NumReactors
        self.ReactorArrangement  = Reactors
        self.InletMatrix         = InletMatrix
        self.OutletMatrix        = OutletMatrix
        
    def SimulateEquations(self):
        """
        Simulate equations with state space format where inlet flow rates are inputs 

        Xdot = A*X + B*U
            Xdot = Derivative of state variables
            A    = State coefficients
            X    = State vector
            B    = Input coefficients
            U    = Inputs
        """

        return



Sim = ReactionKinetics(ReactionPhase = 'Liquid')

Sim.DefineReactions(['R1: A + 2B + C <=> D', 'R2: A + B <=> C'])

Sim.DefineReactorSystem(ReactorType= ['CSTR', 'PFR'], Inlets = [['A', 'B'], ['All']], Outlets = [['All'], ['All']])

print(Sim.ReactorArrangement)
print(Sim.InletMatrix)
print(Sim.OutletMatrix)

