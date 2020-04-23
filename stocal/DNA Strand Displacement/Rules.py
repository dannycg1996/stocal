from stocal import *


class SplitToehold(TransitionRule) :
    TransitionRule = MassAction

    def novel_reactions(self, strand) :
        yield self.Transition()


class Dilution(TransitionRule) :
    Transition = MassAction

    def novel_reactions(self, species) :
        yield self.Transition([species], [], 0.001)


class Polymerization(TransitionRule) :
    Transition = MassAction

    def novel_reactions(self, k, l) :
        yield self.Transition([k,l], [k+l], 10.)


class Hydrolysis(TransitionRule) :
    Transition = MassAction

    def novel_reactions(self, k) :
        for i in range(1, len(k)) :
            c = 10.*i*(len(k)-i)
            yield self.Transition([k], [k[:i], k[i:]], c)