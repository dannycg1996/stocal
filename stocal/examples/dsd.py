"""Exact autocatalytic polymer model

This model implements the system published in

    H Fellermann, S Tanaka, R Rasmussen, Sequence selection by
    dynamical symmetry breaking in an autocatalytic binary polymer
    model, Phys. Rev. E 96, 062407 (2017)
    https://doi.org/10.1103/PhysRevE.96.062407

The model is based on three rules that determine the behavior of
self-replicating heteropolymers:

  * The DegradationRule allows a polymer to break into two
    substrands (k.l --> k + l)

  * The LigationRule allows two random polymers to form a
    non-commutative ligation product (k+ l --> k.l)

  * The AutoCatalysisRule allows a polymer to be replicated
    catalytically from two exactly matching substrands
    (k + l + k.l --> 2*k.l)

The rules operate over any monomer alphabet, but the model defines
an initial state with two types of monomers, making this a binary exact
autocatalytic polymer model.

With the default parameters, simulations will take a significant
amount of time.
"""
import stocal
import re

alpha = 1.e-10
beta = 1000**-2

initial_state = {"{S' N^* L' R'}": 60, "<L N^ M N^>": 50, "<N^* L' R'>": 60, "{L N^ M N^ R}": 50}
format_check = re.compile('<[^\[\]{}]*?>|\[[^<>{}]*?]|{[^<>\[\]]*?}|\s*:*')

upper_th = re.compile('(\w)(?=\^\s)(?=[^<>]*>)|(\w)(?=\^>)(?=[^<>]*>)')
upper_th_c = re.compile('(\w)(?=\^\*)(?=[^<>]*>)')
lower_th = re.compile('(\w)(?=\^\s)(?=[^{}]*})|(\w)(?=\^})(?=[^{}]*})')
lower_th_c = re.compile('(\w)(?=\^\*)(?=[^{}]*})')
empty_bracket = re.compile('(<(?:\s)*>)|({(?:\s)*})')

class BindingRule(stocal.TransitionRule):
    """Join any two strings into their concatenations"""
    Transition = stocal.MassAction

    def novel_reactions(self, k, l):
        print(k)
        print(l)
        yield from self.toehold_binding(k, l, upper_th, lower_th_c)
        yield from self.toehold_binding(k, l, lower_th, upper_th_c)
        yield from self.toehold_binding(k, l, upper_th_c, lower_th)
        yield from self.toehold_binding(k, l, lower_th_c, upper_th)

    def toehold_binding(self, k, l, regex_1, regex_2):
        for matching_1 in re.finditer(regex_1, k):
            for matching_2 in re.finditer(regex_2, l):
                if matching_1.group() == matching_2.group():
                    if regex_1 == upper_th or regex_1 == upper_th_c:
                        print("upper_th")
                        part_A = k[:matching_1.start()] + ">" + l[:matching_2.start()] + "}"
                        part_B ="{" + l[matching_2.end()+2:] + "<" + k[matching_1.end()+2:]
                    elif regex_1 == lower_th or regex_1 == lower_th_c:
                        print("lower")
                        part_A = k[:matching_1.start()] + "}" + l[:matching_2.start()] + ">"
                        part_B = "<" + l[matching_2.end()+2:] + "{" + k[matching_1.end()+2:]
                    draft_strand = part_A + "[" + matching_2.group() + "^]" + part_B
                    final_strand = re.sub(empty_bracket, '', draft_strand)
                    print(final_strand)
                    yield self.Transition([k, l], [final_strand], alpha)

process = stocal.Process(
    rules=[BindingRule()]
)

if __name__ == '__main__':
    format_correct = True
    for key in initial_state:
        if re.fullmatch(format_check, key) is None:
            format_correct = False
            print("Incorrect Input") #Some of the brackets must not match up. Raise a proper error in time.

    if format_correct:
        traj = process.sample(initial_state, tmax=100.)
        for _ in traj:
            print(traj.time, traj.state)
