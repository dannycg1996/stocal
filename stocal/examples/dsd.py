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

#initial_state = {"{S' N^* L' R'}": 60, "<L N^ M N^>": 50}
#initial_state = {"{ N^* L' R'}": 60, "<L N^ M N^>": 50}
#initial_state = {"{N^* S' N^*}[C^]": 60, "<N^ M N^>": 50, "{L'}<L>[N^]<R>[M^]{P'}<S'>[A^]{B}" : 50}
initial_state = {"{L'}<L>[N^]<R>{R'}" : 50}
format_check = re.compile('<[^\[\]{}]*?>|\[[^<>{}]*?]|{[^<>\[\]]*?}|\s*:*')

#bound_th = re.compile('(\w)(?=\^)(?=[^\[\]]*])')
#double_toehold = re.compile('\[\W*?(\w)\^\W*?\]')
double_toehold = re.compile('(?:\[\W*?(\w)(?:\^\W*?\]))')
double_th_label = re.compile('(\w)(?=\^)(?=[^<>{}]*])')
upper_th = re.compile('(\w)(?=\^\s)(?=[^<>]*>)|(\w)(?=\^>)(?=[^<>]*>)')
lower_th_c = re.compile('(\w)(?=\^\*)(?=[^{}]*})')
start = re.compile('[\[{]*?(<)|[\[<]*?({)|[<{]*?(])')
end = re.compile('[\]}]*?(>)|[\]>]*?(})|[>}]*?(])')
empty_bracket = re.compile('(<(?:\s)*>)|({(?:\s)*})')

#last_upper_th = re.compile('<([^>]*)>*(?!(?:.*<|.*\[))')
last_upper_th = re.compile('<([^>]*)>')
last_lower_th = re.compile('{([^}]*)}')
#last_lower_th = re.compile('({[^}]*)}*(?!(?:.*{|.*\[))')
next_upper_th = re.compile('<([^>]*)>*')
next_lower_th = re.compile('{([^}]*)}*')

open_spaces = re.compile('<(\s)+?\w|{(\s)+?\w|\[(\s)+?\w')
close_spaces = re.compile('(?:\w|\^|\*)(\s)+?>|(?:\w|\^|\*)(\s)+?\]|(?:\w|\^|\*)(\s)+?}')

class BindingRule(stocal.TransitionRule):
    """Join any two strings into their concatenations"""
    Transition = stocal.MassAction

    #TODO: Stop this binding rule from allowing any binding on to toeholds on overhangs.
    #Either don't match on them (how?) or regex check for toeholds on overhangs in if __name__ == '__main__':

    def novel_reactions(self, k, l):
        #print("k:", k)
        #print("l:", l)
        #yield from self.toehold_binding(k, l, upper_th, lower_th_c)
        #yield from self.toehold_binding(k, l, lower_th_c, upper_th)
        yield from self.toehold_unbinding(k)
        yield from self.toehold_unbinding(l)

    def toehold_binding(self, k, l, regex_1, regex_2):
        for matching_1 in re.finditer(regex_1, k):
            for matching_2 in re.finditer(regex_2, l):
                if matching_1.group() == matching_2.group():
                    if regex_1 == upper_th:
                        print("upper_th")
                        part_A = k[:matching_1.start()] + re.search(end,k[matching_1.start():]).group() + l[:matching_2.start()] + re.search(end,l[matching_2.start():]).group()
                        part_B = re.search(start, l[:matching_2.end()]).group() + l[matching_2.end()+2:] + re.search(start, k[:matching_1.end()+1]).group() + k[matching_1.end()+1:]
                    else:
                        print("lower_th")
                        part_A = k[:matching_1.start()] + re.search(end,k[matching_1.start():]).group() + l[:matching_2.start()] + re.search(end,l[matching_2.start():]).group()
                        part_B = re.search(start, l[:matching_2.end()]).group() + l[matching_2.end()+1:] + re.search(start, k[:matching_1.end()+1]).group() + k[matching_1.end()+2:]

                    draft_strand = part_A + "[" + matching_2.group() + "^]" + part_B
                    final_strand = re.sub(empty_bracket, '', draft_strand)
                    print("final", final_strand)
                    yield self.Transition([k, l], [final_strand], alpha)

    def toehold_unbinding(self, k):
        k = re.sub(open_spaces, '', k)
        k = re.sub(close_spaces, '', k)
        print('k:', k)
        for double_th in re.finditer(double_toehold, k):
            print(k[:double_th.start()])
            label = re.search(double_th_label, double_th.group()).group()
            prefix, suffix, upper_1, lower_1, upper_2, lower_2 = "", "", "", "", "", ""
            bracket_open = k[:double_th.start()].rfind(']')
            bracket_close = k[double_th.end():].find('[')
            if bracket_open !=-1:
                prefix = k[:bracket_open]
            else:
                bracket_open = 0

            if bracket_close != -1:
                suffix = k[bracket_close:]
            else:
                bracket_close = len(double_th.group())

            if re.search(last_upper_th, k[bracket_open:double_th.start()]) is not None:
                upper_1 = re.search(last_upper_th, k[bracket_open:double_th.start()]).group()
            if re.search(last_lower_th, k[bracket_open:double_th.start()]) is not None:
                lower_1 = re.search(last_lower_th, k[bracket_open:double_th.start()]).group()
            if re.search(next_upper_th, k[double_th.end()+1:bracket_close]) is not None:
                upper_2 = re.search(next_upper_th, k[double_th.end()+1:bracket_close]).group()
            if re.search(last_lower_th, k[double_th.end():bracket_close]) is not None:
                lower_2 = re.search(next_lower_th, k[double_th.end():bracket_close]).group()

            #TODO: STOP the matchings (upper_th next etc) from including the brackets!!!! HOW?!
            print("prefix" ,prefix)
            print("search range", re.search(last_upper_th, k[bracket_open:double_th.start()]))
            print(bracket_open)
            print(double_th.start())
            print("upper 1", upper_1)
            print("label", label)
            print("upper 2", upper_2)
            print("lower 1", lower_1)
            print("lower 2", lower_2)
            print("suffix", suffix)
            part_A = prefix + "<" + upper_1[] + " " + label + "^ " + upper_2 + ">"
            part_A_fin = re.sub(empty_bracket, '', part_A)
            part_B = suffix + "{" + lower_1 + " " + label + "^*" + lower_2 + "}"
            part_B_fin = re.sub(empty_bracket, '', part_B)
            print("Final A:    ", part_A_fin)
            print("Final B:    ", part_B_fin)
            yield self.Transition([k], [part_A_fin, part_B_fin], alpha)


process = stocal.Process(
    rules=[BindingRule()]
)

if __name__ == '__main__':
    format_correct = True
    #TODO:  re.fullmatch() does not work as I thought it did, so error checking needs to be updated so as to make sure every <, { and [ has a corresponding >, }, ].
    #for key in initial_state:

        # if re.fullmatch(format_check, key) is None:
        #     format_correct = False
        #     print("Incorrect Input:", key) #Some of the brackets must not match up. Raise a proper error in time.

    #if format_correct:
    traj = process.sample(initial_state, tmax=100.)
    for _ in traj:
        print(traj.time, traj.state)

#Extension (Allow toehold complements on the upper strand and toeholds on the lower strand)
#upper_th_c = re.compile('(\w)(?=\^\*)(?=[^<>]*>)')
#lower_th = re.compile('(\w)(?=\^\s)(?=[^{}]*})|(\w)(?=\^})(?=[^{}]*})')
#         yield from self.toehold_binding(k, l, lower_th, upper_th_c)
#         yield from self.toehold_binding(k, l, upper_th_c, lower_th)

    # def toehold_binding(self, k, l, regex_1, regex_2):
    #     for matching_1 in re.finditer(regex_1, k):
    #         for matching_2 in re.finditer(regex_2, l):
    #             if matching_1.group() == matching_2.group():
    #                 if regex_1 == upper_th or regex_1 == upper_th_c:
    #                     print("upper_th")
    #                     part_A = k[:matching_1.start()] + ">" + l[:matching_2.start()] + "}"
    #                     part_B ="{" + l[matching_2.end()+2:] + "<" + k[matching_1.end()+2:]
    #                 elif regex_1 == lower_th or regex_1 == lower_th_c:
    #                     print("lower")
    #                     part_A = k[:matching_1.start()] + "}" + l[:matching_2.start()] + ">"
    #                     part_B = "<" + l[matching_2.end()+2:] + "{" + k[matching_1.end()+2:]
    #                 draft_strand = part_A + "[" + matching_2.group() + "^]" + part_B
    #                 final_strand = re.sub(empty_bracket, '', draft_strand)
    #                 print(final_strand)
    #                 yield self.Transition([k, l], [final_strand], alpha)


     #     print(re.findall(double_toehold_label,matching.group()))
        #     print("experiment", matching)
        # if len(re.findall(double_toehold, k)) == 1 and k.count('[') == 1:
        #     for matching in re.finditer(double_toehold, k):
        #         print(matching.group())
        #     print("findall", re.search(double_toehold, k))
        #     print("findall2", re.search(double_toehold, k).group())
        #     print(len(re.findall(double_toehold, k)))
        #    # print("finditer", re.finditer(double_toehold, k).group())
        #     yield self.Transition([k], [k], alpha)
