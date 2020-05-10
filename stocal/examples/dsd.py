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

format = re.compile('<[^\[\]{}]*?>|\[[^<>{}]*?]|{[^<>\[\]]*?}|\s*:*') #TODO: Fix this - speak to Harold.
# Needs to check that each part of a strand meets one of these requirements. Fullmatch does not work.

double_toehold = re.compile('(?:\[\W*?(\w)(?:\^\W*?\]))') #Matches on double toeholds of the form [A^] not [A^ B]
double_label = re.compile('(\w)(?=\^)(?=[^<>{}]*])') #Returns the label of a double toehold regex. Better way of doing this?
upper_label = re.compile('(\w)(?=\^\s)(?=[^<>]*>)|(\w)(?=\^>)(?=[^<>]*>)') #Returns the labels of upper toeholds.
lower_label = re.compile('(\w)(?=\^\*)(?=[^{}]*})') #Returns labels of lower toeholds/
open_bracket = re.compile('[\[{]*?(<)|[\[<]*?({)|[<{]*?(])') #Matches on open brackets [, { and <
close_bracket = re.compile('[\]}]*?(>)|[\]>]*?(})|[>}]*?(])') #Matches on close brackets ], } and >
empty_bracket = re.compile('(<(?:\s)*>)|({(?:\s)*})') #Matches on empty brackets like <>, {} and [ ].

upper_sequence = re.compile('((?<=\<).+?(?=\>))') #Pulls through characters between < and > chars, excluding the brackets themselves.
lower_sequence = re.compile('((?<=\{).+?(?=\}))') #Pulls through characters between { and } chars, excluding the brackets themselves.

spaces_start = re.compile('(?<=\<)(\s+?)(?=\w)|(?<=\{)(\s+?)(?=\w)|(?<=\[)(\s+?)(?=\w)') #Matches on spaces between brackets and words, like < A or { B
spaces_end = re.compile('(?<=\S)(\s)+?(?=\>)|(?<=\S)(\s)+?(?=\])|(?<=\S)(\s)+?(?=\})') #Matches on spaces before bracket closes


def find_sub_sequence(regex, seq):
    """Takes a regex and a sub sequence, and either returns the regex match or a blank string '' """
    if re.search(regex, seq) is not None:
        return re.search(regex, seq).group()
    return ""


def format_sequence(seq):
    """Remove unnecessary whitespaces and empty brackets"""
    seq = re.sub(spaces_end, '', re.sub(spaces_start, '', seq)) #Remove unneccesary whitespaces
    return re.sub(empty_bracket, '', seq)


class BindingRule(stocal.TransitionRule):
    """Join any two strings into their concatenations"""
    Transition = stocal.MassAction

    #TODO: Stop this binding rule from allowing any binding on to toeholds on overhangs.
    #Either don't match on them (how?) or regex check for toeholds on overhangs in if __name__ == '__main__':

    def novel_reactions(self, k, l):
        yield from self.toehold_binding(k, l, upper_label, lower_label)
        yield from self.toehold_binding(k, l, lower_label, upper_label)

    def toehold_binding(self, k, l, regex_1, regex_2):
        for match_1 in re.finditer(regex_1, k):
            for match_2 in re.finditer(regex_2, l):
                if match_1.group() == match_2.group():
                    if regex_1 == upper_label:
                        print("upper_th")
                        part_a = k[:match_1.start()] + re.search(close_bracket, k[match_1.start():]).group() + l[:match_2.start()] + re.search(close_bracket, l[match_2.start():]).group()
                        part_b = re.search(open_bracket, l[:match_2.end()]).group() + l[match_2.end() + 2:] + re.search(open_bracket, k[:match_1.end() + 1]).group() + k[match_1.end() + 1:]
                    else:
                        print("lower_th")
                        part_a = k[:match_1.start()] + re.search(close_bracket, k[match_1.start():]).group() + l[:match_2.start()] + re.search(close_bracket, l[match_2.start():]).group()
                        part_b = re.search(open_bracket, l[:match_2.end()]).group() + l[match_2.end() + 1:] + re.search(open_bracket, k[:match_1.end() + 1]).group() + k[match_1.end() + 2:]

                    draft_strand = part_a + "[" + match_2.group() + "^]" + part_b
                    final_strand = format_sequence(draft_strand)
                    print("final", final_strand)
                    yield self.Transition([k, l], [final_strand], alpha)


class UnbindingRule(stocal.TransitionRule):
    """Splits two strings when a toehold unbinds"""
    Transition = stocal.MassAction

    #TODO: Does this rule still work when toeholds are involved?
    def novel_reactions(self, kl):
        yield from self.toehold_unbinding(kl)

    def toehold_unbinding(self, kl):
        kl = format_sequence(kl)
        for double_th in re.finditer(double_toehold, kl):
            label = re.search(double_label, double_th.group()).group() #Retrieve the label of the toehold we are unbinding
            prefix, suffix = "", ""
            bracket_open = kl[:double_th.start()].rfind(']') #Possibly modify this to be min(.rfind(']'),.rfind(':') for overhangs
            bracket_close = kl[double_th.end():].find('[') #Possibly modify this to be max(.rfind('['),.rfind(':') for overhangs

            if bracket_open !=-1:
                prefix = kl[:bracket_open]
            else:
                bracket_open = 0
            if bracket_close != -1:
                suffix = kl[bracket_close:]
            else:
                bracket_close = len(double_th.group())

            upper_1 = find_sub_sequence(upper_sequence, kl[bracket_open:double_th.start()])
            lower_1 = find_sub_sequence(lower_sequence, kl[bracket_open:double_th.start()])
            upper_2 = find_sub_sequence(upper_sequence, kl[double_th.end() + 1:bracket_close])
            lower_2 = find_sub_sequence(lower_sequence, kl[double_th.end():bracket_close])

            print("prefix" ,prefix)
            print("upper 1", upper_1)
            print("label", label)
            print("upper 2", upper_2)
            print("lower 1", lower_1)
            print("lower 2", lower_2)
            print("suffix", suffix)
            #Todo: Fix this part. The upper and lower strands don't necessarily connect in this way.
            part_a = prefix + "<" + upper_1 + " " + label + "^ " + upper_2 + ">"
            part_a_fin = format_sequence(part_a)
            part_b = suffix + "{" + lower_1 + " " + label + "^*" + lower_2 + "}"
            part_b_fin = format_sequence(part_b)
            print("Final A:    ", part_a_fin)
            print("Final B:    ", part_b_fin)
            yield self.Transition([kl], [part_a_fin, part_b_fin], alpha)


process = stocal.Process(
    rules=[BindingRule(), UnbindingRule()]
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

#Unused code is below, in case of extensions being needed or previous regex being needed again.

# Extension (Allow toehold complements on the upper strand and toeholds on the lower strand)
# upper_th_c = re.compile('(\w)(?=\^\*)(?=[^<>]*>)')
# lower_th = re.compile('(\w)(?=\^\s)(?=[^{}]*})|(\w)(?=\^})(?=[^{}]*})')
#     yield from self.toehold_binding(k, l, lower_th, upper_th_c)
#     yield from self.toehold_binding(k, l, upper_th_c, lower_th)

#Original toehold_binding function
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

#bound_th = re.compile('(\w)(?=\^)(?=[^\[\]]*])')
#double_toehold = re.compile('\[\W*?(\w)\^\W*?\]')
#last_upper_th = re.compile('<([^>]*)>*(?!(?:.*<|.*\[))')
#last_upper_th = re.compile('<([^>]*)>')
#last_lower_th = re.compile('{([^}]*)}')
#last_lower_th = re.compile('({[^}]*)}*(?!(?:.*{|.*\[))')
#spaces_at_start = re.compile('<(\s)+?\w|{(\s)+?\w|\[(\s)+?\w')
#spaces_end = re.compile('(?:\w|\^|\*)(\s)+?>|(?:\w|\^|\*)(\s)+?\]|(?:\w|\^|\*)(\s)+?}') #Matches on spaces before bracket closes.


#Predecessor to find_sequence function. Can be deleted once testing shows the unbinding function works.
#if re.search(upper_sequence, kl[bracket_open:double_th.start()]) is not None:
#    upper_1 = re.search(upper_sequence, kl[bracket_open:double_th.start()]).group()
#    print("upper_1", upper_1)
#if re.search(lower_sequence, kl[bracket_open:double_th.start()]) is not None:
#    lower_1 = re.search(lower_sequence, kl[bracket_open:double_th.start()]).group()
            # if re.search(upper_sequence, kl[double_th.end() + 1:bracket_close]) is not None:
#     upper_2 = re.search(upper_sequence, kl[double_th.end() + 1:bracket_close]).group()
            # if re.search(lower_sequence, kl[double_th.end():bracket_close]) is not None:
#     lower_2 = re.search(lower_sequence, kl[double_th.end():bracket_close]).group()
