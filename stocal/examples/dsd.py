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
beta = 1000 ** -2

re_double = re.compile('(?:\[.*?\])')  # Matches on any double strand (includes brackets).
re_upper = re.compile('(?:\<.*?\>)')  # Matches on any upper strand (includes brackets).
re_lower = re.compile('(?:\{.*?\})')  # Matches on any lower strand (includes the brackets).

re_short_double_th = re.compile('(?:\[\W*?(\w)(?:\^\W*?\]))')  # Matches on double toeholds of the form [A^] not [A^ B]
re_double_lab = re.compile('(\w)(?=\^)(?=[^<>{}]*])')  # Returns the label of a double toehold regex.
re_upper_lab = re.compile('(\w)(?=\^)(?=[^<>]*>)') #Returns the labels of upper toeholds.
re_lower_lab = re.compile('(\w)(?=\^\*)(?=[^{}]*})')  # Returns labels of lower toeholds

re_opening = re.compile('[\[{]*?(<)|[\[<]*?({)|[<{]*?(])')  # Matches on open brackets [, { and <
re_closing = re.compile('[\]}]*?(>)|[\]>]*?(})|[>}]*?(])')  # Matches on close brackets ], } and >
re_empty = re.compile('(<(?:\s)*>)|({(?:\s)*})|\{\*\}')  # Matches on empty brackets like <>, {} and [ ].
re_spaces = re.compile(
    '(?<=\<)(\s+)|(?<=\{)(\s+)|(?<=\[)(\s+)|(\s)+(?=\>)|(\s)+(?=\])|(\s)+(?=\})')  # Matches on spaces at start or end of parts.

re_gate = re.compile(f"({re_lower.pattern})?({re_upper.pattern})?({re_double.pattern})({re_upper.pattern})?({re_lower.pattern})?")
re_strand = re.compile(fr"({re_lower.pattern})|({re_upper.pattern})")
re_cover_start = re.compile(r'(\w)(?=\^\*\s*\}\s*\<.*\1\^\s*\>)')
re_cover_end = re.compile(r'(?<=\<)\s*?(\w)(?=\^.*>\s*\{\s*\1\^\*)')


def find_sub_sequence(regex, seq):
    """Takes a regex and a sub sequence, and either returns the regex match (without the first and last chars) or a blank string '' """
    seq = re.search(regex, seq).group()
    if seq is not None:
        return seq[1:len(seq) - 1]
    return ""


def format_sequence(seq):
    """Remove unnecessary whitespaces and empty brackets"""
    seq = re.sub(re_spaces, '', seq)
    return re.sub(re_empty, '', seq)


class BindingRule(stocal.TransitionRule):
    """Join any two strings into their concatenations"""
    Transition = stocal.MassAction

    # TODO: Stop this binding rule from allowing any binding on to toeholds on overhangs.
    # Either don't match on them (how?) or regex check for toeholds on overhangs in if __name__ == '__main__':

    def novel_reactions(self, k, l):
        yield from self.toehold_binding(k, l, re_upper_lab, re_lower_lab)
        yield from self.toehold_binding(k, l, re_lower_lab, re_upper_lab)

    def toehold_binding(self, k, l, regex_1, regex_2):
        for match_1 in re.finditer(regex_1, k):
            for match_2 in re.finditer(regex_2, l):
                if match_1.group() == match_2.group():
                    if regex_1 == re_upper_lab:
                        print("upper_th")
                        part_a = k[:match_1.start()] + re.search(re_closing, k[match_1.start():]).group() + l[
                                                                                                            :match_2.start()] + re.search(
                            re_closing, l[match_2.start():]).group()
                        part_b = re.search(re_opening, l[:match_2.end()]).group() + l[match_2.end() + 2:] + re.search(
                            re_opening, k[:match_1.end() + 1]).group() + k[match_1.end() + 1:]
                    else:
                        print("lower_th")
                        part_a = k[:match_1.start()] + re.search(re_closing, k[match_1.start():]).group() + l[
                                                                                                            :match_2.start()] + re.search(
                            re_closing, l[match_2.start():]).group()
                        part_b = re.search(re_opening, l[:match_2.end()]).group() + l[match_2.end() + 1:] + re.search(
                            re_opening, k[:match_1.end() + 1]).group() + k[match_1.end() + 2:]
                    draft_strand = part_a + "[" + match_2.group() + "^]" + part_b
                    final_strand = format_sequence(draft_strand)
                    print("final", final_strand)
                    yield self.Transition([k, l], [final_strand], alpha)


class UnbindingRule(stocal.TransitionRule):
    """Splits two strings when a toehold unbinds"""
    Transition = stocal.MassAction

    # TODO: Does this rule still work when toeholds are involved?
    def novel_reactions(self, kl):
        print(kl)
        yield from self.toehold_unbinding(kl)

    def toehold_unbinding(self, kl):
        kl = format_sequence(kl)
        for double_th in re.finditer(re_short_double_th, kl):
            print("kl:", kl)
            label = re.search(re_double_lab,
                              double_th.group()).group()  # Retrieve the label of the toehold we are unbinding
            prefix, suffix = "", ""
            bracket_open = kl[:double_th.start()].rfind(
                ']')  # Possibly modify this to be min(.rfind(']'),.rfind(':') for overhangs
            bracket_close = kl[double_th.end():].find(
                '[')  # Possibly modify this to be max(.rfind('['),.rfind(':') for overhangs
            if bracket_open != -1:
                prefix = kl[:bracket_open + 1]
            else:
                bracket_open = 0
            if bracket_close != -1:
                suffix = kl[double_th.end() + bracket_close:]
            else:
                bracket_close = len(kl)

            # Find the upper strands before and after the double toehold. Likewise with lower strands.
            upper_1 = find_sub_sequence(re_upper, kl[bracket_open:double_th.start()])
            lower_1 = find_sub_sequence(re_lower, kl[bracket_open:double_th.start()])
            upper_2 = find_sub_sequence(re_upper, kl[double_th.end():double_th.end() + bracket_close])
            lower_2 = find_sub_sequence(re_lower, kl[double_th.end():double_th.end() + bracket_close])

            part_a = "<" + upper_1 + " " + label + "^ " + upper_2 + ">"
            part_b = "{" + lower_1 + " " + label + "^* " + lower_2 + "}"

            # Attach the prefix and/or suffix to the correct strand (upper or lower):
            if upper_1 == "":
                if upper_2 == "":
                    part_b = prefix + part_b + suffix
                else:
                    part_a = part_a + suffix
                    part_b = prefix + part_b
            else:
                if upper_2 == "":
                    part_a = prefix + part_a
                    part_b = part_b + suffix
                else:
                    part_a = prefix + part_a + suffix

            print("Final A:    ", format_sequence(part_a))
            print("Final B:    ", format_sequence(part_b))
            yield self.Transition([kl], [format_sequence(part_a), format_sequence(part_b)], alpha)


class CoveringRule(stocal.TransitionRule):
    """Splits two strings when a toehold unbinds"""
    Transition = stocal.MassAction

    # TODO: Does this rule still work when toeholds are involved?
    def novel_reactions(self, k):
        yield from self.toehold_covering(k)

    def toehold_covering(self, k):
        k = format_sequence(k)
        for gate in re.finditer(re_gate, k):
            for d_s in re.finditer(re_double, gate.group()): #DOES NOT NEED TO BE FINDITER AS THERE IS ONLY ONE DOUBLE STRAND
                cover_pos = re.search(re_cover_end, gate.group())
                if cover_pos is not None:
                    th_pos = gate.group()[cover_pos.end() + 1:].find(cover_pos.group() + "^*")
                    updated_gate = format_sequence(gate.group()[:d_s.end() - 1] + " " + cover_pos.group() + "]<" + gate.group()[cover_pos.end() + 1:cover_pos.end() + 1 + th_pos] + gate.group()[cover_pos.end() + th_pos + 4:])
                    print("updated_gate", updated_gate)
                    updated_seq = k[:gate.start()] + updated_gate + k[gate.end():]
                    print(updated_seq, "updated seq")
                    yield self.Transition([k], [updated_seq], alpha)
                cover_pos = re.search(re_cover_start, gate.group())
                if cover_pos is not None:
                    th_pos = gate.group()[cover_pos.start()+2:d_s.start()].find(cover_pos.group() + "^")
                    updated_gate = format_sequence(gate.group()[:cover_pos.start()] + gate.group()[cover_pos.end()+2:cover_pos.end()+th_pos]+">[" + cover_pos.group() + "^ " + gate.group()[d_s.start()+1:])
                    updated_seq = k[:gate.start()] + updated_gate + k[gate.end():]
                    print(updated_seq, "UPDATED SEQ")
                    yield self.Transition([k], [updated_seq], alpha)
                    #TODO: TEST ALL THIS AND MAKE IT WORK NICELY (CONDENSE IT)


                #print("gate", gate)
                #print(gate.start(), "GATES")
                updated_gate = self.toehold_covering_matching(gate.group())
        #else:
            #seq = self.toehold_covering_matching(k)
            #yield self.Transition([k], [seq], alpha)

    def toehold_covering_matching(self, gate):
        for d_s in re.finditer(re_double, gate):
            cover_pos = re.search(re_cover_end, gate)
            if cover_pos is not None:
                th_pos = gate[cover_pos.end() + 1:].find(cover_pos.group() + "^*")
                print(th_pos, "lower label")
                updated_gate = format_sequence(gate[:d_s.end() - 1] + " " + cover_pos.group() + "]<" + gate[cover_pos.end() + 1:cover_pos.end() + 1 + th_pos] + gate[cover_pos.end() + th_pos + 4:])
                print("updated_gate", updated_gate)
            return (cover_pos)



process = stocal.Process(
    rules=[CoveringRule()]
    # rules=[UnbindingRule()]
)

if __name__ == '__main__':
    # initial_state = {"{S' N^* L' R'}": 60, "<L N^ M N^>": 50}
    # initial_state = {"{ N^* L' R'}": 60, "<L N^ M N^>": 50}
    #initial_state = {"{N^* S' N^*}[C^]": 60, "<N^ M N^>": 50, "{L'}<L>[N^]<R>[M^]<S'>[A^]{B}" : 50}
    #initial_state = { "{L' N^ R'}": 60}

    # initial_state = {"<A>{B}[D^]<C^ F>{C^* G}": 60}
    initial_state = {"{A L^*}<B L^>[N^]<   R^>{B^*}:{L'}<L>[N^]<  F^>{B^*}": 500}
    # TODO:  re.fullmatch() does not work as I thought it did, so error checking needs to be updated so as to make sure every <, { and [ has a corresponding >, }, ].

    traj = process.sample(initial_state, tmax=100.)
    for _ in traj:
        print(traj.time, traj.state)

# Unused code is below, in case of extensions being needed or previous regex being needed again.

# Extension (Allow toehold complements on the upper strand and toeholds on the lower strand)

# re_upper_lab = re.compile('(\w)(?=\^\s)(?=[^<>]*>)|(\w)(?=\^>)(?=[^<>]*>)') #Returns the labels of upper toeholds.
# upper_th_c = re.compile('(\w)(?=\^\*)(?=[^<>]*>)')
# lower_th = re.compile('(\w)(?=\^\s)(?=[^{}]*})|(\w)(?=\^})(?=[^{}]*})')

#     yield from self.toehold_binding(k, l, lower_th, upper_th_c)
#     yield from self.toehold_binding(k, l, upper_th_c, lower_th)

#re_gate = re.compile('(.{2,}?(?=\:))|(\:)(?!.*\:)(.*)')

# Original toehold_binding function
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

# Unused regex. Kept temporarily for safekeeping.
# bound_th = re.compile('(\w)(?=\^)(?=[^\[\]]*])') #Matches on a bound double toehold label only.
# re_upper = re.compile('((?<=\<).+?(?=\>))')  # Pulls through characters between < and > chars, excluding the brackets themselves.
# re_lower = re.compile('((?<=\{).+?(?=\}))')  # Pulls through characters between { and } chars, excluding the brackets themselves.
# last_upper_th = re.compile('<([^>]*)>*(?!(?:.*<|.*\[))')
# last_upper_th = re.compile('<([^>]*)>')
# last_lower_th = re.compile('{([^}]*)}')
# last_lower_th = re.compile('({[^}]*)}*(?!(?:.*{|.*\[))')
# re_upper_lab_start = re.compile('((?<=\<)\s*?(\w)\^)')
# re_lower_lab_start = re.compile('((?<=\{)\s*?(\w)\^)')
# re_upper_lab_end = re.compile('(\w)\^(?=\s*\>)')
# re_lower_lab_end = re.compile('(\w)\^\*(?=\s*\})')
# re_spaces_end = re.compile('(?<=\S)(\s)+?(?=\>)|(?<=\S)(\s)+?(?=\])|(?<=\S)(\s)+?(?=\})')  # Matches on spaces before bracket closes
# re_spaces_start = re.compile('(?<=\<)(\s+?)(?=\w)|(?<=\{)(\s+?)(?=\w)|(?<=\[)(\s+?)(?=\w)')  # Matches on spaces between brackets and words, like < A or { B

# Predecessor to find_sequence function. Can be deleted once testing shows the unbinding function works.
# if re.search(upper_sequence, kl[bracket_open:double_th.start()]) is not None:
#    upper_1 = re.search(upper_sequence, kl[bracket_open:double_th.start()]).group()
#    print("upper_1", upper_1)
# if re.search(lower_sequence, kl[bracket_open:double_th.start()]) is not None:
#    lower_1 = re.search(lower_sequence, kl[bracket_open:double_th.start()]).group()
# if re.search(upper_sequence, kl[double_th.end() + 1:bracket_close]) is not None:
#     upper_2 = re.search(upper_sequence, kl[double_th.end() + 1:bracket_close]).group()
# if re.search(lower_sequence, kl[double_th.end():bracket_close]) is not None:
#     lower_2 = re.search(lower_sequence, kl[double_th.end():bracket_close]).group()

# Original failed checking:
# for key in initial_state:

# if re.fullmatch(format_check, key) is None:
#     format_correct = False
#     print("Incorrect Input:", key) #Some of the brackets must not match up. Raise a proper error in time.

# if format_correct:
# format_correct = True

# for upper_th in re.finditer(re_upper_lab_end, gate[:seq.start()]):
#     print("upper_th", upper_th)
#     for lower_thc in re.finditer(re_lower_lab_end, gate[:seq.start()]):
#         print("lower_th", lower_thc )
#         if upper_th.group() == lower_thc.group():
#             print(gate[:seq.start()], "gate[:seq.start()]")
#             print("PRE BIND", gate)


# for match_1 in re.finditer(regex_1, k):
#     for match_2 in re.finditer(regex_2, l):
#         if match_1.group() == match_2.group():
#             if regex_1 == upper_label:
#                 print("upper_th")
#                 part_a = k[:match_1.start()] + re.search(close_bracket, k[match_1.start():]).group() + l[                                                                                                            :match_2.start()] + re.search(
#                     close_bracket, l[match_2.start():]).group()
#                 part_b = re.search(open_bracket, l[:match_2.end()]).group() + l[match_2.end() + 2:] + re.search(
#                     open_bracket, k[:match_1.end() + 1]).group() + k[match_1.end() + 1:]
#             else:
#                 print("lower_th")
#                 part_a = k[:match_1.start()] + re.search(close_bracket, k[match_1.start():]).group() + l[:match_2.start()] + re.search(close_bracket, l[match_2.start():]).group()
#                 part_b = re.search(open_bracket, l[:match_2.end()]).group() + l[match_2.end() + 1:] + re.search(
#                     open_bracket, k[:match_1.end() + 1]).group() + k[match_1.end() + 2:]
#             draft_strand = part_a + "[" + match_2.group() + "^]" + part_b
#             final_strand = format_sequence(draft_strand)
#             print("final", final_strand)

# format = re.compile('<[^\[\]{}]*?>|\[[^<>{}]*?]|{[^<>\[\]]*?}|\s*:*')  # TODO: Fix this - speak to Harold.
# # Needs to check that each part of a strand meets one of these requirements. Fullmatch does not work.

# def format_sequence(seq):
#     #seq = re.sub(re_spaces_end, '', re.sub(re_spaces_start, '', seq))  # Remove unneccesary whitespaces
#     return re.sub(re_empty, '', seq)

            # x = re.search(r"(?<=\<)\s*?(\w)(?=\^.*>\s*\{\s*\1\^\*)",gate)
            # return("x")
            # for upper_seq in re.finditer(re_upper, gate[d_s.end():]):
            #     for lower_seq in re.finditer(re_lower, gate[d_s.end():]):
            #         for upper_th in re.finditer(re_upper_lab_start, upper_seq.group()):
            #             for lower_thc in re.finditer(re_lower_lab_start, lower_seq.group()):
            #                 if upper_th.group() == lower_thc.group():
            #                     updated_gate = format_sequence(gate[:d_s.end()-1] + " " + upper_th.group() + "]" + "<" + upper_seq.group()[upper_th.end():] + "{" + lower_seq.group()[lower_thc.end()+1:])
            #                     print(updated_gate)
            #                     return(updated_gate)

    # def toehold_covering_matching(self, gate):
    #     for d_s in re.finditer(re_double, gate):
    #         cover_pos = re.search(re_cover_end, gate)
    #         if cover_pos is not None:
    #             th_pos = gate[cover_pos.end() + 1:].find(cover_pos.group() + "^*")
    #             print(th_pos, "lower label")
    #             updated_gate = format_sequence(gate[:d_s.end() - 1] + " " + cover_pos.group() + "]<" + gate[cover_pos.end() + 1:cover_pos.end() + 1 + th_pos] + gate[cover_pos.end() + th_pos + 4:])
    #             print("updated_gate", updated_gate)
    #         return (cover_pos)
