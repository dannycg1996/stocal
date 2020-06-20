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

re_double = re.compile(r'(\[[^<{\[\]}>]*?\])')  # Matches on any double strand (includes brackets).
re_upper = re.compile(r'(<[^<\[{]*?>)')  # Matches on any upper strand (includes the brackets).
re_lower = re.compile(r'({[^<\[{]*?\})')  # Matches on any lower strand (includes the brackets).
re_short_double_th = re.compile(r'(?:\[\W*?(\w)(?:\^\W*?\]))')  # Matches on double toeholds of the form [A^] not [A^ B]
re_gate = re.compile(
    f"{re_lower.pattern}?{re_upper.pattern}?{re_double.pattern}{re_upper.pattern}?{re_lower.pattern}?")  # Matches on gates

re_double_lab = re.compile(r'(\w)(?=\^)(?=[^<>{}]*])')  # Returns the label of a double toehold regex.
re_upper_lab = re.compile(r'(\w)(?=\^)(?=[^<>]*>)')  # Returns the labels of upper toeholds.
re_lower_lab = re.compile(r'(\w)(?=\^\*)(?=[^{}]*})')  # Returns labels of lower toeholds

re_open = re.compile(r'<|\[|{')  # Matches on open brackets [, { and <
re_close = re.compile(r'(>|\]|})')  # Matches on close brackets ], } and >
re_empty = re.compile(r'(<(?:\s)*>)|({(?:\s)*})|(\[(?:\s)*])')  # Matches on empty brackets like <>, {} and [ ].
re_large_spaces = re.compile(r'(\s{2,})')  # Matches on spaces of length > 1
re_spaces = re.compile(r'(?<=[:>}\]<{\[])(\s+)|(\s+)(?=[:>\]}])')  # Matches on unnecessary spaces.

# The below 4 patterns match on different variants of gates which contain just a single upper or lower strand.
re_lone_upper_1 = re.compile(f"^{re_upper.pattern}::{re_gate.pattern}|(?<=::){re_upper.pattern}::{re_gate.pattern}")
re_lone_upper_2 = re.compile(f"{re_gate.pattern}::({re_upper.pattern})$")
re_lone_lower_1 = re.compile(f"^{re_lower.pattern}:{re_gate.pattern}|(?<=[^:]:){re_lower.pattern}:{re_gate.pattern}")
re_lone_lower_2 = re.compile(f"{re_gate.pattern}:{re_lower.pattern}$")

re_pre_cover = re.compile(r'(\w+)(?=\^\*\s*\}\s*\<.*(\1)\^\s*\>)')  # Matches where the Covering rule can be applied on a gate, before the d_s
re_post_cover = re.compile(r'(?<=\<)\s*?(\w+)(?=\^.*>\s*\{\s*(\1)\^\*)')  # Matches where the Covering rule can be applied on a gate, after the d_s

#    initial_state = {"{L'}<L>[S1]<S R2>:<L1>[S S2]<R>{R'}":60}
re_upper_migrate = re.compile(
    fr"{re_double.pattern}(<(\w+)[^<>:]*?>):{re_upper.pattern}(\[(\3)[^<>]*?\w\s*\])")   # Matches where upper strand migration can occur.
re_lower_migrate = re.compile(
    fr"{re_double.pattern}({{(\w+)[^<>:]*?\}})::{re_lower.pattern}(\[(\3)[^<>]*?\w\s*\])")  # Matches where lower strand migration can occur.
re_upper_migrate_r = re.compile(
    fr"(\[\w[^<>]*?\s(\w+)\s*\]){re_upper.pattern}:(<[^<>:]*?(\2)\s*>){re_double.pattern}")  # Matches where upper strand rev migration can occur.
re_lower_migrate_r = re.compile(
    fr"(\[\w[^<>]*?\s(\w+)\s*\]){re_lower.pattern}::({{[^<>:]*?(\2)\s*}}){re_double.pattern}") # Matches where lower strand rev migration can occur.

re_reduce_upper = re.compile(
    fr"{re_double.pattern}<(\w+)([^<>:]*?)>:{re_upper.pattern}\[(\2)\]{re_upper.pattern}?{re_lower.pattern}?")
re_reduce_lower = re.compile(
    fr"{re_double.pattern}{{(\w+)([^{{}}:]*?)}}::{re_lower.pattern}\[(\2)\]{re_upper.pattern}?{re_lower.pattern}?")
re_reduce_upper_r = re.compile(
    fr"{re_lower.pattern}?{re_upper.pattern}?\[(\w)\]{re_upper.pattern}?:<([^<>:]*?)(\3)>{re_double.pattern}"
)
re_reduce_lower_r = re.compile(
    fr"{re_lower.pattern}?{re_upper.pattern}?\[(\w)\]{re_lower.pattern}?::{{([^<>:]*?)(\3)}}{re_double.pattern}"
)
#"{L'}<L>[S]<R S>:<L2>[S1]<L>{L'}" : 60
#TODO: Finish the reduce rule.

re_format_1 = re.compile(
    f"({re_double.pattern}{re_upper.pattern}{re_lower.pattern}?::{re_lower.pattern}?{re_upper.pattern}{re_double.pattern})")
re_format_2 = re.compile(
    f"({re_double.pattern}{re_upper.pattern}?{re_lower.pattern}:{re_lower.pattern}{re_upper.pattern}?{re_double.pattern})")
re_format_3 = re.compile(
    f"({re_double.pattern}{re_upper.pattern}{re_lower.pattern}?::{re_lower.pattern}?{re_double.pattern})")
re_format_4 = re.compile(
    f"({re_double.pattern}{re_upper.pattern}?{re_lower.pattern}:{re_upper.pattern}?{re_double.pattern})")


def find_sub_sequence(regex, seq):
    """Takes a regex and a sub sequence, and either returns the regex match (without the first and last chars) or a blank string '' """
    seq = re.search(regex, seq)
    if seq is not None:
        return seq.group()[1:len(seq.group()) - 1]
    return ""


def find_sub_seq(seq, trunc=True):
    """Takes a sub sequence, and either returns the regex match (without the first and last chars) or a blank string '' """
    if seq is not None:
        if trunc is True:
            return seq[1:len(seq) - 1]
        else:
            return seq
    return ""


def tidy(sys):
    """Remove unnecessary whitespaces and empty brackets"""
    sys = re.sub(re_large_spaces, " ", sys)  # Replaces spaces of length 2 or more with single spaces.
    sys = re.sub(re_spaces, '', sys)  # Remove unnecessary spaces
    sys = re.sub(re_empty, '', sys)  # Remove empty brackets
    return sys


def fix_upper_gate(sys, match_obj, i):
    """This function takes a system sys, a match object and a starting index. The match object identifies gates which consist solely of
     an upper strand, merges it with a gate to the right, and then returns the updated system"""
    if match_obj.group(3+i) is not None:  # Match object has 6 groups: (< >)::({ })(< >)([ ])(< >)({ })
        print("A")
        strand = sys[:match_obj.start()] + sys[match_obj.end(1+i)+2:match_obj.start(3+i)+1] +\
            match_obj.group(1+i)[1:len(match_obj.group(1+i))-1] + " " + sys[match_obj.start(3+i)+1:]
        print("Strand", strand)
    elif match_obj.group(2+i) is not None:
        print("B")
        strand = sys[:match_obj.start()] + match_obj.group(2+i) + match_obj.group(1+i) + sys[match_obj.start(4+i):]
    else:
        print("C")
        strand = sys[:match_obj.start()] + match_obj.group(1+i) + sys[match_obj.start(4+i):]
    print(strand, "strand 2")
    return strand


def fix_lower_gate(sys, match_obj, i):
    """This function takes a system sys, a match object and a starting index. The match object identifies gates which consist solely of
     a lower strand, merges it with a gate to the right, and then returns the updated system"""
    if match_obj.group(2+i) is not None: # Match object has 6 groups: ({ })::({ })(< >)([ ])(< >)({ })
        strand = sys[:match_obj.start()] + sys[match_obj.end(1+i)+1:match_obj.start(2+i)+1] +\
            match_obj.group(1+i)[1:len(match_obj.group(1+i))-1] + " " + sys[match_obj.start(2+i)+1:]
    else:
        strand = sys[:match_obj.start()] + match_obj.group(1+i) + sys[match_obj.end(1+i)+1:]
    return strand


def merge_gates(sys):
    """This function identifies gates which only contain a single upper or lower strand, and merges this strand to an adjacent gate, with
    the following gate taking priority over the previous gate"""
    upper_g_1 = re.search(re_lone_upper_1, sys)  # Matches on ^< >::{gate} or ::< >::{gate}
    upper_g_2 = re.search(re_lone_upper_2, sys)  # Matches on {gate}::< >$
    lower_g_1 = re.search(re_lone_lower_1, sys)  # Matches on ^{ }:{gate} or :{ }:{gate}
    lower_g_2 = re.search(re_lone_lower_2, sys)  # Matches on {gate}:{ }$

    if upper_g_1 is not None:
        if upper_g_1.group(4) is not None:  # If 1st match condition of upper_g_1 is met.
            return merge_gates(fix_upper_gate(sys, upper_g_1, 0))
        else: # If 2nd match condition of upper_g_1 is met.
            return merge_gates(fix_upper_gate(sys, upper_g_1, 6))
    elif upper_g_2 is not None:
        if upper_g_2.group(4) is not None:  # If gate before the upper strand gate had an upper strand after the double strand
            strand = sys[:upper_g_2.end(4)-1] + " " + upper_g_2.group(6)[1:] + sys[upper_g_2.end(4):upper_g_2.start(6)-2]
        else:
            strand = sys[:upper_g_2.end(3)] + upper_g_2.group(6) + sys[upper_g_2.end(3):upper_g_2.start(6)-2]
        return merge_gates(strand)
    elif lower_g_1 is not None:  # If 1st match condition of lower_g_1 is met.
        if lower_g_1.group(4) is not None:
            return merge_gates(fix_lower_gate(sys, lower_g_1, 0))
        else:  # If 2nd match condition of lower_g_1 is met.
            return merge_gates(fix_lower_gate(sys, lower_g_1, 6))
    elif lower_g_2 is not None:
        if lower_g_2.group(5) is not None: # If gate before the lower strand gate had a lower strand after the double strand
            strand = sys[:lower_g_2.end(5)-1] + " " + lower_g_2.group(6)[1:]
        else:
            strand = sys[:lower_g_2.start(6)-1] + lower_g_2.group(6)
        return merge_gates(strand)
    else:
        return sys


def reformat(sys):
    """This function identifies non-standard patterns and re-formats it. For example, {A}<B>[C]<D>{E}::{F}<G>[H] must be rewritten as
    {A}<B>[C]{E}::{F}<D G>[H] to ensure that the reaction is reversible and the results are clear"""
    format_1 = re.search(re_format_1, sys)
    format_2 = re.search(re_format_2, sys)
    format_3 = re.search(re_format_3, sys)
    format_4 = re.search(re_format_4, sys)

    if format_1 is not None:
        upper = format_1.group(3)[1:len(format_1.group(3)) - 1] + " "
        return reformat(sys[:format_1.start(3)] + sys[format_1.end(3):format_1.start(6) + 1] + upper + sys[format_1.start(6) + 1:])
    elif format_2 is not None:
        lower = format_2.group(4)[1:len(format_2.group(4)) - 1] + " "
        new = sys[:format_2.start(4)] + sys[format_2.end(4):format_2.start(5) + 1] + lower + sys[format_2.start(5) + 1:]
        return reformat(new)
    elif format_3 is not None:
        new = sys[:format_3.start(3)] + sys[format_3.end(3):format_3.start(6)] + format_3.group(3) + sys[format_3.start(6):]
        return reformat(new)
    elif format_4 is not None:
        new = sys[:format_4.start(4)] + ":" + format_4.group(4) + sys[format_4.end(4) + 1:]
        return reformat(new)
    else:
        return sys


def standardise(sys):
    """This function calls three other functions, which act to standardise a system. The format_seq function removes unnecessary spaces and empty
    brackets, the merge_lone_strands function appropriately merges gates which contain only a single strand, and the final function standardises
    sequences which can be described in several ways with Lakin's syntax."""
    sys = tidy(sys)
    sys = merge_gates(sys)
    sys = reformat(sys)
    return sys


class BindingRule(stocal.TransitionRule):
    """Join any two strings into their concatenations"""
    Transition = stocal.MassAction

    def novel_reactions(self, k, l):
        if re.search(re_gate, k) is None or re.search(re_gate, l) is None:
            if re.search(re_gate, k) is not None or re.search(re_gate, l) is not None:
                # TODO: Can I avoid calling this function 4 times? Maybe 2 times instead?
                yield from self.strand_to_gate_binding(k, l, re_upper_lab, re_lower_lab)
                yield from self.strand_to_gate_binding(l, k, re_upper_lab, re_lower_lab)
                yield from self.strand_to_gate_binding(k, l, re_lower_lab, re_upper_lab)
                yield from self.strand_to_gate_binding(l, k, re_lower_lab, re_upper_lab)
            else:
                yield from self.strand_to_strand_binding(k, l, re_upper_lab, re_lower_lab)
                yield from self.strand_to_strand_binding(k, l, re_lower_lab, re_upper_lab)

    def strand_to_gate_binding(self, k, l, regex_1, regex_2):
        if re.search(regex_1, l) is not None:
            for gate in re.finditer(re_gate, k):
                for match in re.finditer(regex_2, gate.group()):
                    for match_2 in re.finditer(regex_1, l):
                        if match.group() == match_2.group():
                            d_s = "[" + match.group() + "^]"
                            i = gate.start()
                            if regex_1 == re_lower_lab:
                                seq_start = "{" + l[1:match_2.start()] + "}"
                                seq_end = "{" + l[match_2.end() + 2:len(l) - 1] + "}"
                                if match.start() > gate.start(2) - i and match.end() < gate.end(2) - i:
                                    u_s_1 = "<" + k[gate.start(2) + 1:match.start() + i] + ">"
                                    u_s_2 = "<" + k[match.end() + 1 + i:gate.end(2) - 1] + ">"
                                    seq = k[:gate.start()] + seq_start + u_s_1 + d_s + seq_end + "::" + gate.group(1) + u_s_2 + k[gate.start(3):]
                                    yield self.Transition([k, l], [standardise(seq)], alpha)
                                elif match.start() > gate.start(4) - i and match.end() < gate.end(4) - i:
                                    u_s_1 = "<" + k[gate.start(4) + 1:match.start() + i] + ">"
                                    u_s_2 = "<" + k[match.end() + i + 1:gate.end(4) - 1] + ">"
                                    seq = k[:gate.end(3)] + gate.group(5) + "::" + seq_start + u_s_1 + d_s + u_s_2 + seq_end + k[gate.end():]
                                    yield self.Transition([k, l], [standardise(seq)], alpha)
                            else:
                                seq_start = "<" + l[1:match_2.start()] + ">"
                                seq_end = "<" + l[match_2.end() + 2:len(l) - 1] + ">"
                                if match.start() > gate.start(1) - i and match.end() < gate.end(1) - i:
                                    l_s_1 = "{" + k[gate.start(1) + 1:match.start() + i] + "}"
                                    l_s_2 = "{" + k[match.end() + i + 2:gate.end(1) - 1] + "}"
                                    seq = k[:gate.start(1)] + l_s_1 + seq_start + d_s + seq_end + l_s_2 + ":" + k[gate.start(2):]
                                    yield self.Transition([k, l], [standardise(seq)], alpha)
                                elif match.start() > gate.start(5) - i and match.end() < gate.end(5) - i:
                                    l_s_1 = "{" + k[gate.start(5) + 1:match.start() + i] + "}"
                                    l_s_2 = "{" + k[match.end() + i + 2:gate.end(5) - 1] + "}"
                                    seq = k[:gate.end(4)] + ":" + l_s_1 + seq_start + d_s + seq_end + l_s_2 + k[gate.end():] #("SEQ SECOND LOWER", seq)
                                    yield self.Transition([k, l], [standardise(seq)], alpha)

    def strand_to_strand_binding(self, k, l, regex_1, regex_2):
        #TODO: Tidy this up; the two part As are the same, just with the halves in different orders. Likewise for Part B.
        for match_1 in re.finditer(regex_1, k):
            for match_2 in re.finditer(regex_2, l):
                if match_1.group() == match_2.group():
                    d_s = "[" + match_2.group() + "^]"
                    if regex_1 == re_upper_lab:
                        part_a = l[:match_2.start()] + re.search(re_close, l[match_2.start():]).group() + \
                             k[:match_1.start()] + re.search(re_close, k[match_1.start():]).group()
                        part_b = re.search(re_open, k[:match_1.end() + 1]).group() + k[match_1.end() + 1:] + \
                                 re.search(re_open, l[:match_2.end()]).group() + l[match_2.end() + 2:]
                    else:
                        part_a = k[:match_1.start()] + re.search(re_close, k[match_1.start():]).group() +\
                                 l[:match_2.start()] + re.search(re_close, l[match_2.start():]).group()
                        part_b = re.search(re_open, l[:match_2.end()]).group() + l[match_2.end() + 1:] + \
                                 re.search(re_open, k[:match_1.end() + 1]).group() + k[match_1.end() + 2:]
                    # print("final", format_seq(part_a + d_s + part_b))
                    yield self.Transition([k, l], [tidy(part_a + d_s + part_b)], alpha)


class UnbindingRule(stocal.TransitionRule):
    """Splits two systems when a toehold unbinds"""
    Transition = stocal.MassAction

    def novel_reactions(self, kl):
        yield from self.toehold_unbinding(kl)

    def toehold_unbinding(self, kl):
        kl = tidy(kl)
        for gate in re.finditer(re_gate, kl):
            d_s = re.search(re_short_double_th, gate.group())
            if d_s is not None:
                label = re.search(re_double_lab, d_s.group()).group()
                upper_1 = find_sub_seq(gate.group(2))
                lower_1 = find_sub_seq(gate.group(1))
                upper_2 = find_sub_seq(gate.group(4))
                lower_2 = find_sub_seq(gate.group(5))
                #upper_1 = find_sub_sequence(re_upper, gate.group()[:d_s.start()])
                #lower_1 = find_sub_sequence(re_lower, gate.group()[:d_s.start()])
                # upper_2 = find_sub_sequence(re_upper, gate.group()[d_s.end():])
                # lower_2 = find_sub_sequence(re_lower, gate.group()[d_s.end():])
                part_a = "<" + upper_1 + " " + label + "^ " + upper_2 + ">"
                part_b = "{" + lower_1 + " " + label + "^* " + lower_2 + "}"
                if gate.start() > 0:
                    if kl[gate.start() - 2:gate.start()] == "::":
                        part_a = kl[:gate.start()] + part_a
                    else:
                        part_b = kl[:gate.start()] + part_b
                if gate.end() < len(kl):
                    if kl[gate.end():gate.end() + 2] == "::":
                        part_a = part_a + kl[gate.end():]
                    else:
                        part_b = part_b + kl[gate.end():]

                print("FINAL A:", standardise(part_a), "FINAL B:", standardise(part_b))
                yield self.Transition([kl], [standardise(part_a), standardise(part_b)], alpha)


class CoveringRule(stocal.TransitionRule):
    """This rule carries out the toehold covering reaction, where an exposed toehold in a lower strand is covered by a complementary
     exposed toehold in the upper strand"""
    Transition = stocal.MassAction

    def novel_reactions(self, k):
        yield from self.toehold_covering(k)

    def toehold_covering(self, k):
        k = tidy(k)
        print("No")
        for gate in re.finditer(re_gate, k):
            pre_cover = re.search(re_pre_cover, gate.group())
            post_cover = re.search(re_post_cover, gate.group())
            if pre_cover is not None:
                updated_gate = gate.group()[:pre_cover.start()] + gate.group()[pre_cover.end() + 2: pre_cover.start(2)] + ">[" + \
                    pre_cover.group() + "^ " + gate.group()[gate.start(3)+1:]
                updated_sys = k[:gate.start()] + tidy(updated_gate) + k[gate.end():]
                print(updated_sys, "updated seq")
                yield self.Transition([k], [updated_sys], alpha)
            if post_cover is not None:
                updated_gate = gate.group()[:gate.end(3) - 1] + " " + post_cover.group() + "^]<" + \
                               gate.group()[post_cover.end()+1:gate.end(4)] + "{" + gate.group()[post_cover.end(2)+2:]
                updated_sys = k[:gate.start()] + tidy(updated_gate) + k[gate.end():]
                print(updated_sys, "updated_sys seq")
                yield self.Transition([k], [updated_sys], alpha)


class MigrationRule(stocal.TransitionRule):
    """Migrates an upper or lower overhang up/down a strand via branch migration"""
    Transition = stocal.MassAction

    def novel_reactions(self, k):
        #  TODO: Does the overhang have to be a double overhang (i.e. what if L1 is missing in the example?)
        k = tidy(k)
        yield from self.migrate(k, re_lower_migrate, re_lower)
        yield from self.migrate(k, re_upper_migrate, re_upper)
        yield from self.migrate_rev(k, re_lower_migrate_r, re_lower)
        yield from self.migrate_rev(k, re_upper_migrate_r, re_upper)

    def migrate(self, k, regex_1, regex_2):
        for match in re.finditer(regex_1, k):
            i = match.start()
            d_s_1 = match.group(1)[:len(match.group(1))-1] + " " + match.group(3) + "]"
            d_s_2 = "[" + match.group()[match.end(6)-i:match.end(5)-i]
            if regex_2 == re_lower:
                strand_1 = "{" + match.group()[match.end(3)-i:match.end(2)-i]
                strand_2 = match.group(4)[:len(match.group(4))-1] + " " + match.group(3) + "}"
                bracket = "::"
            else:
                strand_1 = "<" + match.group()[match.end(3)-i:match.end(2)-i]
                strand_2 = match.group(4)[:len(match.group(4))-1] + " " + match.group(3) + ">"
                bracket = ":"
            seq = tidy(k[:match.start()] + d_s_1 + strand_1 + bracket + strand_2 + d_s_2 + k[match.end():])
            print("seq", seq)
            yield self.Transition([k], [seq], alpha)

    def migrate_rev(self, k, regex_1, regex_2):
        for match in re.finditer(regex_1, k):
            i = match.start()
            d_s_1 = match.group()[:match.start(2)-i] + "]"
            d_s_2 = "[" + match.group(2) + " " + match.group(6)[1:]
            if regex_2 == re_lower:
                strand_1 = "{" + match.group(2) + " " + match.group(3)[1:]
                strand_2 = match.group()[match.start(4)-i: match.start(5)-i] + "}"
                bracket = "::"
            else:
                strand_1 = "<" + match.group(2) + " " + match.group(3)[1:]
                strand_2 = match.group()[match.start(4)-i: match.start(5)-i] + ">"
                bracket = ":"
            seq = tidy(k[:match.start()] + d_s_1 + strand_1 + bracket + strand_2 + d_s_2 + k[match.end():])
            print("seq", seq)
            yield self.Transition([k], [seq], alpha)


class ReductionRule(stocal.TransitionRule):
    """Splits two strings when a toehold unbinds"""
    Transition = stocal.MassAction

    def novel_reactions(self, k):
        k = tidy(k)
        yield from self.reduction_fwd(k, re_reduce_upper)
        yield from self.reduction_fwd(k, re_reduce_lower)
        yield from self.reduction_rev(k, re_reduce_upper_r)
        yield from self.reduction_rev(k, re_reduce_lower_r)

    def reduction_fwd(self, k, regex_1):
        print("K reduce:", k)
        for match in re.finditer(regex_1, k):
            strand_1 = find_sub_seq(match.group(4)) + " " + match.group(2) + " "
            start = k[:match.end(1)-1] + " " + match.group(2) + "]"
            if regex_1 == re_reduce_upper:
                strand_1 = tidy("<" + strand_1 + find_sub_seq(match.group(6)) + ">")
                strand_2 = tidy(start + "<" + find_sub_seq(match.group(3), False) + ">" + find_sub_seq(match.group(7), False) + k[match.end():])
            else:
                strand_1 = tidy("{" + strand_1 + find_sub_seq(match.group(7)) + "}")
                strand_2 = tidy(start + " " + find_sub_seq(match.group(6), False) + "{" + find_sub_seq(match.group(3), False) + "}" + k[match.end():])
            print("strand_1", strand_1, "strand_2", strand_2)

            yield self.Transition([k], [strand_1 , strand_2], alpha)

    def reduction_rev(self, k, regex_1):
        for match in re.finditer(regex_1, k):
            if regex_1 == re_reduce_upper_r:
                strand_1 = "<" + find_sub_seq(match.group(2)) + " " + match.group(3) + " " + find_sub_seq(match.group(4)) + ">"
                strand_2 = find_sub_seq(match.group(1), False) + "<" + find_sub_seq(match.group(5), False) + ">[" + match.group(3) + " " + match.group(7)[1:]
            else:
                strand_1 = "{" + find_sub_seq(match.group(1))  + " " + match.group(3) + " " + find_sub_seq(match.group(4)) + "}"
                strand_2 = "{" + find_sub_seq(match.group(5), False) + "}" + find_sub_seq(match.group(2), False) + "[" + match.group(3) + " " + match.group(7)[1:]
            strand_2 = k[:match.start()] + strand_2 + k[match.end():]
            print("A:", tidy(strand_1), "B:", tidy(strand_2))
            yield self.Transition([k], [tidy(strand_1), tidy(strand_2)], alpha)



process = stocal.Process(
    rules=[BindingRule()]
)

if __name__ == '__main__':
    # initial_state = {"<A>{B}[D^]<C^ F>{C^* G}": 60}
    # initial_state = {"<Z Y C>[B]::<E F G>::[K]": 60}
    initial_state = {"{A}<B>[C^]<D>{E}::{F}<G>[H^]<I>{J}::{K}<L>[M^]<N>{O}": 500}  # THIS ONE
    initial_state = {"<L1 N^ S R1>": 60, "{L' N^*}<L>[S R2]<R>{R'}": 60}
    initial_state = {"{L'}<L1>[N^]<S R1>:<L>[S R2]<R>{R'}" : 60}
    initial_state = {"{N^*}<R N^>[S]<A B C>{D E}" : 60}
    #initial_state = {"{L'}<L>[S1]<S R2>:<L1>[S S2]<R>{R'}":60}
    #initial_state = {"{L'}<L>[S1]{S R2}::{L1}[S S2]<R>{R'}":60}
    initial_state = {"{L'}<L>[S1 S]<R2>:<L1 S>[S2]<R>{R'}": 60}
    initial_state = {"{L'}<L>[S1]<S R>:<L2>[S]<R2>{R'}" : 60}
    initial_state = {"{L'}<L>[S]<L2>:<R S>[S1]<R2>{R'}" : 60}
    initial_state = {"{L'}<L>[S]{L2}::{R S}[S1]<R2>{R'}" : 60}
    #initial_state = {"{L'}<L>[S1]{S R}::{L2}[S]<R2>{R'}" : 60}
    #initial_state = {"{A}<B>[C^]<D>:{E F^* G}<H>[I]<J>{K}": 60, "<Z F^ X>": 60}
    # initial_state = {"{A}<B>[C^]{E}::{K}<D G H^ I L>[M^]<O>{Z N^* G}" :600, "<F N^ J>":600}
    # initial_state = {"{L'}<L>[S1]<S R2 R3>:<L1>[S R2 S2]<R>{R'}" : 6000000}
    # initial_state = {"{L'}<L>[S1]<S R>:<L2>[S]<R2>{R'}" : 60}
    # initial_state = {"{L'}<L>[S]<N^ R>{N^* R'}" : 60}

    initial_state = {standardise(key): value for key, value in initial_state.items()}
    #print("init 2", initial_state)
    # re_lone_upper_1 = re.compile(f"^({re_upper.pattern})::|(?<=::)({re_upper.pattern})::")
    #x = "<A B C>::<F G>[H^ I]{L}::<A B C>::<H>[Y^]<N>{Z}::<A B C>"
    x = "{F}:<A B C>[D^]<M>{J}::<A B>::{F}<F>[G^]"
    #print("X", x)
    #z = standardise(x)
    #print("STANDARDISE", z)
    traj = process.sample(initial_state, tmax=1000000000.)
    for _ in traj:
        print(traj.time, traj.state)

# Unused code is below, in case of extensions being needed or previous regex being needed again.

# Extension (Allow toehold complements on the upper strand and toeholds on the lower strand)

# upper_th_c = re.compile('(\w)(?=\^\*)(?=[^<>]*>)')
# lower_th = re.compile('(\w)(?=\^\s)(?=[^{}]*})|(\w)(?=\^})(?=[^{}]*})')

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
# re_upper = re.compile('((?<=\<).+?(?=\>))')  # Pulls through characters between < and > chars, excluding the brackets themselves.
# re_lower = re.compile('((?<=\{).+?(?=\}))')  # Pulls through characters between { and } chars, excluding the brackets themselves.

# if re.fullmatch(format_check, key) is None:
#     format_correct = False
#     print("Incorrect Input:", key) #Some of the brackets must not match up. Raise a proper error in time.


# for upper_th in re.finditer(re_upper_lab_end, gate[:seq.start()]):
#     print("upper_th", upper_th)
#     for lower_thc in re.finditer(re_lower_lab_end, gate[:seq.start()]):
#         print("lower_th", lower_thc )
#         if upper_th.group() == lower_thc.group():
#             print(gate[:seq.start()], "gate[:seq.start()]")
#             print("PRE BIND", gate)
# re_lower = re.compile(r'(?:\{.*?\})')  # Matches on any lower strand (includes the brackets).

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

# def toehold_unbinding(self, kl):
#     kl = format_sequence(kl)
#     for double_th in re.finditer(re_short_double_th, kl):
#         print("kl:", kl)
#         label = re.search(re_double_lab,
#                           double_th.group()).group()  # Retrieve the label of the toehold we are unbinding
#         prefix, suffix = "", ""
#         bracket_open = kl[:double_th.start()].rfind(']')  # Possibly modify this to be min(.rfind(']'),.rfind(':') for overhangs
#         bracket_close = kl[double_th.end():].find('[')  # Possibly modify this to be max(.rfind('['),.rfind(':') for overhangs
#         if bracket_open != -1:
#             prefix = kl[:bracket_open + 1]
#         else:
#             bracket_open = 0
#         if bracket_close != -1:
#             suffix = kl[double_th.end() + bracket_close:]
#         else:
#             bracket_close = len(kl)
#
#         # Find the upper strands before and after the double toehold. Likewise with lower strands.
#         upper_1 = find_sub_sequence(re_upper, kl[bracket_open:double_th.start()])
#         lower_1 = find_sub_sequence(re_lower, kl[bracket_open:double_th.start()])
#         upper_2 = find_sub_sequence(re_upper, kl[double_th.end():double_th.end() + bracket_close])
#         lower_2 = find_sub_sequence(re_lower, kl[double_th.end():double_th.end() + bracket_close])
#
#         part_a = "<" + upper_1 + " " + label + "^ " + upper_2 + ">"
#         part_b = "{" + lower_1 + " " + label + "^* " + lower_2 + "}"
#
#         # Attach the prefix and/or suffix to the correct strand (upper or lower):
#         if upper_1 == "":
#             if upper_2 == "":
#                 part_b = prefix + part_b + suffix
#             else:
#                 part_a = part_a + suffix
#                 part_b = prefix + part_b
#         else:
#             if upper_2 == "":
#                 part_a = prefix + part_a
#                 part_b = part_b + suffix
#             else:
#                 part_a = prefix + part_a + suffix
#
#         print("Final A:    ", format_sequence(part_a))
#         print("Final B:    ", format_sequence(part_b))
#         yield self.Transition([kl], [format_sequence(part_a), format_sequence(part_b)], alpha)

# re_opening = re.compile(r'[\[{]*?(<)|[\[<]*?({)|[<{]*?(])')  # Matches on open brackets [, { and <
# re_closing = re.compile(r'[\]}]*?(>)|[\]>]*?(})|[>}]*?(])')  # Matches on close brackets ], } and >

# def toehold_binding(self, k, l, regex_1, regex_2):
#     print("k",k, "l", l)
#     if re.search(re_lower_lab, l) is not None:
#         for gate in re.finditer(re_gate, k):
#             print("gate", gate.group())
#             for match in re.finditer(re_upper_lab, gate.group()):
#                 print("matching 1", match.group(), match.start(), match.end())
#                 for match_2 in re.finditer(re_lower_lab, l):
#                     print("matching 2", match_2.group(),match_2.start())
#                     if match.group() == match_2.group():
#                         i = gate.start()
#                         seq_start = l[1:match_2.start()]
#                         seq_end = l[match_2.end()+2:len(l)-1]
#                         d_s = "[" + match.group() + "^]"
#                         if match.start() > gate.start(2)-i and match.end() < gate.end(2)-i:
#                             print("First upper")
#                             u_s_1 = "<" + k[gate.start(2)+1:match.start()+i] + ">"
#                             u_s_2 = "<" + k[match.end()+1+i:gate.end(2)-1] + ">"
#                             seq = k[:gate.start()] + "{" + seq_start + "}" + u_s_1 + d_s + "{" + seq_end + "}::" + gate.group(1) + u_s_2 + k[gate.start(3):]
#                             yield self.Transition([k, l], [format_seq(seq)], alpha)
#                         elif match.start() > gate.start(4)-i and match.end() < gate.end(4)-i:
#                             print("second upper")
#                             u_s_1 = "<" + k[gate.start(4)+1:match.start()+i] + ">"
#                             u_s_2 = "<" + k[match.end()+i+1:gate.end(4)-1] + ">"
#                             seq = k[:gate.end(3)] + gate.group(5) + "::" + "{" + seq_start + "}" + u_s_1 + d_s + u_s_2  + "{" + seq_end + "}" + k[gate.end():]
#                             print(k, l, seq)
#                             yield self.Transition([k, l], [format_seq(seq)], alpha)
#     elif re.search(re_upper_lab, l) is not None:
#        for gate in re.finditer(re_gate, k):
#             print("gate", gate.group())
#             for match in re.finditer(re_lower_lab, gate.group()):
#                 print("matching 1", match.group(), match.start(), match.end())
#                 for match_2 in re.finditer(re_upper_lab, l):
#                     print("matching 2", match_2.group(),match_2.start())
#                     if match.group() == match_2.group():
#                         i = gate.start()
#                         seq_start = l[1:match_2.start()]
#                         seq_end = l[match_2.end()+2:len(l)-1]
#                         d_s = "[" + match.group() + "^]"
#                         if match.start() > gate.start(1)-i and match.end() < gate.end(1)-i:
#                             l_s_1 = "{" + k[gate.start(1)+1:match.start()+i] + "}"
#                             l_s_2 = "{" + k[match.end()+i+2:gate.end(1)-1] + "}"
#                             print("First lower")
#                             seq = k[:gate.start(1)] + l_s_1 + "<" + seq_start + ">" + d_s + "<" + seq_end + ">" + l_s_2 + ":" + k[gate.start(2):]
#                             print("SEQ FIRST LOWER", standardise(seq))
#                             yield self.Transition([k, l], [standardise(seq)], alpha)
#                         elif match.start() > gate.start(5)-i and match.end() < gate.end(5)-i:
#                             print("Second lower")
#                             l_s_1 = "{" + k[gate.start(5)+1:match.start()+i] + "}"
#                             l_s_2 = "{" + k[match.end()+i+2:gate.end(5)-1] + "}"
#                             seq = k[:gate.end(4)] + ":" + l_s_1 + "<" + seq_start + ">" + d_s + "<" + seq_end + ">" + l_s_2 + k[gate.end():]
#                             print("SEQ SECOND LOWER", seq)
#                             yield self.Transition([k, l], [format_seq(seq)], alpha)

# def standardise(seq):
#     """Identifies gates which only contain a single upper or lower strand, and adds this strand to an adjacent gate, with
#     the following gate taking priority over the previous gate"""
#     seq = format_seq(seq)
#
#     upper_g_1 = re.search(re_lone_upper_1, seq)
#     upper_g_2 = re.search(re_lone_upper_2, seq)
#     lower_g_1 = re.search(re_lone_lower_1, seq)
#     lower_g_2 = re.search(re_lone_lower_2, seq)
#     format_err_1 = re.search(re_format_err_1, seq)
#     format_err_2 = re.search(re_format_err_2, seq)
#     format_err_3 = re.search(re_format_err_3, seq)
#     format_err_4 = re.search(re_format_err_4, seq)
#
#     if upper_g_1 is not None:
#         pos = seq[upper_g_1.end():].find('<')
#         pos_2 = seq[upper_g_1.end():].find('[')
#         if pos != -1 and pos < pos_2:
#             upper = find_sub_sequence(re_upper, upper_g_1.group()) + " "
#             return standardise(seq[:upper_g_1.start()] + seq[upper_g_1.end():upper_g_1.end() + pos + 1] + upper + seq[upper_g_1.end() + pos + 1:])
#         elif pos_2 != -1:
#             return standardise(seq[:upper_g_1.end() - 2] + seq[upper_g_1.end():])
#     elif upper_g_2 is not None:
#         pos = seq[:upper_g_2.start()].rfind('>')
#         pos_2 = seq[:upper_g_2.start()].rfind(']')
#         if pos != -1 and pos > pos_2:
#             upper = " " + find_sub_sequence(re_upper, upper_g_2.group())
#             return standardise(seq[:pos] + upper + seq[pos:upper_g_2.start()])
#         elif pos_2 != -1 and pos_2 > pos:
#             return standardise(seq[:pos_2 + 1] + upper_g_2.group()[2:] + seq[pos_2 + 1:upper_g_2.start()])
#     elif lower_g_1 is not None:
#         pos = seq[lower_g_1.end():].find('{')
#         pos_2 = seq[lower_g_1.end():].find('[')
#         if pos != -1 and pos < pos_2:
#             lower = find_sub_sequence(re_lower, lower_g_1.group()) + " "
#             return standardise(seq[:lower_g_1.start()] + seq[lower_g_1.end():lower_g_1.end() + pos + 1] + lower + seq[lower_g_1.end() + pos + 1:])
#         elif pos_2 != -1:
#             return standardise(seq[:lower_g_1.end() - 1] + seq[lower_g_1.end():])
#     elif lower_g_2 is not None:
#         pos = seq[:lower_g_2.start()].rfind('}')
#         pos_2 = seq[:lower_g_2.start()].rfind(']')
#         if pos != -1 and pos > pos_2:
#             lower = " " + find_sub_sequence(re_lower, lower_g_2.group())
#             return standardise(seq[:pos] + lower + seq[pos:lower_g_2.start()])
#         elif pos_2 != -1 and pos_2 > pos:
#             return standardise(seq[:lower_g_2.start()] + seq[lower_g_2.start() + 1:])
#     elif format_err_1 is not None:
#         upper = format_err_1.group(3)[1:len(format_err_1.group(3)) - 1] + " "
#         return standardise(
#             seq[:format_err_1.start(3)] + seq[format_err_1.end(3):format_err_1.start(6) + 1] + upper + seq[format_err_1.start(6) + 1:])
#     elif format_err_2 is not None:
#         lower = format_err_2.group(4)[1:len(format_err_2.group(4)) - 1] + " "
#         new = seq[:format_err_2.start(4)] + seq[format_err_2.end(4):format_err_2.start(5) + 1] + lower + seq[format_err_2.start(5) + 1:]
#         return standardise(new)
#     elif format_err_3 is not None:
#         new = seq[:format_err_3.start(3)] + seq[format_err_3.end(3):format_err_3.start(6)] + format_err_3.group(3) + seq[format_err_3.start(6):]
#         return standardise(new)
#     elif format_err_4 is not None:
#         new = seq[:format_err_4.start(4)] + ":" + format_err_4.group(4) + seq[format_err_4.end(4) + 1:]
#         return standardise(new)
#     else:
#         return seq

# re_lone_upper_1 = re.compile(f"^{re_upper.pattern}::|(?<=::){re_upper.pattern}::")
# re_lone_upper_2 = re.compile(f"::({re_upper.pattern})$")
# re_lone_lower_1 = re.compile(f"^({re_lower.pattern}):(?=[^:])|(?<=[^:]:)({re_lower.pattern}):(?=[^:])")

    # def toehold_covering(self, k):
    #     k = tidy(k)
    #     print("k", k)
    #     for gate in re.finditer(re_gate, k):
    #         d_s = gate.group(3)
    #         #d_s = re.search(re_double, gate.group())
    #         print("Here")
    #         pre_cover = re.search(re_pre_cover, gate.group())
    #         print(pre_cover, "pre")
    #         post_cover = re.search(re_post_cover, gate.group())
    #         print("post", post_cover)
    #         if pre_cover is not None:
    #             th_pos = gate.group()[pre_cover.start() + 2: gate.start(3)].find(pre_cover.group() + "^")
    #             print("th_pos", th_pos)
    #             updated_gate = gate.group()[:pre_cover.start()] + gate.group()[
    #                                                               pre_cover.end() + 2: pre_cover.end() + th_pos] + \
    #                            ">[" + pre_cover.group() + "^ " + gate.group()[d_s.start() + 1:]
    #             updated_seq = k[:gate.start()] + tidy(updated_gate) + k[gate.end():]
    #             print(updated_seq, "updated seq")
    #             yield self.Transition([k], [updated_seq], alpha)
    #         if post_cover is not None:
    #             th_c_pos = gate.group()[post_cover.end() + 1:].find(post_cover.group() + "^*")
    #             print("th_c pos", th_c_pos)
    #             updated_gate = gate.group()[:d_s.end() - 1] + " " + post_cover.group() + "]<" + \
    #                            gate.group()[post_cover.end() + 1: post_cover.end() + th_c_pos + 1] + \
    #                            gate.group()[post_cover.end() + th_c_pos + 4:]
    #             updated_seq = k[:gate.start()] + tidy(updated_gate) + k[gate.end():]
    #             print(updated_seq, "updated seq")
    #             yield self.Transition([k], [updated_seq], alpha)
