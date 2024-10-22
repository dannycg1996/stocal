"""DNA Strand Displacement Simulation

The model is based on  rules that determine the behavior of DSD behaviour (Lakin, 2012)

"""
import math
import re
import stocal
from stocal.structures import multiset

domains = {"Z": 4} # Nucleotide lengths of domains. Should be modified appropriately by the user.

re_double = re.compile(r'(\[[^<{\[\]}>]*?])')  # Matches on any double strand (includes brackets).
re_upper = re.compile(r'(<[^<\[{]*?>)')  # Matches on any upper strand (includes the brackets).
re_lower = re.compile(r'({[^<\[{]*?\})')  # Matches on any lower strand (includes the brackets).
re_short_double_th = re.compile(r'(?:\[\W*?(\w)(?:\^\W*?\]))')  # Matches on double toeholds of the form [A^] not [A^ B]
re_gate = re.compile(
    f"{re_lower.pattern}?{re_upper.pattern}?{re_double.pattern}{re_upper.pattern}?{re_lower.pattern}?")  # Matches on gates

re_double_lab = re.compile(r'(\w)(?=\^)(?=[^<>{}]*])')  # Returns the label of a double toehold regex.
re_upper_lab = re.compile(r'(\w)(?=\^)(?=[^<>]*>)')  # Returns the labels of upper toeholds.
re_lower_lab = re.compile(r'(\w)(?=\^\*)(?=[^{}]*})')  # Returns labels of lower toeholds
re_open = re.compile(r'[<\[{]')  # Matches on open brackets [, { and <
re_close = re.compile(r'([>\]}])')  # Matches on close brackets ], } and >
re_empty = re.compile(r'(<(?:\s)*>)|({(?:\s)*})|(\[(?:\s)*])')  # Matches on empty brackets like <>, {} and [ ].
re_large_spaces = re.compile(r'(\s{2,})')  # Matches on spaces of length > 1
re_spaces = re.compile(r'(?<=[:>}\]<{\[])(\s+)|(\s+)(?=[:>\]}])')  # Matches on unnecessary spaces.

# The below 4 patterns match on different variants of gates which contain just a single upper or lower strand.
re_lone_upper_1 = re.compile(f"^{re_upper.pattern}::{re_gate.pattern}|(?<=::){re_upper.pattern}::{re_gate.pattern}")
re_lone_upper_2 = re.compile(f"{re_gate.pattern}::({re_upper.pattern})$")
re_lone_lower_1 = re.compile(f"^{re_lower.pattern}:{re_gate.pattern}|(?<=[^:]:){re_lower.pattern}:{re_gate.pattern}")
re_lone_lower_2 = re.compile(f"{re_gate.pattern}:{re_lower.pattern}$")

#  Matches where the Covering rule can be applied on a gate, before the d_s
re_pre_cover = re.compile(fr"{{([^}}]*?)(\w+)\^\*\s*}}<([^>]*)(\2)\^\s*>")
# Matches where the Covering rule can be applied on a gate, after the d_s (or across two gates)
re_post_cover = re.compile(
    fr"<(\w+)\^\s*([^>]*)>(:?){{(\1)\^\*\s*([^}}]*)}}|(\w+)\^\*([^}}]*)}}::{re_lower.pattern}?<(\6)\^\s*([^>]*)>")
# Matches where upper strand migration can occur (left to right).
re_upper_migrate = re.compile(fr"{re_double.pattern}(<(\w+)[^<>:]*?>):{re_upper.pattern}?(\[(\3)(?!\])[^\]]*?\])")
# Matches where lower strand migration can occur (left to right).
re_lower_migrate = re.compile(fr"{re_double.pattern}({{(\w+)[^{{}}:]*?}})::{re_lower.pattern}?(\[(\3)(?!\])[^\]]*?\])")
# Matches where upper strand rev migration can occur (right to left).
re_upper_migrate_r = re.compile(fr"(\[[^\]]*(?<=\s)(\w+)\]){re_upper.pattern}?:(<[^<>:]*?(\2)>){re_double.pattern}")
# Matches where lower strand rev migration can occur (right to left).
re_lower_migrate_r = re.compile(fr"(\[[^<>]*(?<=\s)(\w+)\]){re_lower.pattern}?::({{[^<>:]*?(\2)}}){re_double.pattern}")

# Matches where upper strand displacement (left to right) can occur.
re_displace_upper = re.compile(
    fr"{re_double.pattern}<(\w+)([^<>:]*?)>:{re_upper.pattern}?\[(\2)\]{re_upper.pattern}?{re_lower.pattern}?")
# Matches where lower strand displacement (left to right) can occur.
re_displace_lower = re.compile(
    fr"{re_double.pattern}{{(\w+)([^{{}}:]*?)}}::{re_lower.pattern}?\[(\2)\]{re_upper.pattern}?{re_lower.pattern}?")
# Matches where upper strand displacement (left to right) can occur.
re_displace_upper_r = re.compile(
    fr"{re_lower.pattern}?{re_upper.pattern}?\[(\w+)\]{re_upper.pattern}?:<([^<>:]*?)(\3)>{re_double.pattern}")
# Matches where lower strand displacement (left to right) can occur.
re_displace_lower_r = re.compile(
    fr"{re_lower.pattern}?{re_upper.pattern}?\[(\w+)\]{re_lower.pattern}?::{{([^<>:]*?)(\3)}}{re_double.pattern}"
)

# The four regexes below match non-normalised patterns. For example, {A}<B>[C]<D>{E}::{F}<G>[H] is not normalised, and must be
# rewritten as {A}<B>[C]{E}::{F}<D G>[H] to ensure that reactions are reversible and the results are clear
re_format_1 = re.compile(
    f"({re_double.pattern}{re_upper.pattern}{re_lower.pattern}?::{re_lower.pattern}?{re_upper.pattern}{re_double.pattern})")
re_format_2 = re.compile(
    f"({re_double.pattern}{re_upper.pattern}?{re_lower.pattern}:{re_lower.pattern}{re_upper.pattern}?{re_double.pattern})")
re_format_3 = re.compile(
    f"({re_double.pattern}{re_upper.pattern}{re_lower.pattern}?::{re_lower.pattern}?{re_double.pattern})")
re_format_4 = re.compile(
    f"({re_double.pattern}{re_upper.pattern}?{re_lower.pattern}:{re_upper.pattern}?{re_double.pattern})")

re_double_start_leak = re.compile(r'\[((\w\^)([^<\[{]*?\w))\]')
re_double_end_leak = re.compile(r'\[((\w[^<\[{]*?)\s(\w\^))\]')

unbinding_rate =  0.1126  # Rate parameter for the unbinding rule
covering_rate = 699  # Rate parameter for the covering rule
leak_rate = 0.000003  # Rate parameter for the two leakage rules


def check_in(seq):
    """Takes a sequence, and if it isn't None it returns the sequence with the first and last character missing (this will be used to remove
     brackets around a group and make code more readable). If seq == None, return a blank string '' """
    if seq is not None:
        return seq[1:len(seq) - 1]
    return ""


def check_out(seq):
    """Takes a sub sequence, and either returns the regex match (if seq!=None) or a blank string '' """
    if seq is not None:
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
        strand = sys[:match_obj.start()] + sys[match_obj.end(1+i)+2:match_obj.start(3+i)+1] +\
            match_obj.group(1+i)[1:len(match_obj.group(1+i))-1] + " " + sys[match_obj.start(3+i)+1:]
    elif match_obj.group(2+i) is not None:
        strand = sys[:match_obj.start()] + match_obj.group(2+i) + match_obj.group(1+i) + sys[match_obj.start(4+i):]
    else:
        strand = sys[:match_obj.start()] + match_obj.group(1+i) + sys[match_obj.start(4+i):]
    return strand


def fix_lower_gate(sys, match_obj, i):
    """This function takes a system sys, a match object and a starting index. The match object identifies gates which consist solely of
     a lower strand, merges it with a gate to the right, and then returns the updated system"""
    if match_obj.group(2+i) is not None:  # Match object has 6 groups: ({ })::({ })(< >)([ ])(< >)({ })
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
        else:  # If 2nd match condition of upper_g_1 is met.
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
        if lower_g_2.group(5) is not None:  # If gate before the lower strand gate had a lower strand after the double strand
            strand = sys[:lower_g_2.end(5)-1] + " " + lower_g_2.group(6)[1:]
        else:
            strand = sys[:lower_g_2.start(6)-1] + lower_g_2.group(6)
        return merge_gates(strand)
    else:
        return sys


def reformat(sys):
    """This function identifies non-standard patterns and re-formats it. For example, {A}<B>[C]<D>{E}::{F}<G>[H] must be rewritten
    as {A}<B>[C]{E}::{F}<D G>[H] to ensure that the reaction is reversible and the results are clear"""
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


def rotate(strand):
    """Takes a single upper or lower strand, and rotates it to be an upper or lower strand, respectively.
    Rotation is performed as defined in Lakin's DSD calculus"""
    if re.search(re_gate, strand) is None:  # Check that the input is not a gate (this program does not rotate gates)
        domains = check_in(strand).split(" ")[::-1]  # Remove brackets and reverse the domain sequence.
        new_strand = ""  # Define empty strand.
        if re.search(re_upper, strand) is not None:  # if upper strand, loop through reversed input and build lower strand.
            for domain in domains:
                new_strand = new_strand + " " + domain
            return "{" + new_strand[1:] + "}"
        elif re.search(re_lower, strand) is not None:  # if lower strand, loop through reversed input and build upper strand.
            for domain in domains:
                new_strand = new_strand + " " + domain
            return "<" + new_strand[1:] + ">"
    else:
        return ""


def convert_upper_to_lower(strand):
    """Take an upper strand sequence, and convert it to be a lower strand sequence.
    This essentially consists of swapping domain names A to be domain name A*"""
    return re.sub(r'(?<=\S)\s', "* ", strand) + "*"


def convert_lower_to_upper(strand):
    """Take a lower strand sequence, and convert it to be an upper strand sequence.
    This essentially consists of swapping domain names A* to be domain name A"""
    return re.sub(r'\*', "", strand)


def get_binding_rate(t_h_label):
    """Calculate the binding rate for a given toehold. Calculates it based on nucleotide length.
    # If nucleotide length is unknown, then the toehold length of 7 is used, which is an average toehold length."""
    for label in domains:
        if t_h_label == label:
            if domains[label] < 5:
                if domains[label] == 4:
                    return math.log10(5)/2500
                else:
                    return math.log10(domains[label])/2500
    return math.log10(6)/2500


def get_migration_rate(domain_label):
    """Calculate the migration/displacement rate for a given domain. Calculates it based on nucleotide length.
    If nucleotide length is unknown, then the domain length of 20 is used, which is the shortest length a
    long domain can be."""
    nuc_length = 20
    for label in domains:
        if domain_label == label:
            nuc_length = domains[label]
    return 8000/(nuc_length**2)


class BindingRule(stocal.TransitionRule):
    """Join any two strings into their concatenations"""
    Transition = stocal.MassAction

    def novel_reactions(self, k, l):
        gate_k = re.search(re_gate, k)
        gate_l = re.search(re_gate, l)
        # Call the appropriate function depending if k and l are both strands, or a gate and a strand.
        if (gate_k is None and gate_l is not None) or (gate_l is None and gate_k is not None):
            yield from self.strand_to_gate_binding(k, l, re_upper_lab, re_lower_lab)
            yield from self.strand_to_gate_binding(l, k, re_upper_lab, re_lower_lab)
            yield from self.strand_to_gate_binding(k, rotate(l), re_upper_lab, re_lower_lab)
            yield from self.strand_to_gate_binding(l, rotate(k), re_upper_lab, re_lower_lab)
            yield from self.strand_to_gate_binding(k, l, re_lower_lab, re_upper_lab)
            yield from self.strand_to_gate_binding(l, k, re_lower_lab, re_upper_lab)
            yield from self.strand_to_gate_binding(k, rotate(l), re_lower_lab, re_upper_lab)
            yield from self.strand_to_gate_binding(l, rotate(k), re_lower_lab, re_upper_lab)
        elif gate_k is None or gate_l is None:
            yield from self.strand_to_strand_binding(k, l, re_upper_lab, re_lower_lab)
            yield from self.strand_to_strand_binding(k, l, re_lower_lab, re_upper_lab)
            yield from self.strand_to_strand_binding(rotate(k), l, re_upper_lab, re_lower_lab)
            yield from self.strand_to_strand_binding(rotate(k), l, re_lower_lab, re_upper_lab)

    def strand_to_gate_binding(self, k, l, regex_1, regex_2):
        """Simulates binding between a gate and a single upper or lower strand"""
        for gate in re.finditer(re_gate, k):   # Loop through the gates in system k.
            # The next two for loops attempt to find matching upper and lower toeholds on the gate and strand.
            for match in re.finditer(regex_1, gate.group()):
                for match_2 in re.finditer(regex_2, l):
                    if match.group() == match_2.group(): # If matching toeholds are found
                        binding_rate = get_binding_rate(match.group())
                        d_s = "[" + match.group() + "^]"
                        i = gate.start()
                        if regex_1 == re_upper_lab:
                            l_s_1 = "{" + l[1:match_2.start()] + "}"
                            l_s_2 = "{" + l[match_2.end() + 2:len(l) - 1] + "}"
                            if match.start() > gate.start(2) - i and match.end() < gate.end(2) - i:
                                u_s_1 = "<" + k[gate.start(2) + 1:match.start() + i] + ">"
                                u_s_2 = "<" + k[match.end() + 1 + i:gate.end(2) - 1] + ">"
                                sys = k[:gate.start()] + l_s_1 + u_s_1 + d_s + l_s_2 + "::" + gate.group(1) + u_s_2 + k[gate.start(3):]
                                yield self.Transition([k, l], [standardise(sys)], binding_rate)
                            elif match.start() > gate.start(4) - i and match.end() < gate.end(4) - i:
                                u_s_1 = "<" + k[gate.start(4) + 1:match.start() + i] + ">"
                                u_s_2 = "<" + k[match.end() + i + 1:gate.end(4) - 1] + ">"
                                sys = k[:gate.end(3)] + check_out(gate.group(5)) + "::" + l_s_1 + u_s_1 + d_s + u_s_2 + l_s_2 + k[gate.end():]
                                yield self.Transition([k, l], [standardise(sys)], binding_rate)
                        else:
                            u_s_1 = "<" + l[1:match_2.start()] + ">"
                            u_s_2 = "<" + l[match_2.end() + 2:len(l) - 1] + ">"
                            if match.start() > gate.start(1) - i and match.end() < gate.end(1) - i:
                                l_s_1 = "{" + k[gate.start(1) + 1:match.start() + i] + "}"
                                l_s_2 = "{" + k[match.end() + i + 2:gate.end(1) - 1] + "}"
                                sys = k[:gate.start()] + l_s_1 + u_s_1 + d_s + u_s_2 + l_s_2 + ":" + check_out(gate.group(2)) + k[gate.start(3):]
                                yield self.Transition([k, l], [standardise(sys)], binding_rate)
                            elif match.start() > gate.start(5) - i and match.end() < gate.end(5) - i:
                                l_s_1 = "{" + k[gate.start(5) + 1:match.start() + i] + "}"
                                l_s_2 = "{" + k[match.end() + i + 2:gate.end(5) - 1] + "}"
                                sys = k[:gate.end(3)] + check_out(gate.group(4)) + ":" + l_s_1 + u_s_1 + d_s + u_s_2 + l_s_2 + k[gate.end():]
                                yield self.Transition([k, l], [standardise(sys)], binding_rate)

    def strand_to_strand_binding(self, k, l, regex_1, regex_2):
        """Simulates an upper and lower strand annealing together"""
        # The next two loops are to loop through matching toeholds found on the two strands.
        for match_1 in re.finditer(regex_1, k):
            for match_2 in re.finditer(regex_2, l):
                if match_1.group() == match_2.group():
                    binding_rate = get_binding_rate(match_1.group())
                    d_s = "[" + match_2.group() + "^]"
                    part_a = l[:match_2.start()] + re.search(re_close, l[match_2.start():]).group()
                    part_b = k[:match_1.start()] + re.search(re_close, k[match_1.start():]).group()
                    part_c = re.search(re_open, k[:match_1.end() + 1]).group()
                    part_d = re.search(re_open, l[:match_2.end()]).group()
                    if regex_1 == re_upper_lab:
                        sys = part_a + part_b + d_s + part_c + k[match_1.end() + 1:] + part_d + l[match_2.end() + 2:]
                    else:
                        sys = part_b + part_a + d_s + part_d + l[match_2.end() + 1:] + part_c + k[match_1.end() + 2:]
                    yield self.Transition([k, l], [tidy(sys)], binding_rate)


class UnbindingRule(stocal.TransitionRule):
    """Splits a system into two systems when a toehold unbinds"""
    Transition = stocal.MassAction

    def novel_reactions(self, kl):
        yield from self.toehold_unbinding(kl)

    def toehold_unbinding(self, kl):
        """This function loops through a system gate by gate, and identifies double strands which can be unbound i.e.
        double strands of the form [A^]. It then yields the two separate parts, which would be produced when that double strand
        (toehold) unbound."""
        for gate in re.finditer(re_gate, kl):  # Loop through the system gate by gate.
            d_s = re.search(re_short_double_th, gate.group())  # If one exists, retrieve the unbindable double strand in the gate.
            if d_s is not None:
                label = re.search(re_double_lab, d_s.group()).group()  # Retrieve label of unbindable toehold.
                part_a = "<" + check_in(gate.group(2)) + " " + label + "^ " + check_in(gate.group(4)) + ">"  # Build upper part of gate.
                part_b = "{" + check_in(gate.group(1)) + " " + label + "^* " + check_in(gate.group(5)) + "}"  # Build lower part pf hate
                # Assemble the gates with the rest of the system, depending on how the gates were connected.
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
                yield self.Transition([kl], [standardise(part_a), standardise(part_b)], unbinding_rate)


class CoveringRule(stocal.TransitionRule):
    """This rule carries out the toehold covering reaction, where an exposed toehold in a lower strand is covered by a complementary
     exposed toehold in the upper strand"""
    Transition = stocal.MassAction

    def novel_reactions(self, k):
        yield from self.toehold_covering(k)

    def toehold_covering(self, k):
        for match in re.finditer(re_post_cover, k):  # Match on <>{} or <>:{} or {}::{}?<> sequences where Covering can be applied.
            if match.group(1) is not None:  # If matching on <>{} or <>:{} then apply covering to system.
                updated_sys = k[:match.start()-1] + " " + match.group(1) + "^]<" + check_out(match.group(2)) + ">" + \
                    check_out(match.group(3)) + "{" + check_out(match.group(5)) + "}" + k[match.end():]
            else:  # If matching on {}::{}?<> then update system.
                updated_sys = k[:match.start()-2] + " " + match.group(6) + "^]{" + check_out(match.group(7)) + "}::" + \
                    check_out(match.group(8)) + "<" + check_out(match.group(10)) + ">" + k[match.end():]
            yield self.Transition([k], [tidy(updated_sys)], covering_rate)
        for match in re.finditer(re_pre_cover, k):  # Match on {}<> sequences where Covering can be applied.
            updated_sys = k[:match.start()] + "{" + check_out(match.group(1)) + "}<" + check_out(match.group(3)) + ">[" + \
                match.group(2) + "^ " + k[match.end()+1:]
            yield self.Transition([k], [tidy(updated_sys)], covering_rate)


class MigrationRule(stocal.TransitionRule):
    """Migrates an upper or lower overhang up/down a strand via branch migration"""
    Transition = stocal.MassAction

    def novel_reactions(self, k):
        k = tidy(k)
        yield from self.migrate(k, re_lower_migrate, re_lower)
        yield from self.migrate(k, re_upper_migrate, re_upper)
        yield from self.migrate_rev(k, re_lower_migrate_r, re_lower)
        yield from self.migrate_rev(k, re_upper_migrate_r, re_upper)

    def migrate(self, k, regex_1, regex_2):
        for match in re.finditer(regex_1, k):
            migration_rate = get_migration_rate(match.group(3))
            i = match.start()
            d_s_1 = match.group(1)[:len(match.group(1))-1] + " " + match.group(3) + "]"
            d_s_2 = "[" + match.group()[match.end(6)-i:match.end(5)-i]
            if regex_2 == re_lower:
                strand_1 = "{" + match.group()[match.end(3)-i:match.end(2)-i]
                strand_2 = "{" + check_in(match.group(4)) + " " + match.group(3) + "}"
                bracket = "::"
            else:
                strand_1 = "<" + match.group()[match.end(3)-i:match.end(2)-i]
                strand_2 = "<" + check_in(match.group(4)) + " " + match.group(3) + ">"
                bracket = ":"
            seq = tidy(k[:match.start()] + d_s_1 + strand_1 + bracket + strand_2 + d_s_2 + k[match.end():])
            yield self.Transition([k], [seq], migration_rate)

    def migrate_rev(self, k, regex_1, regex_2):
        for match in re.finditer(regex_1, k):
            migration_rate = get_migration_rate(match.group(2))
            i = match.start()
            d_s_1 = match.group()[:match.start(2)-i] + "]"
            d_s_2 = "[" + match.group(2) + " " + match.group(6)[1:]
            if regex_2 == re_lower:
                strand_1 = "{" + match.group(2) + " " + check_in(match.group(3)) + "}"
                strand_2 = match.group()[match.start(4)-i: match.start(5)-i] + "}"
                bracket = "::"
            else:
                strand_1 = "<" + match.group(2) + " " + check_in(match.group(3)) + ">"
                strand_2 = match.group()[match.start(4)-i: match.start(5)-i] + ">"
                bracket = ":"
            seq = tidy(k[:match.start()] + d_s_1 + strand_1 + bracket + strand_2 + d_s_2 + k[match.end():])
            yield self.Transition([k], [seq], migration_rate)


class DisplacementRule(stocal.TransitionRule):
    """Splits two strings when one strand displaces another"""
    Transition = stocal.MassAction

    def novel_reactions(self, k):
        k = tidy(k)
        yield from self.displacement_fwd(k, re_displace_upper)
        yield from self.displacement_fwd(k, re_displace_lower)
        yield from self.displacement_rev(k, re_displace_upper_r)
        yield from self.displacement_rev(k, re_displace_lower_r)

    def displacement_fwd(self, k, regex_1):
        for match in re.finditer(regex_1, k):
            displacement_rate = get_migration_rate(match.group(2))
            strand_1 = check_in(match.group(4)) + " " + match.group(2) + " "
            start = k[:match.end(1)-1] + " " + match.group(2) + "]"
            if regex_1 == re_displace_upper:
                if k[match.end():match.end()+2] != "::":
                    strand_1 = tidy("<" + strand_1 + check_in(match.group(6)) + ">")
                    strand_2 = tidy(start + "<" + check_out(match.group(3)) + ">" + check_out(match.group(7)) + k[match.end():])
                else:
                    strand_1 = tidy(start + "<" + check_out(match.group(3)) + ">" + check_out(match.group(7)))
                    strand_2 = tidy("<" + check_in(match.group(4)) + " " + match.group(2) + ">" + k[match.end()+2:])
            else:
                if k[match.end():match.end()+1] == ":" and k[match.end()+1:match.end()+2] != ":":
                    strand_1 = tidy(start + check_out(match.group(6)) + "{" + check_out(match.group(3)) + "}")
                    strand_2 = tidy("{" + check_in(match.group(4)) + " " + match.group(5) + "}" + k[match.end()+1:])
                else:
                    strand_1 = tidy("{" + strand_1 + check_in(match.group(7)) + "}")
                    strand_2 = tidy(start + " " + check_out(match.group(6)) + "{" + check_out(match.group(3)) + "}" + k[match.end():])
            yield self.Transition([k], [strand_1, strand_2], displacement_rate)

    def displacement_rev(self, k, regex_1):
        for match in re.finditer(regex_1, k):
            displacement_rate = get_migration_rate(match.group(3))
            if regex_1 == re_displace_upper_r:
                if k[match.start()-2:match.start()] != "::":
                    strand_1 = "<" + check_in(match.group(2)) + " " + match.group(3) + " " + check_in(match.group(4)) + ">"
                    strand_2 = k[:match.start()] + check_out(match.group(1)) + "<" + check_out(match.group(5)) + ">[" + match.group(3) + " " + match.group(7)[1:] + k[match.end():]
                else:
                    strand_1 = tidy(k[:match.start()-2]) + "<" + match.group(3) + " " + check_in(match.group(4)) + ">"
                    strand_2 = check_out(match.group(1)) + "<" + match.group(5) + ">[" + match.group(3) + " " + check_in(match.group(7)) + "]" + k[match.end(7):]
            else:
                if k[match.start()-1:match.start()] == ":" and k[match.start()-2:match.start()-1] != ":":
                    strand_1 = tidy(k[:match.start()-1] + "{" + match.group(3) + " " + check_in(match.group(4)) + "}")
                    strand_2 =  "{" + match.group(5) + "}" + check_out(match.group(2)) +"[" + match.group(3) + " " + check_in(match.group(7)) + "]" + k[match.end(7):]
                else:
                    strand_1 = "{" + check_in(match.group(1)) + " " + match.group(3) + " " + check_in(match.group(4)) + "}"
                    strand_2 = k[:match.start()] + "{" + check_out(match.group(5)) + "}" + check_out(match.group(2)) + "[" + \
                        match.group(3) + " " + match.group(7)[1:] + k[match.end():]
            yield self.Transition([k], [tidy(strand_1), tidy(strand_2)], displacement_rate)


class StrandLeakageRule(stocal.TransitionRule):
    """Simulate leak reactions on a double-stranded complex"""
    Transition = stocal.MassAction

    def novel_reactions(self, k, l):
        gate_k = re.search(re_gate, k)
        gate_l = re.search(re_gate, l)
        if (gate_k is None and gate_l is not None) or (gate_l is None and gate_k is not None):
            yield from self.strand_leak(k, l)
            yield from self.strand_leak(l, k)

    def upper_strand_leakage(self, k, l, mod_l, gate):
        leaked_u_s = "<" + check_in(gate.group(2)) + " " + check_in(gate.group(3)) + " " + check_in(gate.group(4)) + ">"
        re_strand = re.sub(r'\^', "\\^", check_in(gate.group(3)))
        re_strand_2 = re_strand + "$|" + re_strand + " "
        for match in re.finditer(fr'{re_strand_2}', mod_l):  # Yield suitable (upper) leaks.
            new_sys = k[:gate.start()] + check_out(gate.group(1)) + "<" + mod_l[:match.start()] + ">" + gate.group(3) + "<" + \
                      mod_l[match.end():] + ">" + check_out(gate.group(5)) + k[gate.end():]
            yield self.Transition([k, l], [tidy(new_sys), tidy(leaked_u_s)], leak_rate)

    def lower_strand_leakage(self, k, l, mod_l, gate):
        re_strand = re.sub(r'\^', "\\^", check_in(gate.group(3)))
        re_strand = re.sub(r'(?<=\S)\s', "\* ", re_strand) + "\*"
        leaked_l_s = "{" + check_in(gate.group(1)) + " " + convert_upper_to_lower(check_in(gate.group(3))) + \
                     " " + check_in(gate.group(5)) + "}"
        for match in re.finditer(fr'{re_strand}', mod_l): # Yield suitable (lower) leaks.
            new_sys = k[:gate.start()] + "{" + mod_l[:match.start()] + "}" + k[gate.start(2):gate.end(4)] +\
              "{" + mod_l[match.end():] + "}" + k[gate.end():]
            yield self.Transition([k, l], [tidy(new_sys), tidy(leaked_l_s)], leak_rate)

    def strand_leak(self, k, l):
        for gate in re.finditer(re_gate, k):
            if re.search(re_short_double_th, gate.group(3)) is None:  # Checks that the d_s in the gate is not of the form [A^]
                upper_gate_join_1 = k[gate.start()-2:gate.start()]  # Used to check if current gate joins last gate via an upper strand.
                upper_gate_join_2 = k[gate.end():gate.end()+2]  # Used to check if current gate joins next gate via an upper strand.
                lower_gate_join_1 = k[gate.start() - 2:gate.start() - 1]  # Used to check if current gate joins last gate via a lower strand.
                lower_gate_join_2 = k[gate.end() + 1:gate.end() + 2]  # Used to check if current gate joins next gate via a lower strand.
                if re.search(re_upper, l) is not None:  # If the strand initiating the leak is an upper strand:
                    if upper_gate_join_1 != "::" and upper_gate_join_2 != "::":  # Check gate isn't joined to others by upper strand.
                        yield from self.upper_strand_leakage(k, l, check_in(l), gate)
                    if lower_gate_join_1 != ":" and lower_gate_join_2 != ":":  # Check gate isn't joined to others by lower strand.
                        yield from self.lower_strand_leakage(k, l, check_in(rotate(l)), gate)
                else:  # If the strand initiating the leak is a lower strand:
                    if lower_gate_join_1 != ":" and lower_gate_join_2 != ":":  # Check gate isn't joined to others by lower strand.
                        yield from self.lower_strand_leakage(k, l, check_in(l), gate)
                    if upper_gate_join_1 != "::" and upper_gate_join_2 != "::":   # Check gate isn't joined to others by upper strand.
                        yield from self.upper_strand_leakage(k, l, check_in(rotate(l)) , gate)


class ToeholdLeakageRule(stocal.TransitionRule):
    """Simulates a leak on a strand where a toehold has to spontaneously unbind for the leak to occur"""
    Transition = stocal.MassAction

    def novel_reactions(self, k, l):
        gate_k = re.search(re_gate, k)
        gate_l = re.search(re_gate, l)
        if (gate_k is None and gate_l is not None) or (gate_l is None and gate_k is not None):
            yield from self.toehold_leak(k, l)
            yield from self.toehold_leak(l, k)

    def lower_toehold_leakage_at_end(self, k, l, end_leak, mod_l, gate):
        re_check_not_l_s = "^" + re.sub(r'\^', "\\^", end_leak.group(3))
        re_end_leak = convert_upper_to_lower(re.sub(r'\^', "\\^", end_leak.group(2)))
        re_leak = re.sub(r'\*', "\\*", re_end_leak)
        for match in re.finditer(re_leak, mod_l):
            if re.search(re_check_not_l_s, l[match.end():]) is None:
                leaked_l_s = "{" + check_in(gate.group(1)) + " " + convert_upper_to_lower(end_leak.group(1)) +\
                                 " " + check_in(gate.group(5)) + "}"
                new_sys = k[:gate.start()] + "{" + mod_l[:match.start()] + "}" + gate.group(2) + "[" + end_leak.group(2) + "]<" + \
                        end_leak.group(3) + " " + check_in(gate.group(4)) + ">{" + mod_l[match.end():] + "}" + k[gate.end():]
                yield self.Transition([k, l], [tidy(leaked_l_s), tidy(new_sys)], leak_rate)

    def upper_toehold_leakage_at_end(self, k, l, end_leak, mod_l, gate):
        re_check_not_l_s = "^" + re.sub(r'\^', "\\^", end_leak.group(3))
        re_end_leak = re.sub(r'\^', "\\^", end_leak.group(2))
        re_end_leak_2 = re_end_leak + "$|" + re_end_leak + " "
        for match in re.finditer(re_end_leak_2, mod_l):
            if re.search(re_check_not_l_s, l[match.end():]) is None:
                leaked_u_s = "<" + check_in(gate.group(2)) + " " + end_leak.group(1) + " " + check_in(gate.group(4)) + ">"
                new_sys = k[:gate.start(2)] + "<" + mod_l[:match.start()] + ">[" + end_leak.group(2) + "]<" + \
                    mod_l[match.end():] + ">{" + end_leak.group(3) + "* " + check_in(gate.group(5)) + "}" + k[gate.end():]
                yield self.Transition([k, l], [tidy(leaked_u_s), tidy(new_sys)], leak_rate)

    def lower_toehold_leakage_at_start(self, k, l, start_leak, mod_l, gate):
        re_check_not_l_s = re.sub(r'\^', "\\^", start_leak.group(2)) + "$"
        re_start_leak = convert_upper_to_lower(re.sub(r'\^', "\\^", start_leak.group(3)))
        re_leak = re.sub(r'\*', "\\*", re_start_leak)
        for match in re.finditer(re_leak, mod_l):
            if re.search(re_check_not_l_s, l[match.end():]) is None:
                leaked_l_s = "{" + check_in(gate.group(1)) + " " + convert_upper_to_lower(start_leak.group(1)) +\
                                 " " + check_in(gate.group(5)) + "}"
                new_sys = k[:gate.start()] + "{" + mod_l[:match.start()] + "}<" + check_in(gate.group(2)) + " " + \
                          start_leak.group(2) + ">[" + start_leak.group(3) + "]<" + check_in(gate.group(4)) + ">" + \
                          "{" + mod_l[match.end():] + "}" + k[gate.end():]
                yield self.Transition([k, l], [tidy(leaked_l_s), tidy(new_sys)], leak_rate)

    def upper_toehold_leakage_at_start(self, k, l, start_leak, mod_l, gate):
        re_check_not_l_s = re.sub(r'\^', "\\^", start_leak.group(2)) + "$"
        re_start_leak = re.sub(r'\^', "\\^", start_leak.group(3))
        re_start_leak_2 = re_start_leak + "$|" + re_start_leak + " "
        for match in re.finditer(re_start_leak_2, mod_l):
            if re.search(re_check_not_l_s, mod_l[:match.start()]) is None:  # TODO: Check this check works
                leaked_u_s = "<" + check_in(gate.group(2)) + " " + start_leak.group(1) + " " + check_in(gate.group(4)) + ">"
                new_sys = k[:gate.start()] + "{" + check_in(gate.group(1)) + " " + start_leak.group(2) + "*}<" +\
                          mod_l[:match.start()] + ">[" + start_leak.group(3) + "]<" + mod_l[match.end():] + ">" + k[gate.end(4):]
                yield self.Transition([k, l], [tidy(leaked_u_s), tidy(new_sys)], leak_rate)

    def toehold_leak(self, k, l):
        for gate in re.finditer(re_gate, k):
            start_leak = re.search(re_double_start_leak, gate.group())
            end_leak = re.search(re_double_end_leak, gate.group())
            upper_gate_join_1 = k[gate.start()-2:gate.start()]  # Used to check if current gate joins last gate via an upper strand.
            upper_gate_join_2 = k[gate.end():gate.end()+2]  # Used to check if current gate joins next gate via an upper strand.
            lower_gate_join_1 = k[gate.start() - 2:gate.start() - 1]  # Used to check if current gate joins last gate via a lower strand.
            lower_gate_join_2 = k[gate.end() + 1:gate.end() + 2]  # Used to check if current gate joins next gate via a lower strand.
            if re.search(re_upper, l) is not None:
                if upper_gate_join_1 != "::" and upper_gate_join_2 != "::":   # Check gate isn't joined to others by upper strand.
                    if start_leak is not None:
                        yield from self.upper_toehold_leakage_at_start(k, l, start_leak, check_in(l), gate)
                        if lower_gate_join_1 != ":" and lower_gate_join_2 != ":":
                            yield from self.lower_toehold_leakage_at_start(k, l, start_leak, check_in(rotate(l)), gate)
                    if end_leak is not None:  # If the strand initiating the leak is an upper strand:
                        yield from self.upper_toehold_leakage_at_end(k, l, end_leak, check_in(l), gate)
                        if lower_gate_join_1 != ":" and lower_gate_join_2 != ":":
                            yield from self.lower_toehold_leakage_at_end(k, l, end_leak, check_in(rotate(l)), gate)
            else:
                if lower_gate_join_1 != ":" and lower_gate_join_2 != ":": # Check gate isn't joined to others by upper strand.
                    if start_leak is not None:
                        yield from self.lower_toehold_leakage_at_start(k, l, start_leak, check_in(l), gate)
                        if upper_gate_join_1 != "::" and upper_gate_join_2 != "::":
                            yield from self.upper_toehold_leakage_at_start(k, l, start_leak, check_in(rotate(l)), gate)
                    if end_leak is not None:  # If the strand initiating the leak is an upper strand:
                        yield from self.lower_toehold_leakage_at_end(k, l, end_leak, check_in(l), gate)
                        if upper_gate_join_1 != "::" and upper_gate_join_2 != "::":
                            yield from self.upper_toehold_leakage_at_end(k, l, end_leak, check_in(rotate(l)), gate)

# process contains the rules which should be applied during the simulation.
process = stocal.Process(
    rules=[MigrationRule(), BindingRule(), UnbindingRule(), DisplacementRule()] # StrandLeakageRule(), ToeholdLeakageRule() CoveringRule()
)


def every(trajectory, dt):
    tmax = trajectory.tmax
    trajectory.tmax = trajectory.time
    while trajectory.time < tmax:
        transitions = multiset()
        if trajectory.steps and trajectory.step >= trajectory.steps:
            break
        trajectory.tmax += dt
        for trans in trajectory:
            transitions += {trans:1}
        yield transitions


def sample(trajectory, species, dt=0):
    """Sample species along trajectory every dt time units

     species is a list of species labels that should be sampled.

     Returns a tuple of two elements, the first is a list of all firing
     times, the second a dictionary that holds for each species the
     list of copy numbers at each corresponding time point. If dt is
     given, it specifies the interval at which the trajectory is sampled.
    """
    times = [trajectory.time]
    numbers = {s:[trajectory.state[s]] for s in species}
    it = every(trajectory, dt) if dt else iter(trajectory)
    for _ in it:
        times.append(trajectory.time)
        for s in species:
            numbers[s].append(trajectory.state[s])
    return times, numbers


class Trajectory(object):
    def __init__(self):
        self.records = []
        self.null = 0

    def append(self, **values):
        self.records.append(values)

    def __iter__(self):
        return iter(self.records)

    def __getattr__(self, value):
        from collections import Mapping, Sequence
        proto = self.records[0][value]

        if isinstance(proto, Mapping):
            keys = set().union(*(rec[value] for rec in self.records))
            return {
                key: [rec[value].get(key, self.null) for rec in self.records]
                for key in keys
            }

        elif isinstance(proto, Sequence):
            return [
                [rec[value][i] for rec in self.records]
                for i in proto
            ]

        else:
            return [rec[value] for rec in self.records]


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # Set initial species within the simulation, as well as the initial no. of molecules of each species.
    initial_state = {"<t^ b>":10, "{t^*}[x t^]:[b t^]:[a t^]:[a]":40, "[x]:[t^ b]:[t^ b]:[t^ a]{t^*}":40,
                     "<t^ x>":40, "<t^ a>":40, "<b t^>":40}
    species = ["<t^ b>", "{t^*}[x t^]:[b t^]:[a t^]:[a]", "[x]:[t^ b]:[t^ b]:[t^ a]{t^*}", "<t^ x>", "<t^ a>", "<b t^>", "<a>"]
    initial_state = {standardise(key): value for key, value in initial_state.items()} # Normalise species descriptions

    # Set duration of simulation and timesteps
    time,species = sample(
        process.trajectory(initial_state, tmax=1200),
        species,
        dt=24.
    )

    for name in species:
        plt.plot(time, species[name], label=str(name))
    plt.legend(loc='upper left')
    plt.show()

#  Below code allows for simulation without any graphing - much quicker runtimes
#     traj = process.sample(initial_state, tmax=1000.)
#     for _ in traj:
#         print(traj.time, traj.state)
