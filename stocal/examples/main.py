from stocal import *
from Rules import *
import re

# Two reactions: Two molecules of A forming a dimer A2 and reverse
from stocal.structures import DNASystem

r1 = MassAction({'A': 2}, {'A2': 1}, 1.0)
r2 = MassAction({'A2': 1}, {'A': 2}, 10.0)  # Last value is the stochastic rate at which the reaction occurs
# process = Process([r1, r2])

# Here, the event feed will occur, and feed an A molecule to the system, at time t=0.0 and then periodically every 1.0 time units.
# Unlike stochastic reactions that occur with an average frequency, nondeterministic events happen at exactly the specified times.
feed = Event([], ['A'], 0.0, 1.0)

process = Process(transitions=[r1, r2], rules=[Dilution()])

def regex_match(dna, category):
    if category == "lower_th":
        return re.findall('{[^<>\[\]]*?\s(\w)\^\s+?[^<>\[\]]*?}', dna)
    elif category == "lower_th_c":
        return re.findall('{[^<>\[\]]*?\s(\w)\^\*+?[^<>\[\]]*?}', dna)
    elif category == "upper_th":
        return re.findall('<[^{}\[\]]*?(\w?)\^\s+?[^{}\[\]]*?>', dna)
    elif category == "upper_th_c":
        return re.findall('<[^{}\]\[]*?\s+?(\w)\^\*+?[^{}\[\]]*?>', dna)
    else:
        print("Erroneous input into strand_regex method")

#def simulate_strand(dict): #a dictionary of dna is inputted

    # if len(systems)>1:
    #     r = list(range(0,len(systems)))
    #     pairs = list(combinations(r,2))
    #     for (x,y) in pairs:
    #system_a = analyse_system(systems[x])
    #system_b = analyse_system(systems[y])

def analyse_system(dna_system):
    lower_toeh = regex_match(dna_system,"lower_th")
    lower_toeh_c = regex_match(dna_system,"lower_th_c")
    upper_toeh = regex_match(dna_system,"upper_th")
    upper_toeh_c = regex_match(dna_system,"upper_th_c")
    return (DNASystem(dna_system, upper_toeh, lower_toeh, upper_toeh_c, lower_toeh_c))

dna = "{L' A^* R'}{L' B^* R'} {L' C^ R'} | <L D^ R> | <L E^* R> "
#simulate_strand(dna)

#print(s1)
#trajectory = process.sample({'A': 100}, steps=1000)
#for dt, transitions in trajectory:
#    print(trajectory.time, trajectory.state['A'], trajectory.state['A2'])

#Call everything

#upper_strand = re.compile('<[^{}\[\]]*?>')
#lower_strand = re.compile('{[^<>\[\]]*?}')
#toehold = re.compile('(\w)(?=\^\s)')
#toehold_c = re.compile('(\w)(?=\^\*)')

#upper_th = re.compile('<(?:[^{}\[\]]*?(\w)(?=(?:\^\s|\^>))+?)+?[^{}\[\]]*?>')
#lower_th = re.compile('{[^<>\[\]]*?\s(\w)\^\s+?[^<>\[\]]*?}')
#lower_th_c = re.compile('{[^<>\[\]]*?\s(\w)\^\*+?[^<>\[\]]*?}')
#upper_th = re.compile('<[^{}\[\]]*?(\w?)\^\s+?[^{}\[\]]*?>')
#upper_th_c = re.compile('<[^{}\]\[]*?\s+?(\w)\^\*+?[^{}\[\]]*?>')

#def regex_match(dna, position, status):
#     matches = []
#     print("LENGTH", len(list(re.finditer(position, dna))))
#     if(len(list(re.finditer(position, dna))))>0:
#         for strand in re.finditer(position, dna):
#             print(strand.string)
#             matches.append(re.finditer(status, strand.string))
#     return(matches)

    # def upper_toehold_binding(self, k, l):
    #     print("SWITCH")
    #     upper_toehold = re.finditer(upper_th, k)
    #     print("upper toehold", list(upper_toehold))
    #     lower_toehold_c = re.finditer(lower_th_c, l)
    #     print("lower toehold", list(lower_toehold_c))
    #     for upper_toehold in re.finditer(upper_th, k):
    #         print("upper %s: %s" % (upper_toehold.start(), upper_toehold.group()))
    #         for lower_toehold_c in re.finditer(lower_th_c, l):
    #             print("lower %s: %s" % (lower_toehold_c.start(), lower_toehold_c.group()))
    #             if upper_toehold.group() == lower_toehold_c.group():
    #                 part_A = k[:upper_toehold.start()] + ">"
    #                 print(upper_toehold.start())
    #                 print(part_A)
    #                 part_B = "<" + k[upper_toehold.end()+2:]
    #                 print(part_B)
    #                 part_C = l[:lower_toehold_c.start()] + "}"
    #                 part_D = "{" + l[lower_toehold_c.end()+3:]
    #                 part_E = "[" + upper_toehold.group() + "^]"
    #                 print(part_C)
    #                 print(part_D)
    #                 final = part_A + part_C + part_E + part_D + part_B
    #                 last = re.sub(empty_bracket, '', final)
    #                 print(last)
    #                 yield self.Transition([k, l], [final], alpha)

    # def generic_toehold_binding(self, k, l, regex_1, regex_2):
    #     for matching_1 in re.finditer(regex_1, k):
    #         for matching_2 in re.finditer(regex_2, l):
    #             if matching_1.group() == matching_2.group():
    #                 if regex_1 == upper_th or regex_1 == upper_th_c:
    #                     print("upper_th")
    #                     part_A = k[:matching_1.start()] + ">" + l[:matching_2.start()] + "}"
    #                     part_B ="{" + l[matching_2.end()+2:] + "<" + k[matching_1.end()+2:]
    #                 #elif regex_1 == upper_th_c:
    #                  #   part_A = k[:matching_1.start()] + ">" + l[:matching_2.start()] + "}"
    #                   #  part_B = "{" + l[matching_2.end()+2:] + "<" + k[matching_1.end()+2:]
    #                 elif regex_1 == lower_th or regex_1 == lower_th_c:
    #                     part_A = k[:matching_1.start()] + "}" + l[:matching_2.start()] + ">"
    #                     part_B = "<" + l[matching_2.end()+3:] + "{" + k[matching_1.end()+2:]
    #                 #elif regex_1 == lower_th_c:
    #                  #   print("lower_th_c")
    #                     #part_A = k[:matching_1.start()] + "}" + l[:matching_2.start()] + ">"
    #                     #part_B = "<" + l[matching_2.end()+2:] + "{" + k[matching_1.end()+2:]
    #                 draft_strand = part_A + "[" + matching_2.group() + "^]" + part_B
    #                 final_strand = re.sub(empty_bracket, '', draft_strand)
    #                 print(final_strand)
    #                 yield self.Transition([k, l], [final_strand], alpha)
