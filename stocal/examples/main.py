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