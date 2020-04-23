from stocal import *
from Rules import *
import re
import itertools
from itertools import combinations

# Two reactions: Two molecules of A forming a dimer A2 and reverse
r1 = MassAction({'A': 2}, {'A2': 1}, 1.0)
r2 = MassAction({'A2': 1}, {'A': 2}, 10.0)  # Last value is the stochastic rate at which the reaction occurs
# process = Process([r1, r2])

# Here, the event feed will occur, and feed an A molecule to the system, at time t=0.0 and then periodically every 1.0 time units.
# Unlike stochastic reactions that occur with an average frequency, nondeterministic events happen at exactly the specified times.
feed = Event([], ['A'], 0.0, 1.0)

process = Process(transitions=[r1, r2], rules=[Dilution()])

def simulate_strand(str):
    separate = str.split("|")
    if (len(separate)>1):
        for i in range(len(separate)):
            print(i)
            lower_toeh = re.findall('{[^<>\[\]]*?\s(\w)\^\s+?[^<>\[\]]*?}', separate[i])
            print(lower_toeh)
            lower_toeh_c = re.findall('{[^<>\[\]]*?\s(\w)\^\*+?[^<>\[\]]*?}', separate[i])
            print(lower_toeh_c)
            upper_toeh = re.findall('<[^{}\[\]]*?(\w?)\^\s+?[^{}\[\]]*?>', separate[i])
            print(upper_toeh)
            upper_toeh_c = re.findall('<[^{}\]\[\]]*?\s+?(\w)\^\*+?[^{}\[\]]*?>', separate[i])
            print(upper_toeh_c)



        #r = list(range(0,len(separate)))
        #pairs = list(combinations(r,2))
        #for (x,y) in pairs:
            #print(x)
            #print(y)

        #for x in len(separate):
        #    p = re.compile('{.*?\s(\w)\^\*?.*?}'
        #    '({.*?\w\^\*?.*?})', separate[0])
        #p = re.compile('ab*')
    #print(separate)

#Test Section
#s1 = MassAction({"{L' N^* R'}": 1, "<L N^ R>": 1}, {"{L' }<L>[N^]<R>{R'}": 1}, 1.0)
#((.)*(\{)+(.)*(\w\^\*)+(.)*(\})+(.)*(/|)+(.)*(<)+(.)*(\w\^)+(.)*(>)+(.)*)|((.)*(<)+(.)*(\w\^)+(.)*(>)+(.)*(/|)+(.)*(/{)+(.)*(\w\^\*)+(.)*(/})+(.)*)|((.)*(\{)+(.)*(\w\^)+(.)*(\})+(.)*(/|)+(.)*(<)+(.)*(\w\^\*)+(.)*(>)+(.)*)|((.)*(<)+(.)*(\w\^\*)+(.)*(>)+(.)*(/|)+(.)*(/{)+(.)*(\w\^)+(.)*(/})+(.)*)
#process = Process(transitions=[r1, r2], rules=[])
#(.*\{+.*(\w\^\*)+.*\}+.*(/|)+.*<+.*(\w\^)+.*>+.*)|(.*<+.*(\w\^)+.*>+.*(/|)+.*(/{)+.*(\w\^\*)+.*(/})+.*)|(.*\{+.*(\w\^)+.*\}+.*(/|)+.*<+.*(\w\^\*)+.*>+.*)|(.*<+.*(\w\^\*)+.*>+.*(/|)+.*(/{)+.*(\w\^)+.*(/})+.*)

str = "{L' A^* R'}{L' B^* R'} {L' C^ R'} | <L D^ R> | <L E^* R> "
simulate_strand(str)

#print(s1)
#trajectory = process.sample({'A': 100}, steps=1000)
#for dt, transitions in trajectory:
#    print(trajectory.time, trajectory.state['A'], trajectory.state['A2'])

#Call everything