#-------------------------------------------------------------------------------
# Name:        Seq_creator_and_mut_anlaysis RUNTIME ANALYSIS
# Purpose:
#
# Author:      James
#
# Created:     07/02/2015
# Copyright:   (c) James 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------


import Sequence_creator_and_mutation_analysis

import datetime, sys
sys.path.append("C:/Users/James/Dropbox/Python Modules/GAT")

start_time = datetime.datetime.now()

def main():
    origin_of_all = Sequence_creator_and_mutation_analysis.sequence_creator()
    members = Sequence_creator_and_mutation_analysis.reproduce([origin_of_all])
    for x in range(3):
        members = Sequence_creator_and_mutation_analysis.reproduce(members)
        print ("early repro complete")

    for y in range(5):
        members = Sequence_creator_and_mutation_analysis.reproduce(members)
        print ("repro complete")
        members = Sequence_creator_and_mutation_analysis.\
        apply_selective_pressure(members, favor_bps = ["A"])
        print ("selection complete")

if __name__ == '__main__':
    main()

end_time = datetime.datetime.now()

print (end_time-start_time)