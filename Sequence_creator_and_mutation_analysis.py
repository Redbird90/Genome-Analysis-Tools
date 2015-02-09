#-------------------------------------------------------------------------------
# Name:        GAT - Sequence creator and Mutation Analysis
# Purpose:
#
# Author:      James
#
# Created:     03/02/2015
# Copyright:   (c) James 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------


"""
    Module will create a sequence of specified length.  Module will then go
    through each base pair of created sequence and copy each base pair, applying
    the specified mutation rate to cause mismatches for each base pair.  Module
    will then repeat this process to create a total of two new sequences.
    Module will then undergo the copying sequence until a specified number of
    cycles, wherein every cycle each available sequence is copied twice.  Module
    will also be able to apply certain selective pressures specified by the user
    to reduce the number of genomes in an a list.
"""


import sys
sys.path.append("C:/Users/James/Dropbox/Python Modules/GAT/find_sequence/\
non_ATCG_letter_converter") # Needed to import dependencies for find_sequence



from find_sequence import find_sequence
### copied from original module:
# find_sequence method, which takes in the genome sequence as a string and
# the query sequence as a string and returns a boolean value for whether a match
# is found, an int value for the number of matches, and a tsil which consists
# of ints representing the starting index of each match
# The function will only search one strand of the genome
""" e.g. arguments:|"ATCGATCGAGGTCAAGATACC","ATCGA"|;returns:|T/F,2,[4,28]| """


import random
import datetime

"""
    sequence_creator takes in a default argument, length, to specify the
    desired length of the randomized genome to return.  The function will
    iter over length and use the import random to get an integer from 0-3,
    which will be used to decide between four possible base pairs to add
    to an empty string.  After iteration of length default argument is
    complete, the randomized genome sequence will be returned in string
    format.
"""
def sequence_creator(length = 100000):
    empty_string = ""

    for x in range(length):
        bp_key = random.randint(0,3)

        if bp_key == 0:
            empty_string = empty_string + "A"
        if bp_key == 1:
            empty_string = empty_string + "T"
        if bp_key == 2:
            empty_string = empty_string + "C"
        if bp_key == 3:
            empty_string = empty_string + "G"

    return empty_string

"""
    reproduce takes in a required list of genomes in string type and the
    default argument mutation_rate of type float.  The function iterates over
    the list of genomes (nested for A).  For every iteration, the function goes
    through two loops (nested for B) to copy the genome.  The function then
    iterates over every base pair in the parent genome (nested for C) and 
    tests to determine whether the copied base pair will be a mutation.  This 
    test is done by comparison of default argument mutation_rate to a random 
    float between 0 and 1.  If the mutant test fails, the original base pair 
    in the parent genome is added to an empty string.  If the mutation test 
    passes, a random integer between 0 and 3 determines the new base pair.
    After the last base pair is copied (end of nested for C), the new genome
    is appended to a list.  The same process is applied for the second copy of
    the parent genome (end of nested for B), and for all genomes in the
    list provided in the required argument(end of nested for C).  The copied
    genomes are then returned as string types in a list.  The function does
    NOT copy the original genome once and keep the parent genome.  Both copied
    genome sequences have the same potential for mutations, given by the
    default argument mutation_rate.  For mutations introduced into the copied
    genome, it should be noted that there is a 25% chance that the same base
    pair will be added into the copied genome, so no mutation would actually
    occur in the copied genome.
"""
def reproduce(genome_tsil, mutation_rate=0.01):
    children_genomes = []

    for parent_genome in genome_tsil: # A
        for each_round in range(0,2): # B
            child_genome = ""
            for each_bp in parent_genome: # C
                is_mutant = random.random()
                if is_mutant <= mutation_rate:
                    bp_key = random.randint(0,3)
                    if bp_key == 0:
                        mutant_bp = "A"
                    if bp_key == 1:
                        mutant_bp = "T"
                    if bp_key == 2:
                        mutant_bp = "C"
                    if bp_key == 3:
                        mutant_bp = "G"
                    new_bp = mutant_bp
                else:
                    new_bp = each_bp
                child_genome = child_genome + new_bp
            children_genomes.append(child_genome) # end B, C

    return children_genomes # end A

    
"""
    apply selective pressure function takes in a required list of genomes in
    string type and several default arguments.  These default args include
    favor_bps, which should consist of a list of singular base pairs in 
    string type, oppose_bps, which are in the same format as favor_bps,
    conserve_sequences, which should consist of a list of genome sequences of
    type string, and ratio_of_genomes_remaining, which consists of a float
    between 0-1.  The function iterates over each genome in the required
    argument (nested for A).  The function then counts the number of favored
    base pairs (if applicable) in the genome and adds 1 to a counter for each 
    favored base pair.  The function then counts the number of opposed base
    pairs (if applicable) in the genome and subtracts 1 from the counter.
    The function then uses the find_sequence method to count the number of
    matches (if applicable) for each conserved sequence and adds 12 for each
    match to the counter.  The counter value for that genome is then appended
    to a list of two-value tuples, with the first value being the number of
    the genome being processed and the second value being the counter value.
    After every genome has had a counter value attached to it, the counter
    list is then sorted from least to greatest.  The function then uses the
    default argument ratio_of_genomes_remaining to calculate how many genomes
    to remove.  The function then iterates over this number (nested for X)
    and iterates over each tuple in the list of tuples (nested for Y) and 
    determines whether that tuple has the lowest score based on the sorted 
    counter list.  When a match for the lowest score is found, that tuple is 
    removed from the list of tuples and the loop is broken (end of nested
    for Y).  After the appropriate number of tuples have been removed, the
    list of tuples is iterated over and the first value of each tuple is used
    to append the correct original genome to a new list.  The list of the 
    remaining original genomes is then returned in the same format as it was 
    first taken.
"""
# ALL arguments (including default) for apply_selective_pressure() will be in
# the form of arrays.
def apply_selective_pressure(children_genomes_tsil, favor_bps=None, \
    oppose_bps=None, conserve_sequences=None, \
    ratio_of_genomes_remaining = 0.2):

    genome_point_tsil = []
    genome_num = -1

    for each_genome in children_genomes_tsil: # A
        current_genome_points = 0
        genome_num += 1

        favor_bps_count = 0
        if favor_bps != None:
            for each_favored_bp in favor_bps:
                for bp in each_genome:
                    if bp == each_favored_bp:
                        favor_bps_count += 1

        oppose_bps_count = 0
        if oppose_bps != None:
            for each_opposed_bp in oppose_bps:
                for bp in each_genome:
                    if bp == each_opposed_bp:
                        oppose_bps_count += 1

        conserve_seq_count = 0
        if conserve_sequences != None:
            for each_conserved_seq in conserve_sequences:
                # find the function that finds the number of
                # occurrences of a desired sequence in a string.
                cons_seq_match, num_of_curr_seq_matches, seq_match_index_tsil =\
                find_sequence.find_sequence(each_genome, each_conserved_seq)
                ### need to garbage collect unused vars?
                conserve_seq_count += (12 * num_of_curr_seq_matches)

        # + conserve_seq_count
        current_genome_points += ((favor_bps_count + conserve_seq_count) - \
        oppose_bps_count)

        genome_point_tsil.append((genome_num, current_genome_points)) # end A

    order_points_tsil = []
    for each_tuple in genome_point_tsil:
        order_points_tsil.append(each_tuple[1])
    order_points_tsil.sort()

    # genome_point_tsil in format [(1, 32), (2, 40)] and 
    # order_points_tsil in format (32, 40).

    num_genomes_dying_off = int(ratio_of_genomes_remaining*\
    len(order_points_tsil))

    # remove appropriate number of genomes from children_genomes_tsil
    for x in range(num_genomes_dying_off): # X
        # iter over genome_point_tsil
        for point_tuple in genome_point_tsil: # Y
            if point_tuple[1] == order_points_tsil[x]:
                genome_point_tsil.remove(point_tuple)
                # entire tuple removed from genome_point_tsil; use remaining
                # tuples with 0 index to repopulate new tsil from
                # children_genomes_tsil
                break # end Y

    remaining_genomes_tsil = []

    for y in range(len(genome_point_tsil)):
        remaining_genomes_tsil.append(\
        children_genomes_tsil[genome_point_tsil[y][0]])

    return remaining_genomes_tsil

