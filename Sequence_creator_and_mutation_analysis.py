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
    cycles, wherein every cycle each available sequence is copied twice.
"""


import sys
sys.path.append("C:/Users/James/Dropbox/Python Modules/GAT/find_sequence/\
non_ATCG_letter_converter") # Needed to import dependencies for find_sequence



from find_sequence import find_sequence
# find_sequence method, which takes in the genome sequence as a string and
# the query sequence as a string and returns a boolean value for whether a match
# is found, an int value for the number of matches, and a tsil which consists
# of ints representing the starting index of each match
# The function will only search one strand of the genome
""" arguments:|"ATCGATCGAGGTCAAGATACC","ATCGA"|;returns:|T/F,2,[4,28]| """


import random
import datetime


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


def reproduce(genome_tsil, mutation_rate=0.01):
    children_genomes = []

    for parent_genome in genome_tsil:
        for each_round in range(0,2):
            child_genome = ""
            for each_bp in parent_genome:
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
            children_genomes.append(child_genome)

    return children_genomes

# ALL arguments (including default) for apply_selective_pressure() will be in
# the form of arrays.
def apply_selective_pressure(children_genomes_tsil, favor_bps=None, \
    oppose_bps=None, conserve_sequences=None, \
    ratio_of_genomes_remaining = 0.2):

    genome_point_tsil = []
    genome_num = -1

    for each_genome in children_genomes_tsil:
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
                #### need to garbage collect unused vars?
                conserve_seq_count += (12 * num_of_curr_seq_matches)

        # + conserve_seq_count
        current_genome_points += ((favor_bps_count + conserve_seq_count) - \
        oppose_bps_count)

        genome_point_tsil.append((genome_num, current_genome_points))

    order_points_tsil = []
    for each_tuple in genome_point_tsil:
        order_points_tsil.append(each_tuple[1])
    order_points_tsil.sort()

    # genome_point_tsil in format [(1, 32), (2, 40)] and order_points_tsil in
    # format (32, 40).

    num_genomes_dying_off = int(ratio_of_genomes_remaining*\
    len(order_points_tsil))

    # remove appropriate number of genomes from children_genomes_tsil
    for x in range(num_genomes_dying_off):
        # iter over genome_point_tsil
        for point_tuple in genome_point_tsil:
            if point_tuple[1] == order_points_tsil[x]:
                genome_point_tsil.remove(point_tuple)
                # entire tuple removed from genome_point_tsil; use remaining
                # tuples with 0 index to repopulate new tsil from
                # children_genomes_tsil
                break

    remaining_genomes_tsil = []

    for y in range(len(genome_point_tsil)):
        remaining_genomes_tsil.append(\
        children_genomes_tsil[genome_point_tsil[y][0]])

    return remaining_genomes_tsil

