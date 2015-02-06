#-------------------------------------------------------------------------------
# Name:        GAT - Sequence similarity_from to
# Purpose:
#
# Author:      James
#
# Created:     02/02/2015
# Copyright:   (c) James 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------


"""
    Module meant to take in two or more whole genome strings and compare them
    to another whole genome string.  Whole genome string arguments are
    expected to be of type string, in all uppercase format, and to start at
    the [5' or 3'] end.  Module will be able to take in optional arguments to
    specify the starting position for comparing every genome.  Module will
    return the percentage of sequence similarity between the target genome
    for comparison and all other genomes.  Module will assume arguments are in
    the correct format.
"""


#  Remember *args returns a tuple of values, here args should be a tuple of
#  arrays, with each array consisting of a genome sequence in string type and
#  an optional second integer value to determine starting position for
#  searching the sequence.
def main(wg_to_compare, *args, startpos_wg_to_compare=0, printout=False):
    matches_dict = {}

    # iter over each genome other than comparison genome.
    for x in range(len(args)): ### x is an integer that corresponds to each
        if len(args[x]) == 2:  ### array of another genome
            startpos_this_other_wg = args[x][1]
        else:
            startpos_this_other_wg = 0

        # check if other genome is longer than genome to compare, if so capture
        # for control of conditional statement below.
        if len(wg_to_compare) < len(args[x][0]):
            length_comparison = "other_genome_greater"
        else:
            length_comparison = "genome_to_compare_greater"
        length_comparison_dict = {"length_value" : length_comparison}
        print (length_comparison_dict)


        num_of_matches = 0
        tot_nucleotides = 0

        # iter over each nucleotide in whole genome to compare and other genome.
        for y in range(startpos_wg_to_compare, len(wg_to_compare)):

            if length_comparison_dict["length_value"] == \
            "genome_to_compare_greater" and y == \
            len(args[x][0][startpos_this_other_wg:len(args[x][0])]):
                break

            tot_nucleotides += 1

            if wg_to_compare[y] == args[x][0][y+startpos_this_other_wg]:
                num_of_matches += 1

        matches_dict[x+1] = [num_of_matches, \
        (num_of_matches/tot_nucleotides)*100]

    if printout == True:
        print ("Results:")
        for z in range(len(args)):
            print (str(matches_dict[z+1][0])+" matches, "+\
            str(matches_dict[z+1][1])+"% similarity.")

    else:
        for z in range(len(args)):
            return matches_dict
