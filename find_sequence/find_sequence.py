#-------------------------------------------------------------------------------
# Name:        find_sequence
# Purpose:      The purpose of this module is to find the number of matches
#               between a desired sequence and another sequence as well as
#               their locations.
#
# Author:      jkougl
#
# Created:     23/08/2014
# Copyright:   (c) jkougl 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------


from non_ATCG_letter_converter import test_for_ATCG
from non_ATCG_letter_converter import non_ATCG_letter_converter



# find_sequence method, which takes in the genome sequence as a string and
# the query sequence as a string and returns a boolean value if a match is
# found, an int value for the number of matches, and a tsil which consists
# of ints representing the starting index of each match
# The function will only search one strand of the genome
def find_sequence(seq_to_search, recognition_seq, overlapping_matches = False):
    current_genome_bp_index = -1
    poss_start_of_match = False
    is_a_match = False
    matched_recently = False

    tsil_of_matching_indexes = []

    if discriminate_find_seq_args(seq_to_search, recognition_seq) == False:
        pass
    else:
        raise TypeError('Error in find_sequence() utility method.')

    num_of_matches_needed = len(recognition_seq)

    for genome_bp in seq_to_search: # Start genome search, only
        non_ATCG_for_genome = False # return here after complete match
        current_genome_bp_index += 1 # is successful or a non-matching bp
        num_of_matches_found = 0
        bp_match_found = False
        break_again = False
        poss_start_of_match = current_genome_bp_index

        if matched_recently == True:
            matched_recently = False

            if overlapping_matches == False:

                if current_genome_bp_index == matched_final_bp_index + 1:
                    allow_matching = True # Allow bp_matching to start
                    #print 'match allowed with ' + str(matched_final_bp_index+1)
                else:
                     allow_matching = False # Move on to genome_bp after match
                     matched_recently = True # Enter current if loop on next iter
                     #print 'match denied'

            else:
                allow_matching = True # Allow bp_matching without the check

        else:
            allow_matching = True

        #print 'current genome bp is ', genome_bp, ' at ', str(current_genome_bp_index)

        if allow_matching == True:

            # Check if match already found, if so update variables
            for each_match in range(num_of_matches_needed + 1):

                if bp_match_found == True:
                    num_of_matches_found += 1
                    bp_match_found = False
                elif break_again == True:
                    break

                #print 'matches needed is ', num_of_matches_needed
                #print 'matches found is ', num_of_matches_found
                # Has a full match been achieved, if so leave loop
                if num_of_matches_needed == num_of_matches_found:
                    tsil_of_matching_indexes.append(poss_start_of_match)
                    is_a_match = True
                    matched_recently = True
                    matched_final_bp_index = current_genome_bp_index +\
                    num_of_matches_found - 1
                    #print 'tsil appended'
                    break

                #print 'current match is ', num_of_matches_found


                for bp_to_match in recognition_seq:
                    # Start matching process
                    non_ATCG_for_recognition_seq = False
                    non_ATCG_for_genome = False

                    # Check if current index is past end of string
                    if current_genome_bp_index + num_of_matches_found >\
                     len(seq_to_search) - 1:
                        #print 'reached end of seq_to_search'
                        break
                    elif num_of_matches_found > len(recognition_seq) - 1:
                        #print 'reached end of query seq'
                        break

                    #print 'genome bp is '+seq_to_search[current_genome_bp_index+\
                    #num_of_matches_found]+' and recognition bp is '+\
                    #recognition_seq[num_of_matches_found]

                    if test_for_ATCG(seq_to_search[current_genome_bp_index +\
                     num_of_matches_found]) == False: # Is genome_bp ATCG?
                        non_ATCG_for_genome = True # If true note it and get bp
                        non_ATCG_genome_bp = test_for_ATCG(seq_to_search\
                        [current_genome_bp_index + num_of_matches_found], \
                        return_non_ATCG_bp = True)
                    if test_for_ATCG(recognition_seq[num_of_matches_found])\
                     == False: # Is bp_to_match ATCG?
                        non_ATCG_for_recognition_seq = True
                        non_ATCG_matching_bp = test_for_ATCG(\
                        recognition_seq[num_of_matches_found], \
                        return_non_ATCG_bp = True) # If true note it and get bp

                    if non_ATCG_for_genome == True or \
                    non_ATCG_for_recognition_seq == True: # Are any bp's non_ATCG?
                        if non_ATCG_letter_converter(seq_to_search\
                        [current_genome_bp_index + num_of_matches_found],\
                         recognition_seq[num_of_matches_found]) == True:
                            bp_match_found = True
                            break # bp's match, move to next corresponding bp's
                        else:
                            break_again = True
                            break # bp's do not match; move on to next genome_bp

                    else: # no non_ATCG bp's present
                        if seq_to_search[current_genome_bp_index +\
                         num_of_matches_found] == \
                         recognition_seq[num_of_matches_found]:
                            bp_match_found = True
                            break
                        else:
                            break_again = True
                            break # bp's do not match; move on to next genome_bp

        else:
            pass # do not match bp's, due to match already made, go to next bp


    return (is_a_match, len(tsil_of_matching_indexes),\
 tsil_of_matching_indexes)



# Utility method for find_sequence to ensure no empty string args or
# non-string args are accepted. arg1 is genome seq and arg2 is query seq.
def discriminate_find_seq_args(arg1, arg2):
    # Are args of type string? If not, raise appropriate error.
    if type(arg1) == str:
        non_string_arg1 = False
    else:
        non_string_arg1 = True
        raise TypeError('Genome sequence is of incorrect type ' +\
     str(type(arg1)))
    if type(arg2) == str:
        non_string_arg2 = False
    else:
        non_string_arg2 = True
        raise TypeError('Recognition sequence is of incorrect type ' +\
    str(type(arg2)))

    # Are string args empty?  If so, raise appropriate error.
    if len(arg1) < 1:
        empty_string_arg1 = True
        raise TypeError('Genome sequence is of insufficient length: ' +\
    str(len(arg1)))
    else:
        empty_string_arg1 = False
    if len(arg2) < 1:
        empty_string_arg2 = True
        raise TypeError('Recognition sequence is of insufficient length: ' +\
     str(len(arg2)))
    else:
        empty_string_arg2 = False

    return False
