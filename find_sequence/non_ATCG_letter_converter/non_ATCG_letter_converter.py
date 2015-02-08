#-------------------------------------------------------------------------------
# Name:        non_ATCG_letter_converter
# Purpose:      To match ambiguous DNA base pairs with their corresponding
#               boolean values.
#
# Author:      Julian
#
# Created:     19/08/2014
# Copyright:   (c) Julian 2014
# Licence:     <pending>
#-------------------------------------------------------------------------------

# non_ATCG_letter_converter() takes non-ATCG base pair as string and next
# correct base pair from recognition sequence as string and returns a boolean
# to determine whether to continue pairing for restriction enzyme, True if
# match is allowed, False otherwise


# Conversions based on following table
'''
Nucleotide Symbols Used:
======================================================================
R = A or G       M = A or C       H = A, C or T       N = A, C, G or T
Y = C or T       K = G or T       V = A, C or G
                 S = C or G       B = C, G or T
                 W = A or T       D = A, G or T
======================================================================
'''

# Tsil of valid bp's used in test_for_bp_symbols
valid_bp_tsil = ['A', 'T', 'C', 'G', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'V', \
'B', 'D', 'N']

def non_ATCG_letter_converter(bp_read, bp_recseq):
    if test_for_ATCG(bp_read) == True and test_for_ATCG(bp_recseq) == True:
        raise TypeError('Either base pair must be non-ATCG.')
    elif len(bp_read) > 1 or len(bp_recseq) > 1:
        raise TypeError('Only one base pair may be compared at a time.')
    elif test_for_bp_symbols(bp_recseq) == False or \
    test_for_bp_symbols(bp_read) == False:
        raise TypeError('Invalid symbol for bp - ' + bp_read + ' or ' +\
    bp_recseq + '.')
    elif bp_recseq == bp_read:
        return True
    elif bp_recseq == 'N':
        return True
    elif bp_recseq == 'R':
        if bp_read == 'A' or bp_read == 'G':
            return True
        else:
            return False
    elif bp_recseq == 'Y':
        if bp_read == 'C' or bp_read == 'T':
            return True
        else:
            return False
    elif bp_recseq == 'M':
        if bp_read == 'A' or bp_read == 'C':
            return True
        else:
            return False
    elif bp_recseq == 'K':
        if bp_read == 'G' or bp_read == 'T':
            return True
        else:
            return False
    elif bp_recseq == 'S':
        if bp_read == 'C' or bp_read == 'G':
            return True
        else:
            return False
    elif bp_recseq == 'W':
        if bp_read == 'A' or bp_read == 'T':
            return True
        else:
            return False
    elif bp_recseq == 'H':
        if bp_read == 'A' or bp_read == 'C' or bp_read == 'T' or \
        bp_read == 'M' or bp_read == 'W' or bp_read == 'Y':
            return True
        else:
            return False
    elif bp_recseq == 'V':
        if bp_read == 'A' or bp_read == 'C' or bp_read == 'G' or \
        bp_read == 'R' or bp_read == 'M' or bp_read == 'S':
            return True
        else:
            return False
    elif bp_recseq == 'B':
        if bp_read == 'C' or bp_read == 'G' or bp_read == 'T' or \
        bp_read == 'Y' or bp_read == 'K' or bp_read == 'S':
            return True
        else:
            return False
    elif bp_recseq == 'D':
        if bp_read == 'A' or bp_read == 'G' or bp_read == 'T' or \
        bp_read == 'R' or bp_read == 'K' or bp_read == 'W':
            return True
        else:
            return False
    elif test_for_ATCG(bp_recseq) == True:
        return False
    else:
        raise TypeError('nonATCGerror - Base pair is of invalid \
        letterhead.')


# Utility method for non_ATCG_to_letter_converter used for taking in a bp and
# returning True if that bp is A or T or C or G or False otherwise
# If kwarg return_non_ATCG_bp is not False, then bp is returned instead of True
def test_for_ATCG(bp_to_test, return_non_ATCG_bp = False):
    if bp_to_test == 'A' or bp_to_test == 'T' or bp_to_test == 'C' or \
    bp_to_test == 'G':
        return True
    else:
        if return_non_ATCG_bp == False:
            return False
        else:
            return bp_to_test

# Utility method for non_ATCG_to_letter_converter used to test whether bp args
# are valid nucleotide symbols, tests each bp in arg string and returns False
# if one of them is invalid
def test_for_bp_symbols(bp_to_test):
    result = False
    for poss_bp in valid_bp_tsil:
        if bp_to_test == poss_bp:
            result = True
    return result
