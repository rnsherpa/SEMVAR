import numpy as np

from semvar.utils import char_to_num, contains_invalid_chars, no_valid_kmers, revdnacomp, get_spdi

def test_char_to_num():
    np.testing.assert_array_equal(char_to_num('ACGT'), np.array([0,1,2,3], dtype=np.int32))
    np.testing.assert_array_equal(char_to_num('ACGTN'), np.array([0,1,2,3,7], dtype=np.int32)) 

def test_contains_invalid_chars():
    assert contains_invalid_chars('ACGT') == False
    assert contains_invalid_chars('ACGTN') == True

def test_no_valid_kmers():
    assert no_valid_kmers('ACGT', 4) == False
    assert no_valid_kmers('ACGN', 4) == True
    assert no_valid_kmers('NACGTN', 4) == False
    assert no_valid_kmers('NACGNACGNACGN', 4) == True

def test_revdnacomp():
    assert revdnacomp('ACGTA') == 'TACGT'

def test_get_spdi(): # Needs to be updated with indel handling when implemented
    assert get_spdi({'chr1':'NC_000001.11'}, 'chr1', 100, 'A', 'T') == 'NC_000001.11:100:A:T'

def test_ref_mismatch():
    pass # Function needs pyfaidx object to test