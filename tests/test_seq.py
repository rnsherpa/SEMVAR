from semvar.seq import get_sequences, get_seq_context

def test_get_sequences():
    pass # Needs pyfaidx object to test

def test_get_seq_context():
    score_idx = 2
    chr = 'chr1'
    pos = 11100
    sem_len = 4
    seq_context = get_seq_context(score_idx, chr, pos, sem_len)
    assert seq_context == 'chr1:11098-11102'