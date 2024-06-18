from genData.utils.statistics import compute_per_sequence_nucleotide_content, compute_per_sequence_dinucleotide_content

def test_compute_per_sequence_nucleotide_content():
    sequences = ['ATCG', 'GCTA', 'CGAT', 'AAAG']
    expected_result = {
        0: {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25},
        1: {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25},
        2: {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25},
        3: {'A': 0.75, 'G': 0.25}
    }
    result = compute_per_sequence_nucleotide_content(sequences)
    assert result == expected_result, "Test failed: compute_per_sequence_nucleotide_content"

from genData.utils.statistics import compute_per_sequence_nucleotide_content

def test_compute_per_sequence_dinucleotide_content():
    sequences = ['ATATATA', 'AAAAAA']

    expected_result = {
        0: {'AT': 0.5, 'TA': 0.5}, 
        1: {'AA': 1.0},
    }
    result = compute_per_sequence_dinucleotide_content(sequences)
    assert result == expected_result, "Test failed: compute_per_sequence_dinucleotide_content"