from aligner import NWA, Sequence, Scheme
import numpy as np
import os

def save_to_file(data, filename):
    """Save data to a text file."""
    with open(filename, 'w') as f:
        for row in data:
            f.write(' '.join(map(str, row)) + '\n')
    print(f"Data saved to {filename} successfully.")

def load_from_file(filename, shape):
    """Load data from a text file."""
    data = np.zeros(shape, dtype=int)
    with open(filename, 'r') as f:
        for i, line in enumerate(f.readlines()):
            data[i] = np.array(list(map(int, line.split())))
    print("Data loaded from file:")
    print(data)
    return data

def compute_alignment_matrix(seqA, seqB, match_score, mismatch_penalty, gap_penalty):
    """Compute the alignment matrix using dynamic programming."""
    m, n = len(seqA), len(seqB)
    matrix = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(1, m + 1):
        matrix[i][0] = i * gap_penalty
    for j in range(1, n + 1):
        matrix[0][j] = j * gap_penalty
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seqA[i - 1] == seqB[j - 1]:
                score_diagonal = matrix[i - 1][j - 1] + match_score
            else:
                score_diagonal = matrix[i - 1][j - 1] + mismatch_penalty
            score_up = matrix[i - 1][j] + gap_penalty
            score_left = matrix[i][j - 1] + gap_penalty
            matrix[i][j] = max(score_diagonal, score_up, score_left)
    return matrix

def process_blocks(seqA, seqB, block_size):
    match_score = 1
    mismatch_penalty = -1
    gap_penalty = -1

    global_alignment_data = compute_alignment_matrix(seqA, seqB, match_score, mismatch_penalty, gap_penalty)

    save_to_file(global_alignment_data, 'global_alignment_data.txt')
    loaded_data = load_from_file('global_alignment_data.txt', (len(seqA) + 1, len(seqB) + 1))

    alignmentA, alignmentB = "", ""
    i, j = len(seqA), len(seqB)
    print("Starting traceback...")
    while i > 0 or j > 0:
        print(f"Traceback at i={i}, j={j}, current score={loaded_data[i][j]}")
        if i > 0 and j > 0 and loaded_data[i][j] == loaded_data[i-1][j-1] + (match_score if seqA[i-1] == seqB[j-1] else mismatch_penalty):
            alignmentA = seqA[i-1] + alignmentA
            alignmentB = seqB[j-1] + alignmentB
            i -= 1
            j -= 1
        elif i > 0 and loaded_data[i][j] == loaded_data[i-1][j] + gap_penalty:
            alignmentA = seqA[i-1] + alignmentA
            alignmentB = '-' + alignmentB
            i -= 1
        elif j > 0 and loaded_data[i][j] == loaded_data[i][j-1] + gap_penalty:
            alignmentA = '-' + alignmentA
            alignmentB = seqB[j-1] + alignmentB
            j -= 1
        else:
            print("No valid move found, breaking out of loop.")
            break

    return alignmentA, alignmentB

def test_process_blocks():
    seqA = "abcdccdc"
    seqB = "abbcadad"
    block_size = 3
    alignmentA, alignmentB = process_blocks(seqA, seqB, block_size)
    print("Alignment A:", alignmentA)
    print("Alignment B:", alignmentB)

if __name__ == '__main__':
    test_process_blocks()
