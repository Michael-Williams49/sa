from aligner import NWA, Sequence, Scheme


def global_traceback(traceback_matrix, seqA, seqB):
    alignmentA, alignmentB = "", ""
    i, j = len(seqA), len(seqB)

    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 1:  # Indel in X
            alignmentA = seqA[i-1] + alignmentA
            alignmentB = '-' + alignmentB
            i -= 1
        elif traceback_matrix[i][j] == 2:  # Match/Mismatch
            alignmentA = seqA[i-1] + alignmentA
            alignmentB = seqB[j-1] + alignmentB
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 3:  # Indel in Y
            alignmentA = '-' + alignmentA
            alignmentB = seqB[j-1] + alignmentB
            j -= 1

    return alignmentA, alignmentB


def process_blocks(seqA, seqB, block_size):
    m = len(seqA)
    n = len(seqB)
    # Initialize a list for each block column and traceback matrix
    former_end_result_x = [[] for _ in range((n + block_size - 1) // block_size)]
    global_alignment_data = [[0] * (n + 1) for _ in range(m + 1)]
    traceback_matrix = [[0] * (n + 1) for _ in range(m + 1)]

    for block_row in range(0, m, block_size):
        former_end_result_y = None

        for block_col in range(0, n, block_size):
            block_index = block_col // block_size
            endX = min(block_row + block_size, m)
            endY = min(block_col + block_size, n)
            seqX = Sequence("seqX", seqA[block_row:endX])
            seqY = Sequence("seqY", seqB[block_col:endY])
            initial_col_data = former_end_result_y if block_col > 0 else None
            initial_row_data = former_end_result_x[block_index][-1] if block_row > 0 and former_end_result_x[block_index] else None

            block_aligner = NWA(seqX, seqY, Scheme(1, -1, -1), initial_col_data=initial_col_data, initial_row_data=initial_row_data, block_start_row=block_row, block_start_col=block_col)
            # Process and fill traceback
            for i in range(1, len(block_aligner.score.seqX)):
                for j in range(1, len(block_aligner.score.seqY)):
                    score_options = [
                        block_aligner.score.data[i-1][j] + block_aligner.scheme.indel,  # Indel in X
                        block_aligner.score.data[i][j-1] + block_aligner.scheme.indel,  # Indel in Y
                        block_aligner.score.data[i-1][j-1] + (block_aligner.scheme.match if seqA[block_row + i - 1] == seqB[block_col + j - 1] else block_aligner.scheme.mismatch)  # Match/Mismatch
                    ]
                    best_option = max(score_options)
                    block_aligner.score.data[i][j] = best_option
                    traceback_matrix[block_row + i][block_col + j] = score_options.index(best_option) + 1  # Adjust for 1-based indexing in traceback

            former_end_result_y = [row[-1] for row in block_aligner.score.data]
            former_end_result_x[block_index].append([block_aligner.score.data[-1][col] for col in range(len(block_aligner.score.data[0]))])

    return global_traceback(traceback_matrix, seqA, seqB)


def test_process_blocks():
    seqA = "ab"
    seqB = "ab"
    block_size = 3
    alignmentA, alignmentB = process_blocks(seqA, seqB, block_size)
    print("Alignment A:", alignmentA)
    print("Alignment B:", alignmentB)


test_process_blocks()
