from aligner import NWA, Sequence, Scheme


def process_blocks(seqA, seqB, block_size):
    scheme = Scheme(1, -1, -1)
    m = len(seqA)
    n = len(seqB)
    # Initialize a list for each block column
    former_end_result_x = [[]
                           for _ in range((n + block_size - 1) // block_size)]

    # Prepare a global matrix to store scores
    global_alignment_data = [[0] * (n + 1) for _ in range(m + 1)]

    for block_row in range(0, m, block_size):
        former_end_result_y = None  # To hold vertical alignment data from the previous block

        for block_col in range(0, n, block_size):
            block_index = block_col // block_size
            endX = min(block_row + block_size, m)
            endY = min(block_col + block_size, n)
            seqX = Sequence("seqX", seqA[block_row:endX])
            seqY = Sequence("seqY", seqB[block_col:endY])

            initial_col_data = former_end_result_y if block_col > 0 else None
            initial_row_data = former_end_result_x[
                block_index][-1] if block_row > 0 and former_end_result_x[block_index] else None

            block_aligner = NWA(seqX, seqY, Scheme(1, -1, -1),
                                initial_col_data=initial_col_data,
                                initial_row_data=initial_row_data,
                                block_start_row=block_row,  # Use actual global starting position
                                block_start_col=block_col)  # Use actual global starting position

        # Store the last row and column for the next block initialization
            former_end_result_y = [
                row[-1]for row in block_aligner.score.data]

            former_end_result_x[block_index].append(
                [block_aligner.score.data[-1][col] for col in range(len(block_aligner.score.data[0]))])

            # Update the global alignment data matrix with the block results
            for i in range(len(seqX.seq) + 1):
                for j in range(len(seqY.seq) + 1):
                    global_row = block_row + i
                    global_col = block_col + j
                    if global_row <= m and global_col <= n:
                        global_alignment_data[global_row][global_col] = block_aligner.score.data[i][j]

    # reconstruct the global alignment from global_alignment_data
    alignmentA, alignmentB = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and global_alignment_data[i][j] == global_alignment_data[i-1][j-1] + (scheme.match if seqA[i-1] == seqB[j-1] else scheme.mismatch):
            alignmentA = seqA[i-1] + alignmentA
            alignmentB = seqB[j-1] + alignmentB
            i -= 1
            j -= 1
        elif i > 0 and global_alignment_data[i][j] == global_alignment_data[i-1][j] + scheme.indel:
            alignmentA = seqA[i-1] + alignmentA
            alignmentB = '-' + alignmentB
            i -= 1
        else:
            alignmentA = '-' + alignmentA
            alignmentB = seqB[j-1] + alignmentB
            j -= 1
    for row in global_alignment_data:
        print("\t".join(map(str, row)))


    return alignmentA, alignmentB

        
    


# test block
def test_process_blocks():
    seqA = "abcdccdc"
    seqB = "abbcadad"
    block_size = 3

    process_blocks(seqA, seqB, block_size)

    alignmentA, alignmentB = process_blocks(seqA, seqB, 3)
    print("Alignment A:", alignmentA)
    print("Alignment B:", alignmentB)


if __name__ == '__main__':
    test_process_blocks()
