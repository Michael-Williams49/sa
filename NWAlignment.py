# Represents a sequence with a name and a sequence string
class Sequence:
    def __init__(self, name: str, seq: str):
        """
        Initialize a Sequence object with a name and a sequence string.

        Args:
            name (str): The name of the sequence.
            seq (str): The sequence string.
        """
        self.name = name
        self.seq = seq

    def __repr__(self):
        return f"{self.name}: {self.seq}"

# Represents a table for storing scores and alignments
class Table:
    def __init__(self, seqX: Sequence, seqY: Sequence):
        """
        Initialize a Table object with two sequences.

        Args:
            seqX (Sequence): The first sequence.
            seqY (Sequence): The second sequence.
        """
        self.seqX = "-" + seqX.seq  # Add a gap character at the beginning
        self.seqY = "-" + seqY.seq  # Add a gap character at the beginning
        self.nameX = seqX.name
        self.nameY = seqY.name
        self.data = [[0 for _ in range(len(self.seqY))] for _ in range(len(self.seqX))]

    def __repr__(self):
        """
        Return a string representation of the table.

        Returns:
            str: A string representation of the table.
        """
        result = ""
        result += " \t"  # Add a tab for the top-left corner
        for i in range(len(self.seqY)):
            result += self.seqY[i] + "\t"  # Add the sequence Y characters
        result += "\n"
        for row in range(len(self.seqX)):
            result += self.seqX[row] + "\t"  # Add the sequence X character
            for col in range(len(self.seqY)):
                result += str(self.data[row][col]) + "\t"  # Add the score or alignment value
            result += "\n"
        result += "\n"
        return result

# Represents a scoring scheme for matches, mismatches, and indels
class Scheme:
    def __init__(self, match: int, mismatch: int, indel: int):
        """
        Initialize a Scheme object with scores for match, mismatch, and indel.

        Args:
            match (int): The score for a match.
            mismatch (int): The score for a mismatch.
            indel (int): The score for an indel (insertion or deletion).
        """
        self.match = match
        self.mismatch = mismatch
        self.indel = indel

# Implements the Needleman-Wunsch algorithm for sequence alignment
class NWA:
    def __init__(self, seqX: Sequence, seqY: Sequence, scoreScheme: Scheme = Scheme(1, -1, -1), preferMismatch: bool = True, allAlignment: bool = False):
        """
        Initialize an NWA object with two sequences and optional parameters.

        Args:
            seqX (Sequence): The first sequence.
            seqY (Sequence): The second sequence.
            scoreScheme (Scheme, optional): The scoring scheme for matches, mismatches, and indels. Defaults to Scheme(1, -1, -1).
            preferMismatch (bool, optional): Whether to prefer mismatches over indels in case of a tie. Defaults to True.
            allAlignment (bool, optional): Whether to find all possible alignments or only the optimal one. Defaults to False.
        """
        self.score = Table(seqX, seqY)
        self.align = Table(seqX, seqY)
        self.scheme = scoreScheme
        if preferMismatch:
            self.pref = [2, 1, 3]  # Preference order: match, indel in X, indel in Y
        else:
            self.pref = [1, 3, 2]  # Preference order: indel in X, indel in Y, match
        self.all = allAlignment
        self.alignment = []

        self.__initialize()
        self.__control()
        stack = []
        self.__trace(len(self.score.seqX) - 1, len(self.score.seqY) - 1, stack)

    def __initialize(self):
        """
        Initialize the score and alignment tables with base cases for indels.
        """
        for i in range(1, len(self.score.seqY)):
            self.score.data[0][i] += self.scheme.indel
            self.align.data[0][i] = 3  # Indel in Y
        for i in range(1, len(self.score.seqX)):
            self.score.data[i][0] += self.scheme.indel
            self.align.data[i][0] = 1  # Indel in X

    def __evaluate(self, i: int, j: int):
        """
        Evaluate the score and alignment for the given indices based on the scoring scheme and preference order.

        Args:
            i (int): The row index (corresponding to seqX).
            j (int): The column index (corresponding to seqY).
        """
        options = [0, 0, 0, 0]
        options[1] = self.score.data[i - 1][j] + self.scheme.indel  # Indel in X
        options[2] = self.score.data[i - 1][j - 1]  # Match or mismatch
        if self.score.seqX[i] == self.score.seqY[j]:
            options[2] += self.scheme.match
        else:
            options[2] += self.scheme.mismatch
        options[3] = self.score.data[i][j - 1] + self.scheme.indel  # Indel in Y
        bestOption = max(options[1:])
        for index in self.pref:
            if bestOption == options[index]:
                self.score.data[i][j] = options[index]
                self.align.data[i][j] = index + self.align.data[i][j] * 4  # Encode the preference order
                if not self.all:
                    break

    def __control(self):
        """
        Control the evaluation of scores and alignments for all positions in the table.
        """
        for sum in range(2, len(self.score.seqX) + len(self.score.seqY) - 1):
            for row in range(max(1, sum - len(self.score.seqY) + 1), min(sum, len(self.score.seqX))):
                self.__evaluate(row, sum - row)
                # print(self.score)  # Uncomment to print the score table after each evaluation

    def __trace(self, i: int, j: int, stack: list):
        """
        Recursively trace back the optimal alignment path and store it in the stack.

        Args:
            i (int): The row index (corresponding to seqX).
            j (int): The column index (corresponding to seqY).
            stack (list): The list to store the alignment path.
        """
        if i == 0 and j == 0:
            self.__result(stack)
            return
        direction = self.align.data[i][j]
        while direction != 0:
            if direction % 4 == 1:
                self.__trace(i - 1, j, stack + [1])  # Indel in X
            elif direction % 4 == 2:
                self.__trace(i - 1, j - 1, stack + [2])  # Match or mismatch
            elif direction % 4 == 3:
                self.__trace(i, j - 1, stack + [3])  # Indel in Y
            direction //= 4

    def __result(self, stack: list):
        """
        Construct the final alignment sequences from the traced alignment path.

        Args:
            stack (list): The list containing the alignment path.
        """
        seqX = ""
        seqY = ""
        indexX = 0
        indexY = 0
        for index in range(len(stack) - 1, -1, -1):
            if stack[index] == 1:  # Indel in X
                indexX += 1
                seqX += self.score.seqX[indexX]
                seqY += "-"
            elif stack[index] == 2:  # Match or mismatch
                indexX += 1
                indexY += 1
                seqX += self.score.seqX[indexX]
                seqY += self.score.seqY[indexY]
            elif stack[index] == 3:  # Indel in Y
                indexY += 1
                seqX += "-"
                seqY += self.score.seqY[indexY]
        self.alignment.append([Sequence(self.score.nameX, seqX), Sequence(self.score.nameY, seqY)])

if __name__ == "__main__":
    """
    This is the main entry point of the program.
    It creates an example NWA object and prints the resulting alignment(s).
    """
    example = NWA(Sequence("1", "abcdccdc"), Sequence("2", "abbcadad"))
    print(example.alignment)
    """
    Output:
    [Sequence('1', 'abcd-ccdc'), Sequence('2', 'abbc-adad')]
    """