from flask import Flask, render_template, request
from aligner import Sequence, NWA, SWA, Scheme
from LSA2 import process_blocks
import os

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        seq1_name = request.form['seq1_name']
        seq1_sequence = request.form['seq1_sequence']
        seq2_name = request.form['seq2_name']
        seq2_sequence = request.form['seq2_sequence']
        show_all_alignments = request.form.get('show_all_alignments', False)
        selected_algorithm = request.form['algorithm']
        file_path = request.form.get('file_path', '')
        try:
            match_score = int(request.form.get('match_score', 2))
            mismatch_score = int(request.form.get('mismatch_score', -1))
            indel_score = int(request.form.get('indel_score', -2))
            ruler = int(request.form.get('ruler', 10))
        except ValueError:
            match_score = 2
            mismatch_score = -1
            indel_score = -2
            ruler = 10

        # match_score = int(request.form.get('match_score', 2))
        # mismatch_score = int(request.form.get('mismatch_score', -1))
        # indel_score = int(request.form.get('indel_score', -2))
        # ruler = int(request.form.get('ruler', 10))

        seq1 = Sequence(seq1_name, seq1_sequence)
        seq2 = Sequence(seq2_name, seq2_sequence)

        score_scheme = Scheme(match_score, mismatch_score, indel_score)

        # check if the sequence is too long
        if len(seq1_sequence) > 1000 and len(seq2_sequence) > 1000:
            # Change to file output
            class EmptyAligner:
                def __init__(self):
                    self.alignment = [[], [], [], []]
            aligner = EmptyAligner()
            alignmentA, AlignmentB, stack = process_blocks(
                seq1_sequence, seq2_sequence, 100)

            alignment_result = format_alignment_result(
                alignmentA, AlignmentB, stack)
            default_path = 'uploads'
            filename = f"alignment.txt"
            if not file_path:
                file_path = os.path.join(default_path, filename)
            else:
                file_path = os.path.join(file_path, filename)
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            with open(file_path, 'w') as file:
                file.write(alignment_result)

        else:

            # Choose the algorithm based on user selection
            if selected_algorithm == 'nwa':
                aligner = NWA(seq1, seq2, score_scheme,
                              allAlignment=show_all_alignments, ruler=ruler)

            elif selected_algorithm == 'swa':
                aligner = SWA(seq1, seq2, score_scheme,
                              allAlignment=show_all_alignments, ruler=ruler)

        # score_scheme = Scheme(match_score, mismatch_score, indel_score)

        # if selected_algorithm == 'nwa':
        #     aligner = NWA(seq1, seq2, score_scheme, allAlignment=show_all_alignments, ruler=ruler)
        # elif selected_algorithm == 'swa':
        #     aligner = SWA(seq1, seq2, score_scheme, allAlignment=show_all_alignments, ruler=ruler)

        return render_template('index.html', seq1_name=seq1_name, seq1_sequence=seq1_sequence, seq2_name=seq2_name, seq2_sequence=seq2_sequence, alignments=aligner.alignment, selected_algorithm=selected_algorithm, show_all_alignments=show_all_alignments, match_score=match_score, mismatch_score=mismatch_score, indel_score=indel_score, ruler=ruler)

    return render_template('index.html', match_score=2, mismatch_score=-1, indel_score=-2, ruler=10)


def format_alignment_result(alignmentA, alignmentB, stack):
    # Example formatting, adjust according to your needs
    return f"Alignment A: {alignmentA}\nAlignment B: {alignmentB}\nStack Trace: {stack}"


if __name__ == '__main__':
    app.run(debug=True)
