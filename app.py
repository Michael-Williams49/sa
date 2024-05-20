from flask import Flask, render_template, request
from aligner import Sequence, NWA, SWA, Scheme

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

        match_score = int(request.form.get('match_score', 2))
        mismatch_score = int(request.form.get('mismatch_score', -1))
        indel_score = int(request.form.get('indel_score', -2))
        ruler = int(request.form.get('ruler', 10))

        seq1 = Sequence(seq1_name, seq1_sequence)
        seq2 = Sequence(seq2_name, seq2_sequence)

        score_scheme = Scheme(match_score, mismatch_score, indel_score)

        if selected_algorithm == 'nwa':
            aligner = NWA(seq1, seq2, score_scheme, allAlignment=show_all_alignments, ruler=ruler)
        elif selected_algorithm == 'swa':
            aligner = SWA(seq1, seq2, score_scheme, allAlignment=show_all_alignments, ruler=ruler)

        return render_template('index.html', seq1_name=seq1_name, seq1_sequence=seq1_sequence, seq2_name=seq2_name, seq2_sequence=seq2_sequence, alignments=aligner.alignment, selected_algorithm=selected_algorithm, show_all_alignments=show_all_alignments, match_score=match_score, mismatch_score=mismatch_score, indel_score=indel_score, ruler=ruler)

    return render_template('index.html', match_score=2, mismatch_score=-1, indel_score=-2, ruler=10)

if __name__ == '__main__':
    app.run(debug=True)