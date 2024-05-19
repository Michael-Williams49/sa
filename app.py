from flask import Flask, render_template, request
from aligner import Sequence, NWA, SWA

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

        seq1 = Sequence(seq1_name, seq1_sequence)
        seq2 = Sequence(seq2_name, seq2_sequence)

        if selected_algorithm == 'nwa':
            aligner = NWA(seq1, seq2, allAlignment = show_all_alignments)
        elif selected_algorithm == 'swa':
            aligner = SWA(seq1, seq2, allAlignment = show_all_alignments)
        else:
            # Handle invalid algorithm selection
            pass

        alignments = aligner.alignment

        return render_template('index.html', seq1_name=seq1_name, seq1_sequence=seq1_sequence, seq2_name=seq2_name, seq2_sequence=seq2_sequence, alignments=alignments, selected_algorithm=selected_algorithm, show_all_alignments=show_all_alignments)

    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)