from flask import Flask, render_template, request
from aligner import Sequence, NWA

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        seq1_name = request.form['seq1_name']
        seq1_sequence = request.form['seq1_sequence']
        seq2_name = request.form['seq2_name']
        seq2_sequence = request.form['seq2_sequence']

        seq1 = Sequence(seq1_name, seq1_sequence)
        seq2 = Sequence(seq2_name, seq2_sequence)

        aligner = NWA(seq1, seq2)
        alignments = aligner.alignment
        
        return render_template('index.html', seq1_name=seq1_name, seq1_sequence=seq1_sequence, seq2_name=seq2_name, seq2_sequence=seq2_sequence, alignments=alignments)

    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)