from flask import Flask, render_template, request, send_from_directory
from aligner import Sequence, NWA, SWA, Scheme
import os

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'  # Define upload folder for storing results

if not os.path.exists(app.config['UPLOAD_FOLDER']):
    os.makedirs(app.config['UPLOAD_FOLDER'])

def process_blocks(seq1, seq2, block_size=100):  # Mock implementation of LSA
    # Simulated output for demonstration
    return f"Processed {seq1.name} vs {seq2.name} with block size {block_size}"

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        seq1_name = request.form['seq1_name']
        seq1_sequence = request.form['seq1_sequence']
        seq2_name = request.form['seq2_name']
        seq2_sequence = request.form['seq2_sequence']
        selected_algorithm = request.form['algorithm']

        match_score = int(request.form.get('match_score', 2))
        mismatch_score = int(request.form.get('mismatch_score', -1))
        indel_score = int(request.form.get('indel_score', -2))
        ruler = int(request.form.get('ruler', 10))

        seq1 = Sequence(seq1_name, seq1_sequence)
        seq2 = Sequence(seq2_name, seq2_sequence)
        score_scheme = Scheme(match_score, mismatch_score, indel_score)

        filepath = None
        if len(seq1_sequence) > 100 or len(seq2_sequence) > 100:
            result = process_blocks(seq1, seq2)
            filename = 'lsa_result.txt'
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            with open(filepath, 'w') as f:
                f.write(result)
        else:
            if selected_algorithm == 'nwa':
                aligner = NWA(seq1, seq2, score_scheme, allAlignment=False, ruler=ruler)
            elif selected_algorithm == 'swa':
                aligner = SWA(seq1, seq2, score_scheme, allAlignment=False, ruler=ruler)
            return render_template('index.html', alignments=aligner.alignment)

        return render_template('index.html', filepath=filepath, seq1_name=seq1_name, seq2_name=seq2_name)

    return render_template('index.html')

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
