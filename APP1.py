from flask import Flask, render_template, request, send_from_directory
import os
import re
from aligner import Sequence, NWA, SWA, Scheme
from LSA2 import process_blocks  # Ensure this module is properly implemented

app = Flask(__name__)

# Directory to save alignment results
alignment_dir = os.path.join(app.root_path, 'alignments')
os.makedirs(alignment_dir, exist_ok=True)


@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        seq1_name = request.form['seq1_name']
        seq1_sequence = request.form['seq1_sequence']
        seq2_name = request.form['seq2_name']
        seq2_sequence = request.form['seq2_sequence']
        show_all_alignments = request.form.get('show_all_alignments') == 'on'
        selected_algorithm = request.form['algorithm']

        match_score = int(request.form.get('match_score', 2))
        mismatch_score = int(request.form.get('mismatch_score', -1))
        indel_score = int(request.form.get('indel_score', -2))
        ruler = int(request.form.get('ruler', 10))

        seq1 = Sequence(seq1_name, seq1_sequence)
        seq2 = Sequence(seq2_name, seq2_sequence)
        score_scheme = Scheme(match_score, mismatch_score, indel_score)

        if len(seq1_sequence) > 1000 and len(seq2_sequence) > 1000:
            alignmentA, alignmentB = process_blocks(
                seq1_sequence, seq2_sequence, 100)
            file_path = save_alignment_to_file(
                alignmentA, alignmentB, seq1_name, seq2_name)
            return render_template('alignment_download.html', file_path=file_path)

        if selected_algorithm == 'nwa':
            aligner = NWA(seq1, seq2, score_scheme,
                          allAlignment=show_all_alignments, ruler=ruler)
        elif selected_algorithm == 'swa':
            aligner = SWA(seq1, seq2, score_scheme,
                          allAlignment=show_all_alignments, ruler=ruler)
        alignments = aligner.alignment

        return render_template('index.html', seq1_name=seq1_name, seq1_sequence=seq1_sequence,
                               seq2_name=seq2_name, seq2_sequence=seq2_sequence, alignments=alignments,
                               selected_algorithm=selected_algorithm, show_all_alignments=show_all_alignments,
                               match_score=match_score, mismatch_score=mismatch_score,
                               indel_score=indel_score, ruler=ruler)

    return render_template('index.html', match_score=2, mismatch_score=-1, indel_score=-2, ruler=10)


@app.route('/download/<path:filename>', methods=['GET'])
def download(filename):
    """Serve a file from the alignment directory to the user."""
    return send_from_directory(directory=alignment_dir, filename=filename, as_attachment=True)


def save_alignment_to_file(alignmentA, alignmentB, seq1_name, seq2_name):
    safe_seq1_name = sanitize_filename(seq1_name)
    safe_seq2_name = sanitize_filename(seq2_name)
    file_name = f"{safe_seq1_name}_vs_{safe_seq2_name}_alignment.txt"
    file_path = os.path.join(alignment_dir, file_name)
    with open(file_path, 'w') as file:
        file.write(f"Alignment for {seq1_name} and {seq2_name}:\n")
        file.write("Sequence A (aligned):\n" + alignmentA + "\n")
        file.write("Sequence B (aligned):\n" + alignmentB + "\n")
    return file_path


def sanitize_filename(input_name):
    sanitized = re.sub(r"[|: ]", "_", input_name)
    sanitized = re.sub(r"[^a-zA-Z0-9_]", "", sanitized)
    return sanitized


if __name__ == '__main__':
    app.run(debug=True)
