from flask import Flask, render_template, request, redirect, url_for, send_from_directory, session
import os
from aligner import Sequence, NWA, SWA, Scheme
from LSA3 import process_blocks

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = r'C:\path\to\upload'
app.secret_key = 'your_secret_key'


@app.route('/', methods=['GET', 'POST'])
    
    
    
def index():
    default_values = {
        'match_score': 2,
        'mismatch_score': -1,
        'indel_score': -2,
        'ruler': 10,
        'seq1_name': 'Example Seq1',
        'seq2_name': 'Example Seq2',
        
        'algorithm': 'nwa',
        'show_all_alignments': False
    }

     
     

    if request.method == 'POST':
        seq1_name = request.form['seq1_name']
        seq1_sequence = request.form['seq1_sequence']
        seq2_name = request.form['seq2_name']
        seq2_sequence = request.form['seq2_sequence']
        show_all_alignments = request.form.get(
            'show_all_alignments', False) == 'on'
        selected_algorithm = request.form['algorithm']

        match_score = int(request.form.get('match_score', 2))
        mismatch_score = int(request.form.get('mismatch_score', -1))
        indel_score = int(request.form.get('indel_score', -2))
        ruler = int(request.form.get('ruler', 10))

        seq1 = Sequence(seq1_name, seq1_sequence)
        seq2 = Sequence(seq2_name, seq2_sequence)
        score_scheme = Scheme(match_score, mismatch_score, indel_score)
        

        if len(seq1_sequence) > 1000 or len(seq2_sequence) > 1000:
            class alignment:
                def __init__(self):
                    self.alignment=[[],[],[],[]]
            aligner=alignment
            result = process_blocks(seq1_sequence, seq2_sequence, 100)
            filename = 'lsa_result.txt'
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            os.makedirs(os.path.dirname(filepath), exist_ok=True)

            if isinstance(result, tuple):
                # Format the output with sequence names and results
                formatted_result = f"{seq1_name}: {result[0]}\n{seq2_name}: {result[1]}"
            else:
                # Assuming result is formatted correctly if not a tuple
                formatted_result = result

            with open(filepath, 'w') as f:
                f.write(formatted_result)
            return redirect(url_for('uploaded_file', filename=filename))
        else:
         
            if selected_algorithm == 'nwa':
                aligner = NWA(seq1, seq2, score_scheme,
                              allAlignment=show_all_alignments, ruler=ruler)
            elif selected_algorithm == 'swa':
                aligner = SWA(seq1, seq2, score_scheme,
                              allAlignment=show_all_alignments, ruler=ruler)
            else:
                return "Invalid algorithm selected", 400
            return render_template('index.html', seq1_name=seq1_name, seq1_sequence=seq1_sequence, seq2_name=seq2_name, seq2_sequence=seq2_sequence, alignments=aligner.alignment, selected_algorithm=selected_algorithm, show_all_alignments=show_all_alignments, match_score=match_score, mismatch_score=mismatch_score, indel_score=indel_score, ruler=ruler)
        return render_template('index.html', seq1_name=seq1_name, seq1_sequence=seq1_sequence, seq2_name=seq2_name, seq2_sequence=seq2_sequence,  selected_algorithm=selected_algorithm, show_all_alignments=show_all_alignments, match_score=match_score, mismatch_score=mismatch_score, indel_score=indel_score, ruler=ruler)
  
        
    
    return render_template('index.html', match_score=2, mismatch_score=-1, indel_score=-2, ruler=10)


@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)



if __name__ == '__main__':
    app.run(debug=True)
