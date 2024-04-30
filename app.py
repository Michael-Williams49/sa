from flask import Flask, request, jsonify, render_template
# Ensure these are correctly defined in NWAlignment.py
from NWAlignment import NWA, Sequence, Scheme

app = Flask(__name__)


@app.route('/align', methods=['POST'])
def align():
    data = request.json
    seq1 = data['seq1']
    seq2 = data['seq2']
    print("Received sequences:", seq1, seq2)  # 打印接收到的序列

    alignment = NWA(Sequence("Seq1", seq1), Sequence(
        "Seq2", seq2), Scheme(1, -1, -1))
    alignment_results = alignment.alignment
    print("Alignment results:", alignment_results)  # 打印对齐结果

    result_as_dicts = [[seq.to_dict() for seq in pair]
                       for pair in alignment_results]
    return jsonify({"result": result_as_dicts})


@app.route('/', endpoint='home')  # 显式指定端点名称
def home():
    return render_template('webpage.html')  # Render the HTML file


if __name__ == '__main__':
    app.run(debug=True)
