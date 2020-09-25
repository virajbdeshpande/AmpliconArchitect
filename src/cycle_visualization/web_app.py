from flask import Flask, render_template, request, session, redirect
import uuid
import os

from draw_episome import EpisomeDrawer
from breakpoint_graph import graph_decomposition



app = Flask(__name__)
app.secret_key = b'7\x0c\xd0\xdd~[\x9d_\xb4\xda\xaf\xae\xd0\x0b6\x8d\xfe\xd9\xf13Z\x1e\xedS'

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(APP_ROOT, 'static/uploads')
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


def draw_structures(true_structure, reconstructed_structures, coverage_image_file_name, min_copy_number=0, merge_error=None, pivot_error=None):
    if 'min_copy_number' in session:
        min_copy_number = session['min_copy_number']
    ed = EpisomeDrawer(image_output=False, min_copy_number_to_show=min_copy_number)

    reconstructed_structure = reconstructed_structures[-1] if reconstructed_structures else None
    input_files = []
    if reconstructed_structure:
        input_files.append(reconstructed_structure)
    if true_structure:
        input_files.append(true_structure)
    ed.draw_episome(input_files)

    session['true_structure'] = true_structure
    session['reconstructed_structures'] = reconstructed_structures
    session['coverage_image'] = coverage_image_file_name
    session['min_copy_number'] = min_copy_number
    show_image = coverage_image_file_name != None
    edited = len(reconstructed_structures) > 1
    image_src = 'uploads/' + coverage_image_file_name if coverage_image_file_name else None
    return render_template('episome_template.html', rectangles=ed.rectangles, horizontal_lines=ed.horizontal_lines, file_names=ed.file_names,
                           vertical_lines=ed.vertical_lines, arrows=ed.arrows, texts=ed.texts, axes_labels_texts=ed.axes_labels_texts,
                           cycle_show=ed.cycle_show, cycles_sections=ed.cycle_sections, reconstructed_cycles=ed.reconstructed_cycles,
                           reconstructed_segments=ed.reconstructed_segments, min_copy_number=min_copy_number, edited=edited, merge_error=merge_error,
                           show_image=show_image, image_src=image_src, pivot_error=pivot_error)


@app.route('/', methods=['GET', 'POST'])
def main():
    reconstructed_structure = None
    true_structure = None
    coverage_image = None
    session.pop('reconstructed_structures', None)
    session.pop('true_structures', None)
    session.pop('coverage_image', None)
    session.pop('min_copy_number', None)
    if request.method == 'POST':
        if 'reconstructed_structure' not in request.files and 'true_structure' not in request.files:
            return render_template('episome_template.html', upload_error=True)
        if request.files['reconstructed_structure'].filename == '' and request.files['true_structure'].filename == '':
            return render_template('episome_template.html', upload_error=True)
        if 'true_structure' in request.files and request.files['true_structure'].filename != '':
            f = request.files['true_structure']
            true_structure_name = f.filename
            content = f.read().decode()
            true_structure = (true_structure_name, content)
        if 'reconstructed_structure' in request.files and request.files['reconstructed_structure'].filename != '':
            f = request.files['reconstructed_structure']
            reconstructed_structure_name = f.filename
            content = str(f.read().decode())
            reconstructed_structure = [(reconstructed_structure_name, content)]
        if 'coverage_image' in request.files and request.files['coverage_image'].filename != '':
            f = request.files['coverage_image']
            coverage_image = f.filename.strip()
            f.save(os.path.join(app.config['UPLOAD_FOLDER'], coverage_image))

        return draw_structures(true_structure, reconstructed_structure, coverage_image)
    else:
        return render_template('episome_template.html', upload_error=True)

@app.route('/update', methods=['POST', 'GET'])
def update():
    if 'reconstructed_structures' not in session:
        return redirect('/')
    true_structure = session['true_structure']
    reconstructed_structures = session['reconstructed_structures']
    coverage_image_file_name = session['coverage_image']

    if 'min_copy_number' in request.form and request.form['min_copy_number'] != '':
        session['min_copy_number'] = float(request.form['min_copy_number'])
    return draw_structures(true_structure, reconstructed_structures, coverage_image_file_name)

@app.route('/merge', methods=['POST'])
def merge():
    print('merge function')
    if 'reconstructed_structures' not in session:
        print('redirecting to home')
        return redirect('/')
    true_structure = session['true_structure']
    reconstructed_structures = session['reconstructed_structures']
    coverage_image_file_name = session['coverage_image']

    if reconstructed_structures is None or len(reconstructed_structures) == 0:
        print('redirect to home 2')
        return redirect('/')

    first_cycle = str(request.form['first_cycle'])
    second_cycle = str(request.form['second_cycle'])
    first_segment = int(request.form['first_segment'])
    second_segment = int(request.form['second_segment'])
    last_structure = reconstructed_structures[-1]

    merge_error = None
    try:
        g = graph_decomposition(file_content=last_structure[1])
        g.merge(first_cycle, second_cycle, first_segment, second_segment)
        new_reconstructed_structure = (last_structure[0], str(g))
        reconstructed_structures.append(new_reconstructed_structure)
    except Exception as e:
        merge_error = str(e)

    return draw_structures(true_structure, reconstructed_structures, coverage_image_file_name, merge_error=merge_error)

@app.route('/pivot', methods=['POST'])
def pivot():
    if 'reconstructed_structures' not in session:
        return redirect('/')
    true_structure = session['true_structure']
    reconstructed_structures = session['reconstructed_structures']
    coverage_image_file_name = session['coverage_image']

    if reconstructed_structures is None or len(reconstructed_structures) == 0:
        return redirect('/')

    first_cycle = str(request.form['first_cycle'])
    first_segment = int(request.form['first_segment'])
    second_segment = int(request.form['second_segment'])
    last_structure = reconstructed_structures[-1]

    pivot_error = None
    try:
        g = graph_decomposition(file_content=last_structure[1])
        g.pivot(first_cycle, first_segment, second_segment)
        new_reconstructed_structure = (last_structure[0], str(g))
        reconstructed_structures.append(new_reconstructed_structure)
    except Exception as e:
        pivot_error = str(e)

    return draw_structures(true_structure, reconstructed_structures, coverage_image_file_name, pivot_error=pivot_error)

@app.route('/undo_edit', methods=['GET'])
def undo_edit():
    if 'reconstructed_structures' not in session:
        return redirect('/')
    true_structure = session['true_structure']
    reconstructed_structures = session['reconstructed_structures']
    coverage_image_file_name = session['coverage_image']

    if reconstructed_structures is None or len(reconstructed_structures) <= 1:
        return redirect('/')

    reconstructed_structures = reconstructed_structures[:-1]

    return draw_structures(true_structure, reconstructed_structures, coverage_image_file_name)
