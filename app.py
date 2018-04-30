from flask import Flask
from flask import request, render_template

app = Flask(__name__)

@app.route('/')
def hello():
    return render_template('search-form.html')

# @app.route('/<star>')
# def hello_star(star):
#     return "Hello {}!".format(star)

@app.route('/', methods=['POST'])
def my_form_post():
    text = request.form['target']
    processed_text = text.upper()
    return processed_text


if __name__ == '__main__':
    app.debug = True
    app.run()