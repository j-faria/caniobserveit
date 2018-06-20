from flask import Flask
from flask import request, render_template

app = Flask(__name__)


import json
obs = json.load(open('data/observatories.txt'))
obs = ','.join([o[1]['name'] for o in obs.items()])

@app.route('/')
def hello():
    return render_template('index.html', obs=obs)

# @app.route('/<star>')
# def hello_star(star):
#     return "Hello {}!".format(star)

@app.route('/', methods=['GET', 'POST'])
def my_form_post():
    text = request.form
    print(text)
    processed_text = text['target'].upper()
    return processed_text


if __name__ == '__main__':
    app.debug = True
    app.run()
