from flask import Flask
from flask import request, render_template
from visibilityPlot import visibilityPlot
import datetime as dt

app = Flask(__name__)


import json
obs = json.load(open('data/observatories.txt'))
obs = {v['name'].replace(',', ';'): k for (k, v) in obs.items()}
obs_names = ','.join(obs.keys())

@app.route('/')
def hello():
    return render_template('index.html', obs=obs_names)


@app.route('/', methods=['GET', 'POST'])
def my_form_post():
    form = request.form
    date = dt.datetime(*[int(i) for i in form['date'].split('-')])
    target = form['target']
    observatory = obs[form['observatory']]

    plot = visibilityPlot(date=date, target=target, observatory=observatory)
    return render_template('plot.html', plot=plot)


if __name__ == '__main__':
    app.debug = True
    app.run()
