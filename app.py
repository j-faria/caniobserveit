from flask import Flask
app = Flask(__name__)


@app.route('/')
def hello():
    return "Hello World!"


@app.route('/<star>')
def hello_star(star):
    return "Hello {}!".format(star)

if __name__ == '__main__':
    app.run()