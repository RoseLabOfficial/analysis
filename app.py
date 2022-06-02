from flask import Flask, render_template
from flask_wtf.csrf import CSRFProtect
from flask_socketio import SocketIO
import random
import string

csrf = CSRFProtect()

app = Flask(__name__, static_folder="./server/static", template_folder="./server/templates")
app.config["UPLOAD_FOLDER"] = "./server/uploads"
_host = "localhost"
_port = 8080
_allowed_extensions=["txt", "csv"]
app.secret_key = "".join(random.choices(string.digits + string.ascii_uppercase + string.ascii_lowercase, k=100))
csrf.init_app(app)
socketio = SocketIO(app)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in _allowed_extensions


@app.route("/analysis")
def analysis():
    return render_template("analysis.html")