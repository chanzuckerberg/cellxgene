import os
from flask import Blueprint, render_template, send_from_directory, current_app


bp = Blueprint("webapp", __name__, template_folder="templates")


@bp.route("/")
def index():
    dataset_title = current_app.config["DATASET_TITLE"]
    scripts = current_app.config["SCRIPTS"]
    return render_template("index.html", datasetTitle=dataset_title, SCRIPTS=scripts)


@bp.route("/favicon.png")
def favicon():
    return send_from_directory(os.path.join(bp.root_path, "static/img/"), "favicon.png")
