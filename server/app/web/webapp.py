from flask import (
    Blueprint, render_template, url_for, current_app
)

bp = Blueprint("webapp", __name__, template_folder="templates")


@bp.route("/")
def index():
    url_base = current_app.config["CXG_API_BASE"]
    dataset_title = current_app.config["DATASET_TITLE"]
    return render_template("index.html", prefix=url_base, datasetTitle=dataset_title)


# renders swagger documentation
@bp.route("/swagger")
def swag():
    return render_template("swagger.html")


# renders swagger documentation
@bp.route("/favicon.png")
def favicon():
    return url_for("static", filename="img/favicon.png")
