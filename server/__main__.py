# Work around bug https://github.com/pallets/werkzeug/issues/461
if __package__ is None:
    import sys
    from pathlib import Path
    PKG_PATH = Path(__file__).parent
    sys.path.insert(0, str(PKG_PATH.parent))
    import server
    __package__ = PKG_PATH.name

# Main thing
from .cli.cli import cli
cli()
