set -eo pipefail
flake8 server
npm run --prefix client/ build
npm run --prefix client/ unit-test
pytest -s server/test
