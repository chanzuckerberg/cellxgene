set -eo pipefail
make lint
npm run --prefix client/ build
npm run --prefix client/ unit-test
pytest -s server/test
