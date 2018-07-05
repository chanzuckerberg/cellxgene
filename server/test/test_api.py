import unittest
import requests
import json


class EndPoints(unittest.TestCase):
    """Test Case for endpoints"""

    def setUp(self):
        # Local
        self.url_base = "http://0.0.0.0:5005/api/" + "v0.1/"
        self.session = requests.Session()

    def test_cells(self):
        url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
            ["louvain=B cells"]))
        result = self.session.get(url)
        assert result.status_code == 200
        result_data = result.json()
        assert "B cells" in result_data["data"]["ranges"]["louvain"]["options"]
        url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
            ["louvain=B cells", "louvain=Megakaryocytes"]))
        result = self.session.get(url)
        assert result.status_code == 200
        result_data = result.json()
        assert "Megakaryocytes" in result_data["data"]["ranges"]["louvain"]["options"]

    def test_initialize(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="initialize")
        result = self.session.get(url)
        assert result.status_code == 200
        result_data = result.json()
        assert result_data["data"]["cellcount"] == 2638
        assert len(result_data["data"]['ranges']['CellName']['options'])  == 2638


    def test_expression_get(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="expression")
        result = self.session.get(url)
        assert result.status_code == 200

    def test_expression_post(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="expression")
        result = self.session.post(url, data=json.dumps({"celllist": ["AAACATACAACCAC-1", "AACCGATGGTCATG-1"], "genelist":	["BACH1", "MIS18A", "ATP5O"]}), headers={'content-type': 'application/json'})
        assert result.status_code == 200
        result_data = result.json()
        assert len(result_data["data"]["cells"]) == 2
        assert len(result_data["data"]["cells"][0]['e']) == 3

    def test_diffexp(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="diffexp")
        result = self.session.post(url, data=json.dumps({"celllist1": ["AAACATACAACCAC-1", "AACCGATGGTCATG-1"], "celllist2": ["CCGATAGACCTAAG-1", "GGTGGAGAAGTAGA-1"]}), headers={'content-type': 'application/json'})
        assert result.status_code == 200
        