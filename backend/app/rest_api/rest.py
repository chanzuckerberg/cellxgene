from flask import (
    Blueprint, request
)
from flask_restful_swagger_2 import Api, swagger, Resource

from ..util.utils import make_payload
from ..util.filter import parse_filter


class InitializeAPI(Resource):
    @swagger.doc({
        'summary': 'get metadata schema, ranges for values, and cell count to initialize cellxgene app',
        'tags': ['initialize'],
        'parameters': [],
        'responses': {
            '200': {
                'description': 'initialization data for UI',
                'examples': {
                    'application/json': {
                        "data": {
                            "cellcount": 3589,
                            "options": {
                                "Sample.type": {
                                    "options": {
                                        "Glioblastoma": 3589
                                    }
                                },
                                "Selection": {
                                    "options": {
                                        "Astrocytes(HEPACAM)": 714,
                                        "Endothelial(BSC)": 123,
                                        "Microglia(CD45)": 1108,
                                        "Neurons(Thy1)": 685,
                                        "Oligodendrocytes(GC)": 294,
                                        "Unpanned": 665
                                    }
                                },
                                "Splice_sites_AT.AC": {
                                    "range": {
                                        "max": 1025,
                                        "min": 152
                                    }
                                },
                                "Splice_sites_Annotated": {
                                    "range": {
                                        "max": 1075869,
                                        "min": 26
                                    }
                                }
                            },
                            "schema": {
                                "CellName": {
                                    "displayname": "Name",
                                    "type": "string",
                                    "variabletype": "categorical"
                                },
                                "Class": {
                                    "displayname": "Class",
                                    "type": "string",
                                    "variabletype": "categorical"
                                },
                                "ERCC_reads": {
                                    "displayname": "ERCC Reads",
                                    "type": "int",
                                    "variabletype": "continuous"
                                },
                                "ERCC_to_non_ERCC": {
                                    "displayname": "ERCC:Non-ERCC",
                                    "type": "float",
                                    "variabletype": "continuous"
                                },
                                "Genes_detected": {
                                    "displayname": "Genes Detected",
                                    "type": "int",
                                    "variabletype": "continuous"
                                }
                            },
                            "genes": ["1/2-SBSRNA4", "A1BG", "A1BG-AS1"]

                        },
                        "status": {
                            "error": False,
                            "errormessage": ""
                        }
                    }
                }
            }
        }
    })
    def get(self):
        from app import data, REACTIVE_LIMIT
        return make_payload({
            "schema": data.schema,
            "cellcount": data.cell_count,
            "reactivelimit": REACTIVE_LIMIT,
            "genes": data.genes(),
            "ranges": data.metadata_ranges(),

        })


class CellsAPI(Resource):
    @swagger.doc({
        'summary': 'filter based on metadata fields to get a subset cells, expression data, and metadata',
        'tags': ['cells'],
        'description': "Cells takes query parameters defined in the schema retrieved from the /initialize enpoint. "
                       "<br>For categorical metadata keys filter based on `key=value` <br>"
                       " For continuous metadata keys filter by `key=min,max`<br> Either value "
                       "can be replaced by a \*. To have only a minimum value `key=min,\*`  To have only a maximum "
                       "value `key=\*,max` <br>Graph data (if retrieved) is normalized"
                       " To only retrieve cells that don't have a value for the key filter by `key`",
        'parameters': [],

        'responses': {
            '200': {
                'description': 'initialization data for UI',
                'examples': {
                    'application/json': {
                        "data": {
                            "badmetadatacount": 0,
                            "cellcount": 0,
                            "cellids": ["..."],
                            "metadata": [
                                {
                                    "CellName": "1001000173.G8",
                                    "Class": "Neoplastic",
                                    "Cluster_2d": "11",
                                    "Cluster_2d_color": "#8C564B",
                                    "Cluster_CNV": "1",
                                    "Cluster_CNV_color": "#1F77B4",
                                    "ERCC_reads": "152104",
                                    "ERCC_to_non_ERCC": "0.562454470489481",
                                    "Genes_detected": "1962",
                                    "Location": "Tumor",
                                    "Location.color": "#FF7F0E",
                                    "Multimapping_reads_percent": "2.67",
                                    "Neoplastic": "Neoplastic",
                                    "Non_ERCC_reads": "270429",
                                    "Sample.name": "BT_S2",
                                    "Sample.name.color": "#AEC7E8",
                                    "Sample.type": "Glioblastoma",
                                    "Sample.type.color": "#1F77B4",
                                    "Selection": "Unpanned",
                                    "Selection.color": "#98DF8A",
                                    "Splice_sites_AT.AC": "102",
                                    "Splice_sites_Annotated": "122397",
                                    "Splice_sites_GC.AG": "761",
                                    "Splice_sites_GT.AG": "125741",
                                    "Splice_sites_non_canonical": "56",
                                    "Splice_sites_total": "126660",
                                    "Total_reads": "1741039",
                                    "Unique_reads": "1400382",
                                    "Unique_reads_percent": "80.43",
                                    "Unmapped_mismatch": "2.15",
                                    "Unmapped_other": "0.18",
                                    "Unmapped_short": "14.56",
                                    "housekeeping_cluster": "2",
                                    "housekeeping_cluster_color": "#AEC7E8",
                                    "recluster_myeloid": "NA",
                                    "recluster_myeloid_color": "NA"
                                },
                            ],
                            "reactive": True,
                            "graph": [
                                [
                                    "1001000173.G8",
                                    0.93836,
                                    0.28623
                                ],

                                [
                                    "1001000173.D4",
                                    0.1662,
                                    0.79438
                                ]

                            ],
                            "status": {
                                "error": False,
                                "errormessage": ""
                            }

                        },
                    }
                },
            },

            '400': {
                'description': 'bad query params',
            }
        }
    })
    def get(self):
        from app import data
        payload = {
            "cellids": [],
            "metadata": [],
            "cellcount": 0,
            "graph": [],
            "ranges": {},
        }
        # get query params
        filter = parse_filter(request.args, data.schema)
        filtered_data = data.filter_cells(filter)
        payload["metadata"] = list(data.metadata(filtered_data))
        payload["ranges"] = data.metadata_ranges(filtered_data)
        payload["cellids"] = filtered_data
        payload["cellcount"] = len(payload["cellids"])
        payload["graph"] = list(data.create_graph(filtered_data))
        return make_payload(payload)



def get_api_resources():
    bp = Blueprint('api', __name__, url_prefix='/api/v2.0')
    api = Api(bp, add_api_spec_resource=False)
    api.add_resource(InitializeAPI, "/initialize")
    api.add_resource(CellsAPI, "/cells")
    return api