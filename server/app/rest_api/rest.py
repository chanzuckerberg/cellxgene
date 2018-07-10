from flask import (
    Blueprint, request
)
from flask_restful_swagger_2 import Api, swagger, Resource

from ..util.utils import make_payload
from ..util.filter import parse_filter


class InitializeAPI(Resource):
    @swagger.doc({
        "summary": "get metadata schema, ranges for values, and cell count to initialize cellxgene app",
        "tags": ["initialize"],
        "parameters": [],
        "responses": {
            "200": {
                "description": "initialization data for UI",
                "examples": {
                    "application/json": {
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
        from ..app import data, REACTIVE_LIMIT
        return make_payload({
            "schema": data.schema,
            "cellcount": data.cell_count,
            "reactivelimit": REACTIVE_LIMIT,
            "genes": data.genes(),
            "ranges": data.metadata_ranges(),

        })


class CellsAPI(Resource):
    @swagger.doc({
        "summary": "filter based on metadata fields to get a subset cells, expression data, and metadata",
        "tags": ["cells"],
        "description": "Cells takes query parameters defined in the schema retrieved from the /initialize enpoint. "
                       "<br>For categorical metadata keys filter based on `key=value` <br>"
                       " For continuous metadata keys filter by `key=min,max`<br> Either value "
                       "can be replaced by a \*. To have only a minimum value `key=min,\*`  To have only a maximum "
                       "value `key=\*,max` <br>Graph data (if retrieved) is normalized"
                       " To only retrieve cells that don't have a value for the key filter by `key`",
        "parameters": [],

        "responses": {
            "200": {
                "description": "initialization data for UI",
                "examples": {
                    "application/json": {
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

            "400": {
                "description": "bad query params",
            }
        }
    })
    def get(self):
        from ..app import data
        payload = {
            "metadata": [],
            "cellcount": 0,
            "graph": [],
            "ranges": {},
        }
        # get query params
        cells_filter = parse_filter(request.args, data.schema)
        filtered_data = data.filter_cells(cells_filter)
        payload["metadata"] = data.metadata(filtered_data)
        payload["ranges"] = data.metadata_ranges(filtered_data)
        payload["graph"] = data.create_graph(filtered_data)
        payload["cellcount"] = data.cell_count
        return make_payload(payload)


class ExpressionAPI(Resource):
    @swagger.doc({
        "summary": "Json with gene list and expression data by cell, limited to first 40 cells",
        "tags": ["expression"],
        "parameters": [
            {
                "name": "include_unexpressed_genes",
                "description": "Include genes that have 0 expression across all cells in set",
                "in": "path",
                "type": "bool",
            }
        ],
        "responses": {
            "200": {
                "description": "Json for heatmap",
                "examples": {
                    "application/json": {
                        "data": {
                            "cells": [
                                {
                                    "cellname": "1/2-SBSRNA4",
                                    "e": [0, 0, 214, 0, 0]
                                },
                            ],
                            "genes": [
                                "1001000173.G8",
                                "1001000173.D4",
                                "1001000173.B4",
                                "1001000173.A2",
                                "1001000173.E2"
                            ],
                            "nonzero_gene_count": 2857
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
        from ..app import data
        expression_data = data.expression()
        return make_payload(expression_data)

    @swagger.doc({
        "summary": "Json with gene list and expression data by cell",
        "tags": ["expression"],
        "parameters": [
            {
                "name": "body",
                "in": "body",
                "schema": {
                    "example": {
                        "celllist": ["1001000173.G8", "1001000173.D4"],
                        "genelist": ["1/2-SBSRNA4", "A1BG", "A1BG-AS1", "A1CF", "A2LD1", "A2M", "A2ML1", "A2MP1",
                                     "A4GALT"],
                        "include_unexpressed_genes": True,
                    }

                }
            },
        ],
        "responses": {
            "200": {
                "description": "Json for expressiondata",
                "examples": {
                    "application/json": {
                        "data": {
                            "cells": [
                                {
                                    "cellname": "1001000173.D4",
                                    "e": [0, 0]
                                },
                                {
                                    "cellname": "1001000173.G8",
                                    "e": [0, 0]
                                }
                            ],
                            "genes": [
                                "ABCD4",
                                "ZWINT"
                            ],
                            "nonzero_gene_count": 2857
                        },
                        "status": {
                            "error": False,
                            "errormessage": ""
                        }

                    }
                }
            },
            "400": {
                "description": "Required parameter missing/incorrect",
            }
        }
    })
    def post(self):
        from ..app import data
        args = request.get_json()
        cell_list = args.get("celllist", [])
        gene_list = args.get("genelist", [])
        if not cell_list and not gene_list:
            return make_payload([], "must include celllist and/or genelist parameter", 400)

        expression_data = data.expression(cell_list, gene_list)

        if cell_list and len(expression_data["cells"]) < len(cell_list):
            return make_payload([], "Some cell ids not available", 400)
        if gene_list and len(expression_data["genes"]) < len(gene_list):
            return make_payload([], "Some genes not available", 400)

        return make_payload(expression_data)


class DifferentialExpressionAPI(Resource):
    @swagger.doc({
        "summary": "Get the top expressed genes for two cell sets. Calculated using t-test",
        "tags": ["expression"],
        "parameters": [
            {
                "name": "body",
                "in": "body",
                "schema": {
                    "example": {
                        "celllist1": ["1001000176.C12", "1001000176.C7", "1001000177.F11"],
                        "celllist2": ["1001000012.D2", "1001000017.F10", "1001000033.C3", "1001000229.D4"],
                        "num_genes": 5,
                        "pval": 0.000001,
                    },
                }
            }
        ],
        "responses": {
            "200": {
                "description": "top expressed genes for cellset1, cellset2",
                "examples": {
                    "application/json": {
                        "data": {
                            "celllist1": {
                                "ave_diff": [
                                    432.0132935431362,
                                    12470.5623982637,
                                    957.0246880086814
                                ],
                                "mean_expression_cellset1": [
                                    438.6185567010309,
                                    13315.536082474227,
                                    1076.5773195876288
                                ],
                                "mean_expression_cellset2": [
                                    6.605263157894737,
                                    844.9736842105264,
                                    119.55263157894737
                                ],
                                "pval": [
                                    3.8906598089944563e-35,
                                    1.9086226376018916e-25,
                                    7.847480544069826e-21
                                ],
                                "topgenes": [
                                    "TMSB10",
                                    "FTL",
                                    "TMSB4X"
                                ]
                            },
                            "celllist2": {
                                "ave_diff": [
                                    -6860.599158979924,
                                    -519.1314432989691,
                                    -10278.328269126423
                                ],
                                "mean_expression_cellset1": [
                                    2.8350515463917527,
                                    0.6185567010309279,
                                    23.09278350515464
                                ],
                                "mean_expression_cellset2": [
                                    6863.434210526316,
                                    519.75,
                                    10301.421052631578
                                ],
                                "pval": [
                                    4.662891833748732e-44,
                                    3.6278087029927103e-37,
                                    8.396825170618402e-35
                                ],
                                "topgenes": [
                                    "SPARCL1",
                                    "C1orf61",
                                    "CLU"
                                ]
                            }
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
    def post(self):
        from ..app import data
        args = request.get_json()
        cell_list_1 = args.get("celllist1", [])
        cell_list_2 = args.get("celllist2", [])
        num_genes = args.get("num_genes", 7)
        pval = args.get("pval", 0.5)
        if not (cell_list_1 and cell_list_2):
            return make_payload([],
                                "must include celllist1 and celllist2 parameters",
                                400)
        data = data.diffexp(cell_list_1, cell_list_2, pval, num_genes)
        return make_payload(data)


def get_api_resources():
    bp = Blueprint("api", __name__, url_prefix="/api/v0.1")
    api = Api(bp, add_api_spec_resource=False)
    api.add_resource(InitializeAPI, "/initialize")
    api.add_resource(CellsAPI, "/cells")
    api.add_resource(ExpressionAPI, "/expression")
    api.add_resource(DifferentialExpressionAPI, "/diffexpression")
    return api
