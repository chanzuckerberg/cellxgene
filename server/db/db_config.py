import os


class CellxGeneDbConfig:
    def __init__(self, *args, **kwargs):
        component_name="backend",
        # secret_name=f"database{'_local' if 'CELLXGENE_LOCAL_DEV' in os.environ else ''}",

    def database_uri(self):
        return "postgresql://postgres:test_pw@localhost:5432"
