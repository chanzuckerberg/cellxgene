# -*- coding: utf-8 -*-

from server import __version__ as cellxgene_version


class AppFeature(object):
    def __init__(self, path, available=False, method="POST", extra={}):
        self.path = path
        self.available = available
        self.method = method
        self.extra = extra
        for k, v in extra.items():
            setattr(self, k, v)

    def todict(self):
        d = dict(
            available=self.available,
            method=self.method,
            path=self.path)
        d.update(self.extra)
        return d


class AppConfig(object):

    def __init__(self, **kw):
        # app inputs
        self.title = ""
        self.about = None
        self.scripts = []
        self.layout = None
        self.max_category_items = 100
        self.diffexp_lfc_cutoff = 0.01
        self.obs_names = None
        self.var_names = None
        self.scanpy_backed = False
        self.disable_diffexp = False

        inputs = ["title", "about", "scripts", "layout",
                  "max_category_items", "diffexp_lfc_cutoff",
                  "obs_names", "var_names",
                  "scanpy_backed", "disable_diffexp"]
        for k, v in kw.items():
            if k in inputs:
                setattr(self, k, v)
            else:
                raise RuntimeError(f"unknown config parameter {k}.")

        # parameters
        self.diffexp_may_be_slow = False

    def get_config(self, data_engine, annotation=None):

        # FIXME The current set of config is not consistently presented:
        # we have camalCase, hyphen-text, and underscore_text

        # features
        features = [f.todict() for f in data_engine.get_features().values()]

        # display_names
        display_names = dict(
            engine=data_engine.get_name(),
            dataset=self.title)

        # library_versions
        library_versions = {}
        library_versions.update(data_engine.get_library_versions())
        library_versions["cellxgene"] = cellxgene_version

        # links
        links = {"about-dataset" : self.about}

        # parameters
        parameters = {
            "layout": self.layout,
            "max_category_items": self.max_category_items,
            "obs_names": None,
            "var_names": None,
            "diffexp_lfc_cutoff": self.diffexp_lfc_cutoff,
            "annotations": False,
            "annotations_file": None,
            "annotations_output_dir": None,
            "annotations_cell_ontology_enabled": False,
            "annotations_cell_ontology_obopath": None,
            "annotations_cell_ontology_terms": None,
            "backed": False,
            "disable_diffexp": False,
            "diffexp_may_be_slow": False,
        }

        data_engine.update_parameters(parameters)
        if annotation:
            annotation.update_parameters(parameters)

        # gather it all together
        c = {}
        config = c["config"] = {}
        config["features"] = features
        config["displayNames"] = display_names
        config["library_versions"] = library_versions
        config["links"] = links
        config["parameters"] = parameters

        return c
