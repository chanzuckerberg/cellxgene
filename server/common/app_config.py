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
        d = dict(available=self.available, method=self.method, path=self.path)
        d.update(self.extra)
        return d


class AppConfig(object):
    def __init__(self, **kw):
        super().__init__()

        # app inputs
        self.datapath = None
        self.dataroot = None
        self.title = ""
        self.about = None
        self.scripts = []
        self.layout = []
        self.max_category_items = 100
        self.diffexp_lfc_cutoff = 0.01
        self.disable_diffexp = False
        self.enable_reembedding = False
        self.anndata_backed = False

        # TODO these options may not apply to all datasets in the multi dataset.
        # may need to invent a way to associate these config parameters with
        # specific datasets.
        self.obs_names = None
        self.var_names = None

        # parameters
        self.diffexp_may_be_slow = False

        inputs = [
            "datapath",
            "dataroot",
            "title",
            "about",
            "scripts",
            "layout",
            "max_category_items",
            "diffexp_lfc_cutoff",
            "obs_names",
            "var_names",
            "anndata_backed",
            "disable_diffexp",
            "enable_reembedding",
        ]

        self.update(inputs, kw)

    def update(self, inputs, kw):

        for k, v in kw.items():
            if k in inputs:
                setattr(self, k, v)
            else:
                raise RuntimeError(f"unknown config parameter {k}.")

    def get_title(self, data_adaptor):
        return self.title if self.title else data_adaptor.get_title()

    def get_about(self, data_adaptor):
        return self.about if self.about else data_adaptor.get_about()

    def get_config(self, data_adaptor, annotation=None):

        # FIXME The current set of config is not consistently presented:
        # we have camalCase, hyphen-text, and underscore_text

        # features
        features = [f.todict() for f in data_adaptor.get_features(annotation)]

        # display_names
        title = self.get_title(data_adaptor)
        about = self.get_about(data_adaptor)

        display_names = dict(engine=data_adaptor.get_name(), dataset=title)

        # library_versions
        library_versions = {}
        library_versions.update(data_adaptor.get_library_versions())
        library_versions["cellxgene"] = cellxgene_version

        # links
        links = {"about-dataset": about}

        # parameters
        parameters = {
            "layout": self.layout,
            "max-category-items": self.max_category_items,
            "obs_names": self.obs_names,
            "var_names": self.var_names,
            "diffexp_lfc_cutoff": self.diffexp_lfc_cutoff,
            "backed": self.anndata_backed,
            "disable-diffexp": self.disable_diffexp,
            "enable-reembedding": self.enable_reembedding,
            "annotations": False,
            "annotations_file": None,
            "annotations_dir": None,
            "annotations_cell_ontology_enabled": False,
            "annotations_cell_ontology_obopath": None,
            "annotations_cell_ontology_terms": None,
            "diffexp-may-be-slow": False,
        }

        data_adaptor.update_parameters(parameters)
        if annotation:
            annotation.update_parameters(parameters, data_adaptor)

        # gather it all together
        c = {}
        config = c["config"] = {}
        config["features"] = features
        config["displayNames"] = display_names
        config["library_versions"] = library_versions
        config["links"] = links
        config["parameters"] = parameters

        return c
