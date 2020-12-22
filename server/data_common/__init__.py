# import the data adaptor types so they can be registered

import server.data_cxg.cxg_adaptor  # noqa: F401
import server.data_anndata.anndata_adaptor  # noqa: F401

# FIXME: may want to make this one a plugin so that hosted CZI specific
# details are not in the public repo.  For development it is easier to
# keep here for now.
import server.data_cxg.czi_cxg_adaptor  # noqa: F401

