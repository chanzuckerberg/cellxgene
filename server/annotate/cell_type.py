from dataclasses import dataclass
from pathlib import Path
from typing import Set, Tuple, Dict
import scvi


from anndata import AnnData
from scvi.models import SCANVI


def extract_tissue_types(dataset: AnnData) -> Set[str]:
    """
    Returns distinct set of tissue type ontology term IDs.
    """
    pass


class CellTypeTissueModel(dataclass):
    cell_type_ontology_term_id: str
    tissue_type_ontology_term_id: str
    model: SCANVI


def retrieve_models(tissue_types: Set[str]) -> Dict[str, CellTypeTissueModel]:
    pass


def annotate(dataset: AnnData, models: Dict[str, CellTypeTissueModel], output_h5ad_file: Path):
    # TODO: necessary to partition into AnnData-per-tissue, then invoke model on each partitioned AnnData object, then reconstitute into a single AnnData
    # TODO: record model metadata (version) in AnnData
    return dataset
