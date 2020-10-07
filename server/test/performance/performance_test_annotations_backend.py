import json
import string
from contextlib import contextmanager
from timeit import default_timer
import concurrent.futures
import numpy as np
import requests
import sys
from server.data_common.fbs.matrix import encode_matrix_fbs
import pandas as pd
import random

"""
Before running, sign into the dataportal, copy the cookie and paste it below. To test in staging or prod update the
url base below. It is also possible to configure the number of categories created and the number of unique labels per
category.
"""

cookie = ""

test_datasets = {
    "smallest": {
        "dataset_url": "kampmann_lab_human_AD_snRNAseq_EC_inhibitoryNeurons-53-remixed.cxg",
        "name": "smallest",
        "num_cells": 5270,
    },
    "10k": {
        "dataset_url": "krasnow_lab_human_lung_cell_atlas_smartseq2-2-remixed.cxg",
        "name": "10k",
        "num_cells": 9409,
    },
    "80k": {
        "dataset_url": "Single_cell_gene_expression_profiling_of_SARS_CoV_2_infected_human_cell_lines_H1299-27-remixed.cxg",  # noqa E501
        "name": "80k",
        "num_cells": 81736,
    },
    "140k": {"dataset_url": "Single_cell_drug_screening_a549-42-remixed.cxg", "name": "140k", "num_cells": 143015},
    "largest": {"dataset_url": "human_cell_landscape.cxg", "name": "largest", "num_cells": 599926},
    "1million": {"dataset_url": None, "name": "1million", "num_cells": 1000000},
    "4million": {"dataset_url": None, "name": "4million", "num_cells": 4000000},
}

url_base = "https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/e/"
annotations_category_count = [1, 10, 50]
max_labels = [5, 50, 100]


class PerformanceTestingAnnotations:
    def __init__(
        self,
        datasets=test_datasets,
        annotations_category_count=annotations_category_count,
        max_labels=max_labels,
        url_base=url_base,
    ):
        self.test_datasets = datasets
        self.annotations_category_count = annotations_category_count
        self.max_labels = max_labels
        self.url_base = url_base
        self.test_notes = self.create_info_dict()

    def set_cell_count(self, dataset_name):
        dataset_url = self.test_datasets[dataset_name]["dataset_url"]
        headers = {"Content-Type": "application/octet-stream", "Cookie": cookie}
        response = self.client.get(f"{self.url_base}{dataset_url}/api/v0.2/schema", headers=headers)
        cell_count = json.loads(response._content)["schema"]["dataframe"]["nObs"]
        self.test_datasets[dataset_name]["cell_count"] = cell_count

    def create_info_dict(self):
        request_info = {}
        for dataset in self.test_datasets.keys():
            request_info[dataset] = {}
            for cat_count in self.annotations_category_count:
                request_info[dataset][f"num_categories_{cat_count}"] = {}
                for unique_labels in self.max_labels:
                    request_info[dataset][f"num_categories_{cat_count}"][f"max_label_{unique_labels}"] = {}
        return request_info

    def create_annotations_dict_multi_process(self, dataset_name, category_count, label_max):
        annotation_dict = {}
        futures = []
        categories = [f"Category{i}" for i in range(category_count)]
        if not self.test_datasets[dataset_name]["num_cells"]:
            self.set_cell_count(dataset_name)
        with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
            for category in categories:
                futures.append(
                    executor.submit(
                        self.build_array_for_category,
                        category,
                        self.test_datasets[dataset_name]["num_cells"],
                        label_max,
                    )
                )
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    category_name, cells = result
                    annotation_dict[category_name] = pd.Series(cells, dtype="category")
                except Exception as e:
                    print(f"Issue creating the annotations dict: {e}")
        return annotation_dict

    def build_array_for_category(self, category_name, cell_count, label_max):
        unique_label_count = label_max
        labels = self.generate_labels(unique_label_count)
        cells_per_label = int(cell_count / len(labels))
        extra = cell_count % len(labels)
        cells = []
        for label in labels:
            cells.extend([label] * cells_per_label)
        cells.extend(["extra"] * extra)
        rng = np.random.default_rng()
        rng.shuffle(cells)
        return category_name, cells

    @staticmethod
    def convert_to_fbs(annotation_dict):
        df = pd.DataFrame(annotation_dict)
        return encode_matrix_fbs(matrix=df, row_idx=None, col_idx=df.columns)

    @staticmethod
    def generate_labels(unique_label_count):
        labels = ["undefined"]
        for i in range(unique_label_count):
            length = random.randrange(10, 20)
            labels.append(f"{i}__" + "".join(random.choice(string.ascii_letters) for z in range(length)))
        return labels

    @contextmanager
    def elapsed_timer(self):
        start = default_timer()
        elapser = lambda: default_timer() - start  # noqa E731
        yield lambda: elapser()
        end = default_timer()
        elapser = lambda: end - start  # noqa E731

    def create_matrix(self, dataset_name, num_cat, max_labels):
        with self.elapsed_timer() as elapsed:
            annon_dict = self.create_annotations_dict_multi_process(dataset_name, num_cat, max_labels)
            dict_size = sum(sys.getsizeof(value) for value in annon_dict.values()) / 1024 ** 2
            self.test_notes[dataset_name][f"num_categories_{num_cat}"][f"max_label_{max_labels}"]["annotation_dict"] = {
                "creation_time": str(elapsed()),
                "size": f"{dict_size} mb",
            }
            df = pd.DataFrame(annon_dict)
            df_size = sys.getsizeof(df) / 1024 ** 2
            self.test_notes[dataset_name][f"num_categories_{num_cat}"][f"max_label_{max_labels}"]["data_frame"] = {
                "creation_time": str(elapsed()),
                "size": f"{df_size} mb",
            }
            try:
                matrix = encode_matrix_fbs(matrix=df, row_idx=None, col_idx=df.columns)
                matrix_size = sys.getsizeof(matrix) / 1024 ** 2
                self.test_notes[dataset_name][f"num_categories_{num_cat}"][f"max_label_{max_labels}"]["fbs_matrix"] = {
                    "creation_time": str(elapsed()),
                    "size": f"{matrix_size} mb",
                }
                return matrix
            except Exception as e:
                print(f"Issue creating fbs matrix: {e}, for {dataset_name}")
                return []

    def send_put_request(self, dataset_url, data):
        url = self.url_base + f"{dataset_url}/api/v0.2/annotations/obs"
        with self.elapsed_timer() as elapsed:
            try:
                headers = {"Content-Type": "application/octet-stream", "Cookie": cookie}
                response = requests.put(url=url, data=data, headers=headers)
            except Exception as e:
                print(f"Issue with put request: {e}")
                return None, elapsed()
            return response, elapsed()

    def test_categories_max_label_matrix(self, dataset_name):
        for unique_labels in self.max_labels:
            for category_count in self.annotations_category_count:
                print(f"Starting dataset: {dataset_name}, categories: {category_count}, labels: {unique_labels}")
                fbs_matrix = self.create_matrix(dataset_name, category_count, unique_labels)
                if self.test_datasets[dataset_name]["dataset_url"] and fbs_matrix:
                    response, response_time = self.send_put_request(
                        self.test_datasets[dataset_name]["dataset_url"], fbs_matrix
                    )
                    if response is None:
                        self.test_notes[dataset_name][f"num_categories_{category_count}"][f"max_label_{unique_labels}"][
                            "put_request"
                        ] = {"response_status": "failed", "request_time": str(response_time)}
                    else:
                        self.test_notes[dataset_name][f"num_categories_{category_count}"][f"max_label_{unique_labels}"][
                            "put_request"
                        ] = {"response_status": response.status_code, "request_time": str(response_time)}


def test_all_datasets():
    """
    Run time is dependent on number of datasets, dataset size, number of categories/number being tested and number of
     unique label counts being tested. However it generally takes a long time. I recommend running this in tmux
    """
    perf_test = PerformanceTestingAnnotations()
    for dataset_name in perf_test.test_datasets.keys():
        print(f"Testing annotation creation for: {dataset_name}")
        try:
            perf_test.test_categories_max_label_matrix(dataset_name)
        except Exception as e:
            print(f"something went wrong with {dataset_name}: {e}")
    return perf_test.test_notes


def main():
    notes = test_all_datasets()
    print(notes)


if __name__ == "__main__":
    main()
