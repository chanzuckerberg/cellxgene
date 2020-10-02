## Find data sets with
# - [ ] 10k
# - [ ] 100k
# - [ ] 1 millions cells
## Function to Generate X number of categories (with full array)
## Function to generate X different labels of Y length for each cell in the dataset
## Use locust to hit endpoint with different numbers of
# - [ ]users how many??  -- does it need to be with different users?
#  - [ ]how to do auth? -- need to sign in to get token? -- store tokens and use randomly -- how to pass the credentials to lotus?

## create scheduled github action
## clean up db/s3?


## Perf test creation

## Perf test retrieving annotations as well

# perf test get schema route
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

users = ["cellxgene-smoke-test+staging@chanzuckerberg.com", "cellxgene-smoke-test+dev@chanzuckerberg.com"]
password = "Test1111"

# number of cells
datasets = {
    "smallest": {"dataset_url": "kampmann_lab_human_AD_snRNAseq_EC_inhibitoryNeurons-53-remixed.cxg",
                 "name": "smallest", "num_cells": 5270},
    "10k": {
        "dataset_url": "krasnow_lab_human_lung_cell_atlas_smartseq2-2-remixed.cxg",
        "name": "10k",
        "num_cells": 9409
    },
    "80k": {"dataset_url": "Single_cell_gene_expression_profiling_of_SARS_CoV_2_infected_human_cell_lines_H1299-27-remixed.cxg",
            "name": "80k",
            "num_cells": 81736},
    "140k": {"dataset_url": "Single_cell_drug_screening_a549-42-remixed.cxg", "name": "140k",
             "num_cells": 143015},
    "largest": {"dataset_url": "human_cell_landscape.cxg", "name": "largest", "num_cells": 599926},
    "1million": {"dataset_url": None, "name": "1million", "num_cells": 1000000},

    "4million": {"dataset_url": None, "name": "4million", "num_cells": 4000000}
}
# annotations_category_count = [1, 10, 100]
# max_labels = [5, 50, 200]
annotations_category_count = [50, 75, 90]
max_labels = [5, 50, 200]

def signin_to_cellxgene(username, password):
    return {
        "Cookie": "cxguser=eyJhY2Nlc3NfdG9rZW4iOiAiY1lDcDBhQTJ5aGlUVU14R2FIZkNYVExxeEpFNHZxZkoiLCAiaWRfdG9rZW4iOiAiZXlKaGJHY2lPaUpTVXpJMU5pSXNJblI1Y0NJNklrcFhWQ0lzSW10cFpDSTZJazVpUlVOS2VEbFRWRWh6Tm5FeVgxSm9ORVZGTFNKOS5leUp1YVdOcmJtRnRaU0k2SW1ObGJHeDRaMlZ1WlMxemJXOXJaUzEwWlhOMEsyUmxkaUlzSW01aGJXVWlPaUpqWld4c2VHZGxibVV0YzIxdmEyVXRkR1Z6ZEN0a1pYWkFZMmhoYm5wMVkydGxjbUpsY21jdVkyOXRJaXdpY0dsamRIVnlaU0k2SW1oMGRIQnpPaTh2Y3k1bmNtRjJZWFJoY2k1amIyMHZZWFpoZEdGeUx6ZGtNekUyWm1Sak1URTNOV1JqWVRJME9HTmtabUppWXpsbE5UbGlZalF6UDNNOU5EZ3dKbkk5Y0djbVpEMW9kSFJ3Y3lVelFTVXlSaVV5Um1Oa2JpNWhkWFJvTUM1amIyMGxNa1poZG1GMFlYSnpKVEpHWTJVdWNHNW5JaXdpZFhCa1lYUmxaRjloZENJNklqSXdNakF0TURrdE16QlVNVGc2TURZNk1UTXVNREk0V2lJc0ltVnRZV2xzSWpvaVkyVnNiSGhuWlc1bExYTnRiMnRsTFhSbGMzUXJaR1YyUUdOb1lXNTZkV05yWlhKaVpYSm5MbU52YlNJc0ltVnRZV2xzWDNabGNtbG1hV1ZrSWpwbVlXeHpaU3dpYVhOeklqb2lhSFIwY0hNNkx5OXNiMmRwYmk1alpXeHNlR2RsYm1VdVpHVjJMbk5wYm1kc1pTMWpaV3hzTG1ONmFTNTBaV05vYm05c2IyZDVMeUlzSW5OMVlpSTZJbUYxZEdnd2ZEVm1Oekk0TVRnd05EZGlNVGhpTURBM05tRTBZakV6WmlJc0ltRjFaQ0k2SW1NelR6VjJkelpzVkdWMVREVjFhRkozWWpWR2RVc3hTbmR4TWpKdFVqVlNJaXdpYVdGMElqb3hOakF4TkRnNU1UY3pMQ0psZUhBaU9qRTJNREUxTWpVeE56TXNJbTV2Ym1ObElqb2lTSE51WkVvNFVGZFpTRlJ5YTBkTVVWcDBhM0lpZlEuaEZCcDM3VnJUVi1vX3RMR0p1SW1ieURyRXVZU0tBSlBjT0dxdW9YWm5QZ1Q5UW9PX2RGZGdBZzVnQ2RtVTlUa2dQTkxSXzJsczlPX0lxX2ZBX2NMbTNpVW5rY0ZnOFBsUVByMmh6RzVlN1AzLUVmX3RxZVdOTGw5LW9RTWJwekUyUTlaMXJyaTgyQklJMndzcEFnVWNJVHRYcDdLOTloOTVoYmM5a3U3REI0ZEtYUHNmQ0M0ZDFBc1cwcGNtcXJyMElEeHhPOVZ2LVB6cVNNdm5PQ2t2WmRza3ctMEpHOElsTm9ZOXBBWkxNWG9wMTIwNkZKdU95WHBCVWxiTUxJOVE5aVNhbEQ1dzduU0VQUzhjUl9YTHZXTzhiQ1ZxaXRURTZsajhnRlR1Rlc3emxfY005WVpqbVRwbUcweE95Q3R4VXdfQ0NCMi1NUmx2dU9GN3NmOHRnIiwgInJlZnJlc2hfdG9rZW4iOiAiTTlKakdLSnV2eklPMjRIbW9wQkhreEJwem5jcmxTSU4zaVZCMWxZaWxQRFh4IiwgImV4cGlyZXNfYXQiOiAxNjAxNTc1NTczfQ==; session=eyJfYXV0aDBfYXV0aGxpYl9ub25jZV8iOiJIc25kSjhQV1lIVHJrR0xRWnRrciJ9.X3TJFQ.1aZy0qSeclvinL6PESQl6rNuoUY"}


def create_annotations_dict_threaded(dataset, category_count, label_max):
    annotation_dict = {}
    futures = []
    categories = [f"Category{i}" for i in range(category_count)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
        for category in categories:
            futures.append(executor.submit(build_array_for_category, category, dataset["num_cells"], label_max))
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()

                category_name, cells = result
                annotation_dict[category_name] = pd.Series(cells, dtype="category")
            except Exception as e:
                print(F"OH NO SOMETHING WENT WRONG: {e}")

    return annotation_dict


def create_annotations_dict(dataset, category_count, label_max):
    annotation_dict = {}
    categories = [f"Category{i}" for i in range(category_count)]
    for category in categories:
        category_name, cells = build_array_for_category(category, dataset["num_cells"], label_max)
        annotation_dict[category_name] = pd.Series(cells, dtype="category")
    return annotation_dict


def build_array_for_category(category_name, cell_count, label_max):
    unique_label_count = label_max
    labels = generate_labels(unique_label_count)
    cells_per_label = int(cell_count / len(labels))
    extra = cell_count % len(labels)
    cells = []
    for label in labels:
        cells.extend([label] * cells_per_label)
    cells.extend(["extra"] * extra)
    rng = np.random.default_rng()
    rng.shuffle(cells)
    return category_name, cells


def convert_to_fbs(annotation_dict):
    df = pd.DataFrame(annotation_dict)
    return encode_matrix_fbs(matrix=df, row_idx=None, col_idx=df.columns)


def generate_labels(unique_label_count):
    labels = ["undefined"]
    for i in range(unique_label_count):
        length = random.randrange(10, 20)
        labels.append(''.join(random.choice(string.ascii_letters) for i in range(length)))
    return labels


@contextmanager
def elapsed_timer():
    start = default_timer()
    elapser = lambda: default_timer() - start
    yield lambda: elapser()
    end = default_timer()
    elapser = lambda: end - start


def run_with_timer_threaded(dataset, num_cat, max_labels):
    info = ""
    with elapsed_timer() as elapsed:
        info += f"dataset: {dataset['name']}, \nnum_cat: {num_cat}, max_labels: {max_labels}, num_cells: {dataset['num_cells']}"
        hold = create_annotations_dict_threaded(dataset, num_cat, max_labels)
        dict_size = sum(sys.getsizeof(value) for value in hold.values())
        info += f"\ndict creation took {elapsed()}, dict is {dict_size} bytes"
        df = pd.DataFrame(hold)
        import pdb
        pdb.set_trace()
        info += f"\ndataframe creation took {elapsed()}, df is {sys.getsizeof(df) / 1024 ** 2} mb"
        try:
            a = encode_matrix_fbs(matrix=df, row_idx=None, col_idx=df.columns)
            info +=f"\n matrix creation took {elapsed()}, matrix is {sys.getsizeof(a) / 1024 ** 2} mb"
            return a, info
        except Exception as e:
            return []




cookie = "cxguser=eyJhY2Nlc3NfdG9rZW4iOiAiY1lDcDBhQTJ5aGlUVU14R2FIZkNYVExxeEpFNHZxZkoiLCAiaWRfdG9rZW4iOiAiZXlKaGJHY2lPaUpTVXpJMU5pSXNJblI1Y0NJNklrcFhWQ0lzSW10cFpDSTZJazVpUlVOS2VEbFRWRWh6Tm5FeVgxSm9ORVZGTFNKOS5leUp1YVdOcmJtRnRaU0k2SW1ObGJHeDRaMlZ1WlMxemJXOXJaUzEwWlhOMEsyUmxkaUlzSW01aGJXVWlPaUpqWld4c2VHZGxibVV0YzIxdmEyVXRkR1Z6ZEN0a1pYWkFZMmhoYm5wMVkydGxjbUpsY21jdVkyOXRJaXdpY0dsamRIVnlaU0k2SW1oMGRIQnpPaTh2Y3k1bmNtRjJZWFJoY2k1amIyMHZZWFpoZEdGeUx6ZGtNekUyWm1Sak1URTNOV1JqWVRJME9HTmtabUppWXpsbE5UbGlZalF6UDNNOU5EZ3dKbkk5Y0djbVpEMW9kSFJ3Y3lVelFTVXlSaVV5Um1Oa2JpNWhkWFJvTUM1amIyMGxNa1poZG1GMFlYSnpKVEpHWTJVdWNHNW5JaXdpZFhCa1lYUmxaRjloZENJNklqSXdNakF0TURrdE16QlVNVGc2TURZNk1UTXVNREk0V2lJc0ltVnRZV2xzSWpvaVkyVnNiSGhuWlc1bExYTnRiMnRsTFhSbGMzUXJaR1YyUUdOb1lXNTZkV05yWlhKaVpYSm5MbU52YlNJc0ltVnRZV2xzWDNabGNtbG1hV1ZrSWpwbVlXeHpaU3dpYVhOeklqb2lhSFIwY0hNNkx5OXNiMmRwYmk1alpXeHNlR2RsYm1VdVpHVjJMbk5wYm1kc1pTMWpaV3hzTG1ONmFTNTBaV05vYm05c2IyZDVMeUlzSW5OMVlpSTZJbUYxZEdnd2ZEVm1Oekk0TVRnd05EZGlNVGhpTURBM05tRTBZakV6WmlJc0ltRjFaQ0k2SW1NelR6VjJkelpzVkdWMVREVjFhRkozWWpWR2RVc3hTbmR4TWpKdFVqVlNJaXdpYVdGMElqb3hOakF4TkRnNU1UY3pMQ0psZUhBaU9qRTJNREUxTWpVeE56TXNJbTV2Ym1ObElqb2lTSE51WkVvNFVGZFpTRlJ5YTBkTVVWcDBhM0lpZlEuaEZCcDM3VnJUVi1vX3RMR0p1SW1ieURyRXVZU0tBSlBjT0dxdW9YWm5QZ1Q5UW9PX2RGZGdBZzVnQ2RtVTlUa2dQTkxSXzJsczlPX0lxX2ZBX2NMbTNpVW5rY0ZnOFBsUVByMmh6RzVlN1AzLUVmX3RxZVdOTGw5LW9RTWJwekUyUTlaMXJyaTgyQklJMndzcEFnVWNJVHRYcDdLOTloOTVoYmM5a3U3REI0ZEtYUHNmQ0M0ZDFBc1cwcGNtcXJyMElEeHhPOVZ2LVB6cVNNdm5PQ2t2WmRza3ctMEpHOElsTm9ZOXBBWkxNWG9wMTIwNkZKdU95WHBCVWxiTUxJOVE5aVNhbEQ1dzduU0VQUzhjUl9YTHZXTzhiQ1ZxaXRURTZsajhnRlR1Rlc3emxfY005WVpqbVRwbUcweE95Q3R4VXdfQ0NCMi1NUmx2dU9GN3NmOHRnIiwgInJlZnJlc2hfdG9rZW4iOiAiTTlKakdLSnV2eklPMjRIbW9wQkhreEJwem5jcmxTSU4zaVZCMWxZaWxQRFh4IiwgImV4cGlyZXNfYXQiOiAxNjAxNTc1NTczfQ==; session=eyJfYXV0aDBfYXV0aGxpYl9ub25jZV8iOiJIc25kSjhQV1lIVHJrR0xRWnRrciJ9.X3TJFQ.1aZy0qSeclvinL6PESQl6rNuoUY"
def send_put_request(dataset_url, data, info):
    if data == []:
        info += "\n unable to create dataset"
    url = f"https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/e/{dataset_url}/api/v0.2/annotations/obs"
    with elapsed_timer() as elapsed:
        try:
            headers = {"Content-Type":"application/octet-stream", "Cookie": cookie}
            something = requests.put(url=url, data=data, headers=headers)
        except Exception as e:
            print(f"SOMETHING WENT WRONG: {e}")
            return
        info += f"\nrequest took this long {elapsed()}"
        info += f"\nRequest {something.status_code}"
        return info


def here_we_Go():
    for dataset in datasets.keys():
        print(f"DATASET: {dataset}")
        test_category_label_max(dataset)


def test_category_label_max(dataset):
    futures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        for count in annotations_category_count:
            for max in max_labels:
                matrix, info = run_with_timer_threaded(datasets[dataset], count, max)
                if datasets[dataset]["dataset_url"] is not None:
                    futures.append(executor.submit(send_put_request, datasets[dataset]["dataset_url"], matrix, info))
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                print(result)
            except Exception as e:
                info += f"\n SOMETHING WENT WRONG: {e}"
                print(info)




def create_array_and_put_request(dataset, count, max):
    hold = run_with_timer_threaded(datasets[dataset], count, max)
    if datasets[dataset]["dataset_url"] is not None:
        a = send_put_request(datasets[dataset]["dataset_url"], hold)
        if a.status_code == 504:
            return


def df_memory_usage(df):
    return(round(df.memory_usage(deep=True).sum() / 1024 ** 2, 2))


def get_annotations(dataset_url, cookie):
    url = f"https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/e/{dataset_url}/api/v0.2/annotations/obs"
    with elapsed_timer() as elapsed:
        try:
            headers = {"Content-Type": "application/octet-stream", "Cookie": cookie}
            something = requests.get(url=url, headers=headers)
        except Exception as e:
            print(f"SOMETHING WENT WRONG: {e}")
        return something


def get_schema(dataset_url, cookie):
    url = f"https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/e/{dataset_url}/api/v0.2/schema"
    with elapsed_timer() as elapsed:
        try:
            headers = {"Content-Type": "application/octet-stream", "Cookie": cookie}
            something = requests.get(url=url, headers=headers)
        except Exception as e:
            print(f"SOMETHING WENT WRONG: {e}")
        return something
