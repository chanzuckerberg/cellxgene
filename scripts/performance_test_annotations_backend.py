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
import multiprocessing as mp
from timeit import default_timer
import concurrent.futures
import numpy as np

from server.compute.diffexp_cxg import get_thread_executor
from server.data_common.fbs.matrix import encode_matrix_fbs
from server.test import make_fbs
import pandas as pd
import random

users = ["cellxgene-smoke-test+staging@chanzuckerberg.com", "cellxgene-smoke-test+dev@chanzuckerberg.com"]
password = "Test1111"

# number of cells
datasets = {
    "33": {"url": "", "name": "", "num_cells": 33},
    "10k": {"url": "", "name": "", "num_cells": 9409},
    "100k": {"url": "", "name": "", "num_cells": 143015},
    "1million": {"url": "", "name": "", "num_cells": 1000000},
    "4million": {"url": "", "name": "", "num_cells": 4000000}
}
annotations_category_count = [1, 10, 100]


def signin_to_cellxgene(username, password):
    return {
        "Cookie": "cxguser=eyJhY2Nlc3NfdG9rZW4iOiAiY1lDcDBhQTJ5aGlUVU14R2FIZkNYVExxeEpFNHZxZkoiLCAiaWRfdG9rZW4iOiAiZXlKaGJHY2lPaUpTVXpJMU5pSXNJblI1Y0NJNklrcFhWQ0lzSW10cFpDSTZJazVpUlVOS2VEbFRWRWh6Tm5FeVgxSm9ORVZGTFNKOS5leUp1YVdOcmJtRnRaU0k2SW1ObGJHeDRaMlZ1WlMxemJXOXJaUzEwWlhOMEsyUmxkaUlzSW01aGJXVWlPaUpqWld4c2VHZGxibVV0YzIxdmEyVXRkR1Z6ZEN0a1pYWkFZMmhoYm5wMVkydGxjbUpsY21jdVkyOXRJaXdpY0dsamRIVnlaU0k2SW1oMGRIQnpPaTh2Y3k1bmNtRjJZWFJoY2k1amIyMHZZWFpoZEdGeUx6ZGtNekUyWm1Sak1URTNOV1JqWVRJME9HTmtabUppWXpsbE5UbGlZalF6UDNNOU5EZ3dKbkk5Y0djbVpEMW9kSFJ3Y3lVelFTVXlSaVV5Um1Oa2JpNWhkWFJvTUM1amIyMGxNa1poZG1GMFlYSnpKVEpHWTJVdWNHNW5JaXdpZFhCa1lYUmxaRjloZENJNklqSXdNakF0TURrdE16QlVNVGc2TURZNk1UTXVNREk0V2lJc0ltVnRZV2xzSWpvaVkyVnNiSGhuWlc1bExYTnRiMnRsTFhSbGMzUXJaR1YyUUdOb1lXNTZkV05yWlhKaVpYSm5MbU52YlNJc0ltVnRZV2xzWDNabGNtbG1hV1ZrSWpwbVlXeHpaU3dpYVhOeklqb2lhSFIwY0hNNkx5OXNiMmRwYmk1alpXeHNlR2RsYm1VdVpHVjJMbk5wYm1kc1pTMWpaV3hzTG1ONmFTNTBaV05vYm05c2IyZDVMeUlzSW5OMVlpSTZJbUYxZEdnd2ZEVm1Oekk0TVRnd05EZGlNVGhpTURBM05tRTBZakV6WmlJc0ltRjFaQ0k2SW1NelR6VjJkelpzVkdWMVREVjFhRkozWWpWR2RVc3hTbmR4TWpKdFVqVlNJaXdpYVdGMElqb3hOakF4TkRnNU1UY3pMQ0psZUhBaU9qRTJNREUxTWpVeE56TXNJbTV2Ym1ObElqb2lTSE51WkVvNFVGZFpTRlJ5YTBkTVVWcDBhM0lpZlEuaEZCcDM3VnJUVi1vX3RMR0p1SW1ieURyRXVZU0tBSlBjT0dxdW9YWm5QZ1Q5UW9PX2RGZGdBZzVnQ2RtVTlUa2dQTkxSXzJsczlPX0lxX2ZBX2NMbTNpVW5rY0ZnOFBsUVByMmh6RzVlN1AzLUVmX3RxZVdOTGw5LW9RTWJwekUyUTlaMXJyaTgyQklJMndzcEFnVWNJVHRYcDdLOTloOTVoYmM5a3U3REI0ZEtYUHNmQ0M0ZDFBc1cwcGNtcXJyMElEeHhPOVZ2LVB6cVNNdm5PQ2t2WmRza3ctMEpHOElsTm9ZOXBBWkxNWG9wMTIwNkZKdU95WHBCVWxiTUxJOVE5aVNhbEQ1dzduU0VQUzhjUl9YTHZXTzhiQ1ZxaXRURTZsajhnRlR1Rlc3emxfY005WVpqbVRwbUcweE95Q3R4VXdfQ0NCMi1NUmx2dU9GN3NmOHRnIiwgInJlZnJlc2hfdG9rZW4iOiAiTTlKakdLSnV2eklPMjRIbW9wQkhreEJwem5jcmxTSU4zaVZCMWxZaWxQRFh4IiwgImV4cGlyZXNfYXQiOiAxNjAxNTc1NTczfQ==; session=eyJfYXV0aDBfYXV0aGxpYl9ub25jZV8iOiJIc25kSjhQV1lIVHJrR0xRWnRrciJ9.X3TJFQ.1aZy0qSeclvinL6PESQl6rNuoUY"}


def create_annotations_dict_threaded(dataset, category_count, label_max):
    annotation_dict = {}
    futures = []
    categories = [f"Category{i}" for i in range(category_count)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        for category in categories:
            print(f"staring a thread for {category}")
            futures.append(executor.submit(build_array_for_category, category, dataset["num_cells"], label_max))
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()

                category_name, cells = result
                print(f"got the results for {category_name}")
                annotation_dict[category_name] = pd.Series(cells, dtype="category")
            except Exception as e:
                print(F"OH NO SOMETHING WENT WRONG: {e}")

    return annotation_dict


def create_annotations_dict_process_pool(dataset, category_count, label_max):
    pool = mp.Pool(processes=4)
    annotation_dict = {}
    results = [pool.apply(build_array_for_category, args=(f"Category{x}", dataset["num_cells"], label_max)) for x in range(category_count)]
    for x in results:
        annotation_dict[x[0]] = pd.Series(x[1], dtype="category")
    return annotation_dict


def create_annotations_dict(dataset, category_count, label_max):
    annotation_dict = {}
    categories = [f"Category{i}" for i in range(category_count)]
    for category in categories:
        category_name, cells = build_array_for_category(category, dataset["num_cells"], label_max)
        annotation_dict[category_name] = pd.Series(cells, dtype="category")
    return annotation_dict


def build_array_for_category(category_name, cell_count, label_max):
    unique_label_count = random.randrange(1, label_max)
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


def run_with_timer(dataset, num_cat, max_labels):
    with elapsed_timer() as elapsed:
        hold = create_annotations_dict(dataset, num_cat, max_labels)
        print(f"dict creation took {elapsed()}")
        a = convert_to_fbs(hold)
        print(f"total time took {elapsed()}")


def run_with_timer_threaded(dataset, num_cat, max_labels):
    with elapsed_timer() as elapsed:
        hold = create_annotations_dict_threaded(dataset, num_cat, max_labels)
        print(f"dict creation took {elapsed()}")
        a = convert_to_fbs(hold)
        print(f"total time took {elapsed()}")


def run_with_timer_pool(dataset, num_cat, max_labels):
    with elapsed_timer() as elapsed:
        hold = create_annotations_dict_process_pool(dataset, num_cat, max_labels)
        print(f"dict creation took {elapsed()}")
        a = convert_to_fbs(hold)
        print(f"total time took {elapsed()}")
