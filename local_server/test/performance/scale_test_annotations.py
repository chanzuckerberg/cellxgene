import time
import random

from locust import HttpUser, between, task

random.seed(time.time())
"""
To run this script sign into cellxgene in the desired environment and grab the returned cookie, update the cookie
variable below with your cookie and run the following command to see results in the terminal:
locust -f local_server/test/performance/scale_test_annotations.py --headless -u 30 -r 10 --host https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/e/ --run-time 5m 2>&1 | tee locust_dev_stats.txt

Or if you want to use the locust gui run:
locust -f local_server/test/performance/scale_test_annotations.py -u 30 -r 10 --host https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/e/

If you want to test staging you'll need to substitute staging for dev in the host url
To test prod you'll need to replace dev.single-cell.czi.technology with cziscience.com
If you'd like to test additional datasets you'll need to add them to the dataset_urls array

Todo @mdunitz update script to retrieve different annotation categories -- may need to create them to ensure the
categories are shared across datasets for a given user.
"""
cookie = ""


class WebsiteUser(HttpUser):
    wait_time = between(1, 2)
    dataset_urls = [
        "human_cell_landscape.cxg",
        "Single_cell_drug_screening_a549-42-remixed.cxg",
        "krasnow_lab_human_lung_cell_atlas_smartseq2-2-remixed.cxg",
        "Single_cell_gene_expression_profiling_of_SARS_CoV_2_infected_human_cell_lines_H1299-27-remixed.cxg",
    ]

    @task
    def get_annotations(self):
        dataset_url = random.choice(self.dataset_urls)
        url = f"{dataset_url}/api/v0.2/annotations/obs?annotation-name=cell_type"
        headers = {"Content-Type": "application/octet-stream", "Cookie": cookie}
        self.client.get(url, headers=headers)

    @task
    def get_schema(self):
        dataset_url = random.choice(self.dataset_urls)
        headers = {"Content-Type": "application/octet-stream", "Cookie": cookie}
        self.client.get(f"{dataset_url}/api/v0.2/schema", headers=headers)
