import time
import random

from locust import HttpUser, between, task

random.seed(time.time())
"""
To run this script sign into cellxgene in the desired environment and grab the returned cookie, update the cookie 
variable below with your cookie and run the following command to see results in the terminal:
locust -f scripts/scale_test_annotations.py --headless -u 30 -r 10 --host https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/e/ --run-time 5m 2>&1 | tee locust_dev_stats.txt

Or if you want to use the locust gui run:
locust -f scripts/scale_test_annotations.py -u 30 -r 10 --host https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/e/

If you want to test staging you'll need to substitute staging for dev in the host url
To test prod you'll need to replace dev.single-cell.czi.technology with cziscience.com
If you'd like to test additional datasets you'll need to add them to the dataset_urls array

Todo @mdunitz update script to retrieve different annotation categories -- may need to create them to ensure the 
categorys are shared across datasets for a given user 
"""
# Example Cookie
cookie = "cxguser=eyJhY2Nlc3NfdG9rZW4iOiAiY1lDcDBhQTJ5aGlUVU14R2FIZkNYVExxeEpFNHZxZkoiLCAiaWRfdG9rZW4iOiAiZXlKaGJHY2lPaUpTVXpJMU5pSXNJblI1Y0NJNklrcFhWQ0lzSW10cFpDSTZJazVpUlVOS2VEbFRWRWh6Tm5FeVgxSm9ORVZGTFNKOS5leUp1YVdOcmJtRnRaU0k2SW1ObGJHeDRaMlZ1WlMxemJXOXJaUzEwWlhOMEsyUmxkaUlzSW01aGJXVWlPaUpqWld4c2VHZGxibVV0YzIxdmEyVXRkR1Z6ZEN0a1pYWkFZMmhoYm5wMVkydGxjbUpsY21jdVkyOXRJaXdpY0dsamRIVnlaU0k2SW1oMGRIQnpPaTh2Y3k1bmNtRjJZWFJoY2k1amIyMHZZWFpoZEdGeUx6ZGtNekUyWm1Sak1URTNOV1JqWVRJME9HTmtabUppWXpsbE5UbGlZalF6UDNNOU5EZ3dKbkk5Y0djbVpEMW9kSFJ3Y3lVelFTVXlSaVV5Um1Oa2JpNWhkWFJvTUM1amIyMGxNa1poZG1GMFlYSnpKVEpHWTJVdWNHNW5JaXdpZFhCa1lYUmxaRjloZENJNklqSXdNakF0TURrdE16QlVNVGc2TURZNk1UTXVNREk0V2lJc0ltVnRZV2xzSWpvaVkyVnNiSGhuWlc1bExYTnRiMnRsTFhSbGMzUXJaR1YyUUdOb1lXNTZkV05yWlhKaVpYSm5MbU52YlNJc0ltVnRZV2xzWDNabGNtbG1hV1ZrSWpwbVlXeHpaU3dpYVhOeklqb2lhSFIwY0hNNkx5OXNiMmRwYmk1alpXeHNlR2RsYm1VdVpHVjJMbk5wYm1kc1pTMWpaV3hzTG1ONmFTNTBaV05vYm05c2IyZDVMeUlzSW5OMVlpSTZJbUYxZEdnd2ZEVm1Oekk0TVRnd05EZGlNVGhpTURBM05tRTBZakV6WmlJc0ltRjFaQ0k2SW1NelR6VjJkelpzVkdWMVREVjFhRkozWWpWR2RVc3hTbmR4TWpKdFVqVlNJaXdpYVdGMElqb3hOakF4TkRnNU1UY3pMQ0psZUhBaU9qRTJNREUxTWpVeE56TXNJbTV2Ym1ObElqb2lTSE51WkVvNFVGZFpTRlJ5YTBkTVVWcDBhM0lpZlEuaEZCcDM3VnJUVi1vX3RMR0p1SW1ieURyRXVZU0tBSlBjT0dxdW9YWm5QZ1Q5UW9PX2RGZGdBZzVnQ2RtVTlUa2dQTkxSXzJsczlPX0lxX2ZBX2NMbTNpVW5rY0ZnOFBsUVByMmh6RzVlN1AzLUVmX3RxZVdOTGw5LW9RTWJwekUyUTlaMXJyaTgyQklJMndzcEFnVWNJVHRYcDdLOTloOTVoYmM5a3U3REI0ZEtYUHNmQ0M0ZDFBc1cwcGNtcXJyMElEeHhPOVZ2LVB6cVNNdm5PQ2t2WmRza3ctMEpHOElsTm9ZOXBBWkxNWG9wMTIwNkZKdU95WHBCVWxiTUxJOVE5aVNhbEQ1dzduU0VQUzhjUl9YTHZXTzhiQ1ZxaXRURTZsajhnRlR1Rlc3emxfY005WVpqbVRwbUcweE95Q3R4VXdfQ0NCMi1NUmx2dU9GN3NmOHRnIiwgInJlZnJlc2hfdG9rZW4iOiAiTTlKakdLSnV2eklPMjRIbW9wQkhreEJwem5jcmxTSU4zaVZCMWxZaWxQRFh4IiwgImV4cGlyZXNfYXQiOiAxNjAxNTc1NTczfQ==; session=eyJfYXV0aDBfYXV0aGxpYl9ub25jZV8iOiJIc25kSjhQV1lIVHJrR0xRWnRrciJ9.X3TJFQ.1aZy0qSeclvinL6PESQl6rNuoU"


class WebsiteUser(HttpUser):
    wait_time = between(1, 2)
    dataset_urls = [
        "human_cell_landscape.cxg",
        "Single_cell_drug_screening_a549-42-remixed.cxg",
        "kampmann_lab_human_AD_snRNAseq_EC_inhibitoryNeurons-53-remixed.cxg",
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
