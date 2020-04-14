import sys
import argparse
import random
import time
import numpy as np

import server.compute.diffexp_cxg as diffexp_cxg
import server.compute.diffexp_generic as diffexp_generic

from server.common.app_config import AppConfig
from server.data_common.matrix_loader import MatrixDataLoader
from server.data_cxg.cxg_adaptor import CxgAdaptor


def main():
    parser = argparse.ArgumentParser("A command to test diffexp")
    parser.add_argument("dataset", help="name of a dataset to load")
    parser.add_argument("-na", "--numA", type=int, required=True, help="number of rows in group A")
    parser.add_argument("-nb", "--numB", type=int, required=True, help="number of rows in group B")
    parser.add_argument("-t", "--trials", default=1, type=int, help="number of trials")
    parser.add_argument(
        "-a", "--alg", choices=("default", "generic", "cxg"), default="default", help="algorithm to use"
    )
    parser.add_argument("-s", "--show", default=False, action="store_true", help="show the results")
    parser.add_argument(
        "-n", "--new-selection", default=False, action="store_true", help="change the selection between each trial"
    )
    parser.add_argument("--seed", default=1, type=int, help="set the random seed")

    args = parser.parse_args()

    app_config = AppConfig()
    app_config.single_dataset__datapath = args.dataset
    app_config.server__verbose = True
    app_config.complete_config()

    loader = MatrixDataLoader(args.dataset)
    adaptor = loader.open(app_config)

    if args.show:
        if isinstance(adaptor, CxgAdaptor):
            adaptor.open_array("X").schema.dump()

    numA = args.numA
    numB = args.numB
    rows = adaptor.get_shape()[0]

    random.seed(args.seed)

    if not args.new_selection:
        samples = random.sample(range(rows), numA + numB)
        filterA = samples[:numA]
        filterB = samples[numA:]

    for i in range(args.trials):
        if args.new_selection:
            samples = random.sample(range(rows), numA + numB)
            filterA = samples[:numA]
            filterB = samples[numA:]

        maskA = np.zeros(rows, dtype=bool)
        maskA[filterA] = True
        maskB = np.zeros(rows, dtype=bool)
        maskB[filterB] = True

        t1 = time.time()
        if args.alg == "default":
            results = adaptor.compute_diffexp_ttest(maskA, maskB)
        elif args.alg == "generic":
            results = diffexp_generic.diffexp_ttest(adaptor, maskA, maskB)
        elif args.alg == "cxg":
            if not isinstance(adaptor, CxgAdaptor):
                print("cxg only works with CxgAdaptor")
                sys.exit(1)
            results = diffexp_cxg.diffexp_ttest(adaptor, maskA, maskB)

        t2 = time.time()
        print("TIME=", t2 - t1)

    if args.show:
        for res in results:
            print(res)


if __name__ == "__main__":
    main()
