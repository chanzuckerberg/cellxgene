import sys
import argparse
import random
import time
import numpy as np

import local_server.compute.diffexp_cxg as diffexp_cxg
import local_server.compute.diffexp_generic as diffexp_generic

from local_server.common.config.app_config import AppConfig
from local_server.data_common.matrix_loader import MatrixDataLoader
from local_server.data_cxg.cxg_adaptor import CxgAdaptor


def main():
    parser = argparse.ArgumentParser("A command to test diffexp")
    parser.add_argument("dataset", help="name of a dataset to load")
    parser.add_argument("-na", "--numA", type=int, help="number of rows in group A")
    parser.add_argument("-nb", "--numB", type=int, help="number of rows in group B")
    parser.add_argument("-va", "--varA", help="obs variable:value to use for group A")
    parser.add_argument("-vb", "--varB", help="obs variable:value to use for group B")
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
    app_config.update_server_config(single_dataset__datapath=args.dataset)
    app_config.update_server_config(app__verbose=True)
    app_config.complete_config()

    loader = MatrixDataLoader(args.dataset)
    adaptor = loader.open(app_config)

    if args.show:
        if isinstance(adaptor, CxgAdaptor):
            adaptor.open_array("X").schema.dump()

    random.seed(args.seed)
    np.random.seed(args.seed)
    rows = adaptor.get_shape()[0]

    if args.numA:
        filterA = random.sample(range(rows), args.numA)
    elif args.varA:
        vname, vval = args.varA.split(":")
        filterA = get_filter_from_obs(adaptor, vname, vval)
    else:
        print("must supply numA or varA")
        sys.exit(1)

    if args.numB:
        filterB = random.sample(range(rows), args.numB)
    elif args.varB:
        vname, vval = args.varB.split(":")
        filterB = get_filter_from_obs(adaptor, vname, vval)
    else:
        print("must supply numB or varB")
        sys.exit(1)

    for i in range(args.trials):
        if args.new_selection:
            if args.numA:
                filterA = random.sample(range(rows), args.numA)
            if args.numB:
                filterB = random.sample(range(rows), args.numB)

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


def get_filter_from_obs(adaptor, obsname, obsval):
    attrs = adaptor.get_obs_columns()
    if obsname not in attrs:
        print(f"Unknown obs attr {obsname}: expected on of {attrs}")
        sys.exit(1)
    obsvals = adaptor.query_obs_array(obsname)[:]
    obsval = type(obsvals[0])(obsval)

    vfilter = np.where(obsvals == obsval)[0]
    if len(vfilter) == 0:
        u = np.unique(obsvals)
        print(f"Unknown value in variable {obsname}:{obsval}: expected one of {list(u)}")
        sys.exit(1)

    return vfilter


if __name__ == "__main__":
    main()
