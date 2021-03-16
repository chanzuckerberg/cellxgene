import sys
import argparse
import yaml

from backend.czi_hosted.common.config.app_config import AppConfig

def main():
    parser = argparse.ArgumentParser("A script to check hosted configuration files")
    parser.add_argument("config_file", help="the configuration file")
    parser.add_argument(
        "-s",
        "--show",
        default=False,
        action="store_true",
        help="print the configuration. NOTE: this may print secret values to stdout",
    )

    args = parser.parse_args()

    app_config = AppConfig()
    try:
        app_config.update_from_config_file(args.config_file)
        app_config.complete_config()
    except Exception as e:
        print(f"Error: {str(e)}")
        print("FAIL:", args.config_file)
        sys.exit(1)

    if args.show:
        yaml_config = app_config.config_to_dict()
        yaml.dump(yaml_config, sys.stdout)

    print("PASS:", args.config_file)
    sys.exit(0)


if __name__ == "__main__":
    main()
