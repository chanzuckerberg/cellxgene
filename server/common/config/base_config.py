import copy

from flatten_dict import flatten
from server.common.errors import ConfigurationError


class BaseConfig(object):
    """This class handles the mechanics of updating and checking attributes.
    Derived classes are expected to store the actual attributes"""

    def __init__(self, app_config, default_config, dictval_cases={}):
        # reference back to the app_config
        self.app_config = app_config
        # the complete set of attribute and their default values (unflattened)
        self.default_config = default_config
        # attributes where the value may be a dict (and therefore are not flattened)
        self.dictval_cases = dictval_cases
        # used to make sure every attribute value is checked
        self.attr_checked = {key_name: False for key_name in self.create_mapping(default_config).keys()}

    def create_mapping(self, config):
        """Create a mapping from attribute names to (location in the config tree, value)"""
        config_copy = copy.deepcopy(config)
        mapping = {}

        # special cases where the value could be a dict.
        # If its value is not None, the entry is added to the mapping, and not included
        # in the flattening below.
        for dictval_case in self.dictval_cases:
            cur = config_copy
            for part in dictval_case[:-1]:
                cur = cur.get(part, {})
            val = cur.get(dictval_case[-1])
            if val is not None:
                key = "__".join(dictval_case)
                mapping[key] = (dictval_case, val)
                del cur[dictval_case[-1]]

        flat_config = flatten(config_copy)
        for key, value in flat_config.items():
            # name of the attribute
            attr = "__".join(key)
            mapping[attr] = (key, value)

        return mapping

    def check_attr(self, attrname, vtype):
        val = getattr(self, attrname)
        if type(vtype) in (list, tuple):
            if type(val) not in vtype:
                tnames = ",".join([x.__name__ for x in vtype])
                raise ConfigurationError(
                    f"Invalid type for attribute: {attrname}, expected types ({tnames}), got {type(val).__name__}"
                )
        else:
            if type(val) != vtype:
                raise ConfigurationError(
                    f"Invalid type for attribute: {attrname}, "
                    f"expected type {vtype.__name__}, got {type(val).__name__}"
                )

        self.attr_checked[attrname] = True

    def check_config(self):
        mapping = self.create_mapping(self.default_config)
        for key in mapping.keys():
            if not self.attr_checked[key]:
                raise ConfigurationError(f"The attr '{key}' has not been checked")

    def update(self, **kw):
        for key, value in kw.items():
            if not hasattr(self, key):
                raise ConfigurationError(f"unknown config parameter {key}.")
            try:
                if type(value) == tuple:
                    # convert tuple values to list values
                    value = list(value)
                setattr(self, key, value)
            except KeyError:
                raise ConfigurationError(f"Unable to set config parameter {key}.")

            self.attr_checked[key] = False

    def update_from_config(self, config, prefix):
        mapping = self.create_mapping(config)
        for attr, (key, value) in mapping.items():
            if not hasattr(self, attr):
                raise ConfigurationError(f"Unknown key from config file: {prefix}__{attr}")
            try:
                setattr(self, attr, value)
            except KeyError:  # TODO ask brian when this would be raised (instead of being caught in attribute check above)
                raise ConfigurationError(f"Unable to set config attribute: {prefix}__{attr}")

            self.attr_checked[attr] = False

    def changes_from_default(self):
        """Return all the attribute that are different from the default"""
        mapping = self.create_mapping(self.default_config)
        diff = []
        for attrname, (key, defval) in mapping.items():
            curval = getattr(self, attrname)
            if curval != defval:
                diff.append((attrname, curval, defval))
        return diff
