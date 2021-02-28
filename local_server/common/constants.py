from enum import Enum


class AugmentedEnum(Enum):
    def __hash__(self):
        return self.value.__hash__()

    def __eq__(self, other):
        if isinstance(other, type(self)) or isinstance(other, str):
            return self.value == other
        return False

    def __str__(self) -> str:
        return self.value


class Axis(AugmentedEnum):
    OBS = "obs"
    VAR = "var"


class DiffExpMode(AugmentedEnum):
    TOP_N = "topN"
    VAR_FILTER = "varFilter"


JSON_NaN_to_num_warning_msg = "JSON encoding failure - please verify all data are finite values (no NaN or Infinities)"
REACTIVE_LIMIT = 1_000_000

MAX_LAYOUTS = 30
