import threading
from collections.abc import MutableMapping


class ImmutableKVCache(MutableMapping):
    """
    Guarantees that the factory will be called for each key once, and
    only once.
    """

    def __init__(self, factory):
        self.factory = factory
        self.lock = threading.Lock()
        self.cache = {}
        super().__init__()

    def __getitem__(self, key):
        if key not in self.cache:
            with self.lock:
                if key not in self.cache:
                    self.cache[key] = self.factory(key)
        return self.cache[key]

    def __iter__(self):
        """ weak iter, don't call factory """
        return self.cache.__iter__()

    def __len__(self):
        return self.cache.__len__()

    def __contains__(self, key):
        """ weak contain - don't call factory """
        return self.cache.__contains__(key)

    def __delitem__(self, key):
        with self.lock:
            del self.cache[key]

    def __setitem__(self, key, value):
        """ unsupported """
        raise NotImplementedError
