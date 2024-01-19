import threading
from collections.abc import MutableMapping


class ImmutableKVCache(MutableMapping):
    """
    Guarantees that the factory will be called for each key once, and
    only once.
    """

    def __init__(self, factory):
        self.factory = factory  # user-provided factory function
        self.lock = threading.Lock()  # guards factory_calls
        self.factory_calls = {}  # per-key factory condition variables
        self.cache = {}  # result cache, indexed by key
        super().__init__()

    def __getitem__(self, key):
        if key in self.cache:
            return self.cache[key]

        # we need to call factory.  First grab the main lock and the per-key CV.
        factory_calls = None
        creation_thr = False
        with self.lock:
            if key in self.cache:
                return self.cache[key]
            if key not in self.factory_calls:
                creation_thr = True
                self.factory_calls[key] = {"cv": threading.Condition(), "is_done": False, "error": None}
            factory_calls = self.factory_calls[key]

        # with the CV, create the value (or wait for it to be created)
        cv = factory_calls["cv"]
        with cv:
            if creation_thr:
                try:
                    self.cache[key] = self.factory(key)
                except Exception as e:
                    factory_calls["error"] = e

                factory_calls["is_done"] = True
                cv.notify_all()
            else:
                """wait for the value to be available"""
                while not factory_calls["is_done"]:
                    cv.wait()

        with self.lock:
            if key in self.factory_calls:
                del self.factory_calls[key]

        return self.cache[key]

    def __iter__(self):
        """weak iter, don't call factory"""
        return self.cache.__iter__()

    def __len__(self):
        return self.cache.__len__()

    def __contains__(self, key):
        """weak contain - don't call factory"""
        return self.cache.__contains__(key)

    def __delitem__(self, key):
        del self.cache[key]

    def __setitem__(self, key, value):
        """unsupported"""
        raise NotImplementedError
