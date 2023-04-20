import time
from contextlib import ContextDecorator


class timing(ContextDecorator):
    def __init__(self, name="block", enabled=True):
        self.name = name
        self.start = 0
        self.end = 0
        self.enabled = enabled

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, *exc):
        self.end = time.time()
        if self.enabled:
            print(f"{self.name} took {1000*(self.end-self.start)}ms")
