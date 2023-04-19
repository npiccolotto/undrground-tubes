import time
from contextlib import ContextDecorator


class timing(ContextDecorator):
    def __init__(self, name="block"):
        self.name = name
        self.start = 0
        self.end = 0

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, *exc):
        self.end = time.time()
        print(f"{self.name} took {1000*(self.end-self.start)}ms")
