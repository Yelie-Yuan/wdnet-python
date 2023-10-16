# user accessible functions
from ._hello import hello
from .wdnet_class import WDNet

import wdnet.rewire as rewire

del _hello
del wdnet_class

__all__ = ["WDNet", "hello", "rewire"]
