import json

import pprint
pprint = pprint.PrettyPrinter(indent=1).pprint

class options():
    def __init__(self, fname):
        self.props = dict()
        self.fill(fname)

    def fill(self, fname):
        f = open(fname)
        self.props = json.load(f)

class universalOptions():
    def __init__(self):
        self.multip = 1
        self.bigNumber = 300
