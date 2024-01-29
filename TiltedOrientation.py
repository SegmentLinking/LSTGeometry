#!/bin/env python

from tqdm import tqdm

class TiltedOrientation:
    def __init__(self, filepath):
        self.filepath = filepath
        self.data = {}
        self.f = open(self.filepath)
        for line in tqdm(self.f.readlines(), desc="Loading tilted orientation data"):
            if "#" == line[0]:
                continue
            ls = line.split()
            self.data[int(ls[0])] = float(ls[1])

    def getDrDz(self, detid):
        return self.data[detid]
