import pandas as pd

def join(dataFrames, names, truth):
    assert(len(dataFrames) == len(names))
    j = truth
    for df, n in zip(dataFrames, names):
        j = j.join(df, rsuffix='_{}'.format(n))
    return j
