using Conda
using PyCall

Conda.add("scipy")

py"""
from scipy import stats

def ktau(x, y):
    return stats.kendalltau(x, y)

def wktau(x, y):
    return stats.weightedtau(x, y)

def spear(x,y):
    return stats.spearmanr(x, y)
"""

function compute_correlations(x, y)
    sp = py"spear"(x, y)
    println("    Spearman computed")
    kt = py"ktau"(x, y)
    println("    Kendall tau computed")
    wkt = py"wktau"(x, y)
    println("    Weighted Kendall tau computed")
    return sp, kt, wkt
end
