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

function compute_correlations(x, y, verbose::Bool)
    sp = py"spear"(x, y)
    if (verbose)
        println("    Spearman computed")
    end
    kt = py"ktau"(x, y)
    if (verbose)
        println("    Kendall tau computed")
    end
    wkt = py"wktau"(x, y)
    if (verbose)
        println("    Weighted Kendall tau computed")
    end
    return sp, kt, wkt
end
