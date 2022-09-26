from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

res1 = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(0.0099999, 10, 0.010), rs.Mirror(0.300, 0, 0.595))
print(res1.waist_search(1064 * 10 ** -9))

res2 = rs.Resonator(2, rs.Mirror(0, 0, 0.005), rs.Mirror(0.00999956, 0, 0.005))
print(res2.waist_search(1064 * 10 ** -9))
