from robowski.settings import *
import importlib
import numpy as np
import robowski.smoothness_theory.smoothness_v2 as smoothness

def test_right_hand_side_speedup():
    smoothness_vs = smoothness.ReactionNetwork(N=5)
    for i in range(smoothness_vs.N):
        for j in range(smoothness_vs.N):
            for k in range(smoothness_vs.N):
                if j <= k:
                    smoothness_vs.k_tensor[i, j, k] = np.random.uniform(0, 1)
                else:
                    smoothness_vs.k_tensor[i, j, k] = 0

    y = np.random.uniform(0, 1, smoothness_vs.N)
    t = np.linspace(0, 1, 100)
    fast_result = smoothness_vs.right_hand_side_of_ode(t, y)
    slow_result = smoothness_vs.right_hand_side_of_ode_slow(t, y)
    # compare with isclose
    assert np.isclose(fast_result, slow_result).all()