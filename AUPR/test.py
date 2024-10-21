from AUPR import AUPR

expected_abd = [0, 0, 0.5, 0.5]
y_scores = [0.1, 0.4, 0.35, 0.8]
r = AUPR(expected_abd, y_scores)
print(r)
