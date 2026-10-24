"""Independent 90-digit second differences of Bernoulli cloglog log likelihood.

Usage: python glm_cloglog_curvature_reference.py /tmp/glm-newton-curvature.csv
Requires Python's standard-library Decimal exp/ln; no numerical dependency.
"""
import csv
import sys
from decimal import Decimal, localcontext

with localcontext() as context, open(sys.argv[1], "w", newline="") as target:
    context.prec = 90
    writer = csv.writer(target)
    writer.writerow(["eta", "y", "a", "curvature"])
    step = Decimal("1e-15")
    for eta in [-34, -30, -20, -12, -10, -9.211, -9.2104, -9, -5, -3, -1, 0, 1, 2, 3]:
        value = Decimal(str(eta))
        for outcome in [0, 1]:
            def loglik(predictor):
                rate = predictor.exp()
                return (1 - (-rate).exp()).ln() if outcome else -rate

            curvature = -(loglik(value + step) - 2 * loglik(value) + loglik(value - step)) / step**2
            for weight in [1, 3]:
                writer.writerow([eta, outcome, weight, str(weight * curvature)])
