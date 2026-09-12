"""Plume inflow/outflow model component: standalone functions supporting
IceCauldron.get_plume_fluxes.

Translated from MATLAB source: @IceCauldron/getPlumeFluxes.m's embedded local
functions (predictForest, predictTree) - the same implementation also exists,
after several earlier commented-out iterations, in
physics_sandbox_2023/predictForest.m.

TRANSLATION NOTE: per project direction, the plume component's random-forest
emulator path is explicitly an in-development/test mode, not finalized
physics, and is known to need a revamp. It also cannot be exercised
end-to-end in this repo: the emulator JSON file IceCauldron.plume_emulator_json
defaults to ("hydroplume_emulator_2025-04-22.json") is not present here, so
the exact structure of the decoded JSON (a bare list of trees vs. a wrapper
dict, and the tree node key names) is assumed to match the MATLAB code's
field names ('feature', 'threshold', 'left', 'right', 'value') rather than
independently verified.
"""

import numpy as np


def predict_tree(tree, input_vec):
    """Recursively walk a single decision tree (dict, decoded from JSON) to
    predict an output vector for one input sample.
    """
    if "feature" in tree and "threshold" in tree:
        # TRANSLATION NOTE: MATLAB does `featureIndex = tree.feature + 1` to adjust the JSON's
        # 0-based feature index for MATLAB's 1-based array indexing. Python arrays are already
        # 0-based like the JSON, so the index is used directly here, with no +1 - this is an
        # indexing-convention conversion (per CLAUDE.md's translation notes on 0- vs 1-based
        # indexing), not a logic change.
        feature_index = tree["feature"]
        threshold = tree["threshold"]

        if input_vec[feature_index] <= threshold:
            return predict_tree(tree["left"], input_vec)
        else:
            return predict_tree(tree["right"], input_vec)
    else:
        # TRANSLATION NOTE: MATLAB's iscolumn/isrow/reshape branch exists to force MATLAB's
        # row/column-vector distinction into a consistent row-vector shape - NumPy 1-D arrays
        # have no such distinction, so this collapses to a plain flatten.
        return np.ravel(tree["value"])


def predict_forest(forest, X):
    """finalPredictions = predictForest(forest, X)
    Average predictions across all trees in a random-forest emulator.

    forest = list of trees (dicts), decoded from the emulator JSON
    X      = (n_samples x n_features) input array
    """
    X = np.atleast_2d(X)
    num_trees = len(forest)
    num_samples = X.shape[0]

    # Get first prediction to determine output dimension
    first_pred = predict_tree(forest[0], X[0, :])
    output_dim = len(first_pred)

    predictions = np.zeros((num_samples, output_dim, num_trees))

    for t in range(num_trees):
        for s in range(num_samples):
            predictions[s, :, t] = predict_tree(forest[t], X[s, :])

    # Average across trees
    return predictions.mean(axis=2)  # size: (n_samples, output_dim)
