import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.cluster.hierarchy as sch
from radialtree import radialTree, plot


def _test_1(Z2, numleaf=50):
    # optionally leaves can be labeled by colors
    type_num = 12
    _cmp = cm.get_cmap("bwr", type_num)
    _cmp2 = cm.get_cmap("hot", type_num)
    colors_dict = {
        "example_color": _cmp(np.random.rand(numleaf)),
        "example_color2": _cmp2(np.random.rand(numleaf)),
    }
    colors_legends = {
        "example_color": {
            "colors": _cmp(np.linspace(0, 1, type_num)),
            "labels": ["ex1_" + str(i + 1) for i in range(type_num)],
        },
        "example_color2": {
            "colors": _cmp2(np.linspace(0, 1, type_num)),
            "labels": ["ex2_" + str(i + 1) for i in range(type_num)],
        },
    }
    # fig = pylab.figure(figsize=(8,8))

    # Compute and plot the dendrogram.
    # ax2 = fig.add_axes([0.3,0.71,0.6,0.2])

    fig, ax = plt.subplots(figsize=(10, 5))
    # plot(Z2, colorlabels=colors_dict,colorlabels_legend=colors_legends,show=True)
    radialTree(Z2, ax=ax, colorlabels=colors_dict, colorlabels_legend=colors_legends)
    fig.show()


def _test_2(Z2, numleaf=50):
    type_num = 6
    type_list = ["ex" + str(i) for i in range(type_num)]
    sample_classes = {
        "example_color": [np.random.choice(type_list) for i in range(numleaf)]
    }
    fig, ax = plt.subplots(figsize=(10, 5))
    radialTree(Z2, ax=ax, sample_classes=sample_classes)
    fig.show()
    # plot(Z2, sample_classes=sample_classes,show=True)


def _test_3(Z2, numleaf=50):
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    ax = ax.flatten()
    # no arguments
    radialTree(Z2, ax=ax[0])
    ax[0].set_aspect(1)

    type_num = 12
    _cmp = cm.get_cmap("bwr", type_num)
    _cmp2 = cm.get_cmap("hot", type_num)
    colors_dict = {
        "example_color": _cmp(np.random.rand(numleaf)),
        "example_color2": _cmp2(np.random.rand(numleaf)),
    }
    colors_legends = {
        "example_color": {
            "colors": _cmp(np.linspace(0, 1, type_num)),
            "labels": ["ex1_" + str(i + 1) for i in range(type_num)],
        },
        "example_color2": {
            "colors": _cmp2(np.linspace(0, 1, type_num)),
            "labels": ["ex2_" + str(i + 1) for i in range(type_num)],
        },
    }
    # fig = pylab.figure(figsize=(8,8))

    # Compute and plot the dendrogram.
    # ax2 = fig.add_axes([0.3,0.71,0.6,0.2])

    # like in test_1
    radialTree(Z2, ax=ax[1], colorlabels=colors_dict, colorlabels_legend=colors_legends)

    type_num = 6
    type_list = ["ex" + str(i) for i in range(type_num)]
    sample_classes = {
        "example_color": [np.random.choice(type_list) for i in range(numleaf)]
    }
    radialTree(Z2, ax=ax[2], sample_classes=sample_classes)
    ax[3].axis("off")
    fig.show()


if __name__ == "__main__":
    # Generate random features and distance matrix.

    test = [0, 1, 2, 3]
    np.random.seed(1)
    numleaf = 50
    _alphabets = [chr(i) for i in range(97, 97 + 24)]
    labels = sorted(
        ["".join(list(np.random.choice(_alphabets, 10))) for i in range(numleaf)]
    )

    x = np.random.rand(numleaf)
    D = np.zeros([numleaf, numleaf])
    for i in range(numleaf):
        for j in range(numleaf):

            D[i, j] = abs(x[i] - x[j])
    Y = sch.linkage(D, method="single")
    Z2 = sch.dendrogram(Y, labels=labels, no_plot=True)
    if 3 in test:
        _test_3(Z2)

    if 0 in test:
        plot(Z2, show=True)

    if 1 in test:
        _test_1(Z2)

    if 2 in test:
        _test_2(Z2)

    plt.show()
