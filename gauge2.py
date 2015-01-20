#!/usr/bin/env python
"""
Speed-Gauge for matplotlib
See http://nbviewer.ipython.org/gist/nicolasfauchereau/794df533eca594565ab3
Adapted to be more typical ratio display.
"""

from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np

from matplotlib.patches import Circle, Wedge, Rectangle


def degree_range(n):
    start = np.linspace(0, 180, n+1, endpoint=True)[0:-1]
    end = np.linspace(0, 180, n+1, endpoint=True)[1::]
    mid_points = start + ((end-start)/2.)
    return np.c_[start, end], mid_points


def rot_text(ang):
    rotation = np.degrees(np.radians(ang) * np.pi / np.pi - np.radians(90))
    return rotation


def gauge(colors='hot', title='', value=0.5, labs=['0%', -1, '100%']):
    """
Gauge display of a single value against a range.

Parameters
----------

hot : string
    matplotlib colour map specification

title : string
    To put at the base of the gauge. If falsy, no label area.

value : float
    Fractional value to plot, 0 for no colour, 1 for whole bar

labs : 3-element list of tuple
    Label the start, arrow and end of the gauge with these. If the middle one
    is -1, replace with percentage of value.
    """

    N = 100  # number of elements

    # if colors is a colormap
    cmap = cm.get_cmap(colors, N)
    cmap = cmap(np.arange(N))
    colors = cmap[::-1, :].tolist()

    fig, ax = plt.subplots()
    ang_range, mid_points = degree_range(N)
    pos = mid_points[int(N*(1-value))]

    patches = []
    for ang, c in zip(ang_range, colors):
        # sectors
        # patches.append(Wedge((0., 0.), .4, *ang, facecolor='w', lw=2))
        # arcs
        if sum(ang) < pos*2:
            c = 'w'
        patches.append(Wedge((0., 0.), .4, *ang, width=0.10,
                             facecolor=c, lw=2, alpha=0.5, edgecolor='None'))

    for mid, lab in zip([180, pos, 0], labs):
        if lab == -1:
            lab = "{:.1%}".format(value)
        if lab:
            ax.text(0.41 * np.cos(np.radians(mid)), 0.41 * np.sin(np.radians(mid)),
                    lab, ha='center', va='center', fontsize=14, fontweight='bold',
                    rotation=rot_text(mid))

    [ax.add_patch(p) for p in patches]

    if title:
        r = Rectangle((-0.4, -0.1), 0.8, 0.1, facecolor='w', lw=2)
        ax.add_patch(r)
        ax.text(0, -0.05, title, horizontalalignment='center',
                verticalalignment='center', fontsize=22, fontweight='bold')

    ax.arrow(0, 0, 0.225 * np.cos(np.radians(pos)), 0.225 *
             np.sin(np.radians(pos)), width=0.04, head_width=0.09,
             head_length=0.1, fc='k', ec='k')

    ax.add_patch(Circle((0, 0), radius=0.02, facecolor='k'))
    ax.add_patch(Circle((0, 0), radius=0.01, facecolor='w', zorder=11))

    ax.set_frame_on(False)
    ax.axes.set_xticks([])
    ax.axes.set_yticks([])
    ax.axis('equal')
    plt.tight_layout()
