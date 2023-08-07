"""
This module contains an implemeentation of Phenotype Phase Plane Analysis (PhPP).
(Edwards et al. 2001, Characterizing the metabolic phenotype: A phenotype phase plane analysis)
Adapted from code made by Kai Zhuang from the python package framed, as well as cobrapy v0.5.10

November 2021, Vetle Simensen
"""

from builtins import object
import numpy
import matplotlib.pyplot as plt
from palettable.colorbrewer import get_map

class PhenotypePhasePlane(object):
    def __init__(self, rxn_x, rxn_y, rxn_x_range, rxn_y_range):
        self.rxn_x = rxn_x
        self.rxn_y = rxn_y

        # converting reaction ranges to numpy array and storing it inside self
        self.x_range = rxn_x_range
        self.y_range = rxn_y_range

        # find length of reaction ranges
        len_x = len(self.x_range)
        len_y = len(self.y_range)

        # creating empty arrays for storing analysis results
        self.f_objective = numpy.zeros((len_x, len_y))
        self.shadow_price_x = numpy.zeros((len_x, len_y))
        self.shadow_price_y = numpy.zeros((len_x, len_y))
        self.segments = numpy.zeros(self.f_objective.shape, dtype=numpy.int32)
        self.phases = []

    def plot_PhPP(self, theme='Dark2', new_figure=True, show_plot=True):
        """
        theme: color theme used for distinguishing the different phases.
        new_figure: if set to True, a new matplotlib figure will be created.
        show_plot: if set to True, current figure will be shown
        """

        if new_figure:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")

        f = self.f_objective
        x = self.x_range
        y = self.y_range

        if 'EX_' in self.rxn_x:
            x = x * -1

        if 'EX_' in self.rxn_y:
            y = y * -1

        xgrid, ygrid = numpy.meshgrid(x, y)
        xgrid = xgrid.transpose()
        ygrid = ygrid.transpose()
        xgrid_scale, ygrid_scale = (1, 1)

        # Define theme colors
        colors = numpy.empty(self.f_objective.shape, dtype=numpy.dtype((str, 7)))
        n_segments = self.segments.max()
        color_list = get_map(theme, 'Qualitative', n_segments).hex_colors
        for i in range(n_segments):
            colors[self.segments == (i + 1)] = color_list[i]

        # Make surface plots, and add wireframe and axis labels
        ax.plot_surface(xgrid, ygrid, f, facecolors=colors, rstride=1, cstride=1, linewidth=0, antialiased=False)
        ax.plot_wireframe(xgrid, ygrid, f, color="black", rstride=xgrid_scale, cstride=ygrid_scale)
        ax.set_xlabel(self.rxn_x)
        ax.set_ylabel(self.rxn_y)
        ax.set_zlabel("Growth rate")
        ax.view_init(elev=30, azim=-135)
        fig.set_tight_layout(True)

        if show_plot:
            plt.show()

    def plot_shadow_price_x(self, new_figure=True, show_plot=True):
        """
        this method plots the shadow price of metabolites that are associated with reaction x
        new_figure: if set to True, a new matplotlib figure will be created.
        show_plot: if set to True, current figure will be shown
        """

        if new_figure:
            fig = plt.figure()
            ax = fig.add_subplot()
        sp_x = self.shadow_price_x
        x = self.x_range
        y = self.y_range

        if 'EX_' in self.rxn_x:
            x = x * -1

        if 'EX_' in self.rxn_y:
            y = y * -1

        plt.pcolormesh(x, y, numpy.transpose(sp_x))
        plt.colorbar()
        ax.set_xlabel(self.rxn_x)
        ax.set_ylabel(self.rxn_y)

        if show_plot:
            plt.show()

    def plot_shadow_price_y(self, new_figure=True, show_plot=True):
        """
        this method plots the shadow price of metabolites that are associated with reaction x
        new_figure: if set to True, a new matplotlib figure will be created.
        show_plot: if set to True, current figure will be shown
        """

        if new_figure:
            fig = plt.figure()
            ax = fig.add_subplot()
        sp_y = self.shadow_price_y
        x = self.x_range
        y = self.y_range

        if 'EX_' in self.rxn_x:
            x = x * -1

        if 'EX_' in self.rxn_y:
            y = y * -1

        plt.pcolormesh(x, y, numpy.transpose(sp_y))
        plt.colorbar()
        ax.set_xlabel(self.rxn_x)
        ax.set_ylabel(self.rxn_y)

        if show_plot:
            plt.show()

    def segment(self, threshold=0.001):
        """
        this method attempts to segment the data and identify the various phases
        of the phenotype phase plane
        """
        self.segments *= 0
        # each entry in phases will consist of the following tuple
        # ((x, y), shadow_price1, shadow_price2)
        self.phases = []
        segment_id = 0

        while self.segments.min() == 0:
            segment_id += 1
            # i and j are indices for a current point which has not been
            # assigned a segment yet
            i, j = numpy.unravel_index(self.segments.argmin(), self.segments.shape)
            # update the segment id for any point with a similar shadow price
            # to the current point
            d1 = abs(self.shadow_price_x - self.shadow_price_x[i, j])
            d2 = abs(self.shadow_price_y - self.shadow_price_y[i, j])
            self.segments[(d1 < threshold) * (d2 < threshold)] += segment_id
            # add the current point as one of the phases
            self.phases.append(((self.x_range[i], self.y_range[j]), self.shadow_price_x[i, j], self.shadow_price_y[i, j]))

def PhPP(model, rxn_x, rxn_y, rxn_x_range, rxn_y_range, target=None, maximize=True):
    """
    Phenotype Phase Plane Analysis
    analyze the changes in the objective function and the shadow prices
    Arguments:
        model (CBModel): the metabolic model
        rxn_x (str): reaction id to be plotted along x axis.  must be of a type convertable to numpy.array
        rxn_y (str): reaction id to be plotted along y axis.  must be of a type convertable to numpy.array
        rxn_x_range (list or array): the range of the reaction x
        rxn_y_range (list or array): the range of the reaction y
        target (str): the  reaction id of the optimization target.
                       if None is included, it will attempt to detect the biomass function
        maximize: True or False. the sense of the optimization

    Returns:
        phaseplane
    """

    # save initial flux bounds (avoids issues when lb > ub)
    old_lb_x = model.reactions.get_by_id(rxn_x).lower_bound
    old_ub_x = model.reactions.get_by_id(rxn_x).upper_bound
    old_lb_y = model.reactions.get_by_id(rxn_y).lower_bound
    old_ub_y = model.reactions.get_by_id(rxn_y).upper_bound

    # find metabolite ids corresponding to reactions x and y
    met_x = list(model.reactions.get_by_id(rxn_x).metabolites)[0].id
    met_y = list(model.reactions.get_by_id(rxn_y).metabolites)[0].id

    # create a PhenotypePhasePlane instance for storing results
    phase_plane = PhenotypePhasePlane(rxn_x, rxn_y, rxn_x_range, rxn_y_range)

    for i, v_x in enumerate(rxn_x_range):
        model.reactions.get_by_id(rxn_x).lower_bound = v_x
        model.reactions.get_by_id(rxn_x).upper_bound = v_x

        for j, v_y in enumerate(rxn_y_range):
            model.reactions.get_by_id(rxn_y).lower_bound = v_y
            model.reactions.get_by_id(rxn_y).upper_bound = v_y

            try:
                solution = model.optimize()
                if solution.status == 'optimal':
                    phase_plane.f_objective[i, j] = solution.objective_value
                    phase_plane.shadow_price_x[i, j] = solution.shadow_prices[met_x]
                    phase_plane.shadow_price_y[i, j] = solution.shadow_prices[met_y]
            except UserWarning:
                pass

        model.reactions.get_by_id(rxn_y).upper_bound = 0
        model.reactions.get_by_id(rxn_y).lower_bound = 0

    phase_plane.segment()

    # set bounds back to initial values
    model.reactions.get_by_id(rxn_x).upper_bound = old_ub_x
    model.reactions.get_by_id(rxn_x).lower_bound = old_lb_x
    model.reactions.get_by_id(rxn_y).upper_bound = old_ub_y
    model.reactions.get_by_id(rxn_y).lower_bound = old_lb_y

    return phase_plane
