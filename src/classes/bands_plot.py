import sys
import os
from load_data import load_palette
import matplotlib.pyplot as plt
import matplotlib as mpl


class BandsPlot:
    """
    Class for plotting band structures with multiple segments.

    Attributes:
        Np (int): Number of band segments.
        width (list): Widths of each band segment.
        palette (dict): Color palette from the data.
        fig: Figure object for the plot.
        ax: List of Axes objects for each band segment. 
        x (list): x-coordinates for the band structure.
        energy (list): Band energies for the band structure.
    """

    def __init__(self, x_ticks, x_labels, y_range=[0.0, 0.0]):

        self.Np = len(x_ticks)
        self.width = [0.0] * self.Np
        for ip in range(self.Np):
            self.width[ip] = x_ticks[ip][-1]

        self.palette = load_palette()
        mpl.rcParams["font.sans-serif"].insert(0, "Noto Sans")
        mpl.rcParams.update({'font.size': 14})

        self.fig = plt.figure(figsize=(5, 3.75))
        gs0 = self.fig.add_gridspec(1, self.Np, wspace=0.0, hspace=0.0, left=0.14, right=0.98,
                                    top=0.97, bottom=0.07, width_ratios=self.width)

        self.ax = []
        for ip in range(self.Np):
            self.ax.append(self.fig.add_subplot(gs0[ip]))

        for ip in range(self.Np):
            self.ax[ip].grid(axis="x", linewidth=1, color=self.palette["gray"], zorder=0)
            self.ax[ip].axhline(linewidth=1, color=self.palette["gray"], zorder=0)
            self.ax[ip].set_ylim(y_range[0], y_range[1])
            self.ax[ip].set_xlim(x_ticks[ip][0], x_ticks[ip][-1])
            self.ax[ip].set_xticks(x_ticks[ip], x_labels[ip], color=self.palette["black"])
            self.ax[ip].tick_params(axis="x", direction="in", length=0)
            self.ax[ip].tick_params(axis="y", left=False, right=False, direction="in", width=1, zorder=0,
                                    color=self.palette["gray"], labelcolor=self.palette["black"])
            if ip != 0:
                self.ax[ip].yaxis.set_ticklabels([])
            for edge in ["bottom", "top", "left", "right"]:
                self.ax[ip].spines[edge].set_color(self.palette["black"])
                self.ax[ip].spines[edge].set_linewidth(1)
                self.ax[ip].spines[edge].set_zorder(4)
        self.ax[0].set_ylabel("Energy (eV)", labelpad=-2, color=self.palette["black"])
        self.ax[0].tick_params(axis="y", left=True)
        self.ax[-1].tick_params(axis="y", right=True)

    def __str__(self):
        str_out = "BandsPlot:\n"
        str_out += " width = " + str(self.width) + "\n"
        return str_out

    def add_plot(self, x, energy, color=None, label=None, zorder=3):

        if color is None:
            color = self.palette["darkblue"]

        ix_edges = [0] * (self.Np+1)
        for ip in range(self.Np):
            ix_edges[ip+1] = ix_edges[ip] + len(x[ip])

        Nb = len(energy)
        for ip in range(self.Np):
            for ib in range(Nb):
                self.ax[ip].plot(x[ip], energy[ib][ix_edges[ip]:ix_edges[ip+1]], color=color, label=label, zorder=zorder,
                                 linewidth=1)
                label = None

    def add_scatter(self, x, energy, size, color=None, zorder=2.5):

        if color is None:
            color = self.palette["orange"]

        ix_edges = [0] * (self.Np+1)
        for ip in range(self.Np):
            ix_edges[ip+1] = ix_edges[ip] + len(x[ip])

        Nb = len(energy)
        for ip in range(self.Np):
            for ib in range(Nb):
                if not isinstance(color, str):
                    current_color = color[ib][ix_edges[ip]:ix_edges[ip+1]]
                else:
                    current_color = color
                self.ax[ip].scatter(x[ip], energy[ib][ix_edges[ip]:ix_edges[ip+1]],
                                    s=size[ib][ix_edges[ip]:ix_edges[ip+1]], c=current_color, zorder=zorder)

    def plot_bands(self, eigenval1, kpoints1, rlc):

        self.x = eigenval1.eig_x(kp=kpoints1, rlc=rlc)
        self.energy = eigenval1.eigtrans()

        spin_label = ["majority spin", "minority spin"]
        line_color = ["darkblue", "orange"]

        for ispin in range(eigenval1.Ns):
            self.add_plot(self.x, self.energy[ispin],
                          color=self.palette[line_color[ispin]], label=spin_label, zorder=3-ispin)

        f1 = open("eigenval.dat", "w")
        if eigenval1.Ns == 2:
            for ib in range(len(self.energy[0])):
                for ip in range(self.Np):
                    ik0 = 0
                    for ik in range(len(self.x[ip])):
                        f1.write(str(self.x[ip][ik])+" "+str(self.energy[0][ib][ik0]) +
                                 " "+str(self.energy[1][ib][ik0])+"\n")
                        ik0 += 1
                f1.write("\n")
        else:
            for ib in range(len(self.energy[0])):
                for ip in range(self.Np):
                    ik0 = 0
                    for ik in range(len(self.x[ip])):
                        f1.write(str(self.x[ip][ik])+" "+str(self.energy[0][ib][ik0])+" "+"\n")
                        ik0 += 1
                f1.write("\n")
        f1.close()
