from typing import *
import numpy as np


class Lattice:
    def __init__(self, Nx: int, Ny: int, Nz: int):
        self.N = Nx * Ny * Nz
        self.index = np.arange(Nx * Ny * Nz).reshape(Nx, Ny, Nz)
        self.hop = [self.index.flatten()]

        self.hop.append(np.zeros_like(self.index))
        self.hop[-1][:-1] = self.index[1:]
        self.hop[-1][-1] = self.index[0]
        self.hop[-1] = self.hop[-1].flatten()

        self.hop.append(np.zeros_like(self.index))
        self.hop[-1][1:] = self.index[:-1]
        self.hop[-1][0] = self.index[-1]
        self.hop[-1] = self.hop[-1].flatten()

        self.hop.append(np.zeros_like(self.index))
        self.hop[-1][:, :-1] = self.index[:, 1:]
        self.hop[-1][:, -1] = self.index[:, 0]
        self.hop[-1] = self.hop[-1].flatten()

        self.hop.append(np.zeros_like(self.index))
        self.hop[-1][:, 1:] = self.index[:, :-1]
        self.hop[-1][:, 0] = self.index[:, -1]
        self.hop[-1] = self.hop[-1].flatten()

        self.hop.append(np.zeros_like(self.index))
        self.hop[-1][:, :, :-1] = self.index[:, :, 1:]
        self.hop[-1][:, :, -1] = self.index[:, :, 0]
        self.hop[-1] = self.hop[-1].flatten()

        self.hop.append(np.zeros_like(self.index))
        self.hop[-1][:, :, 1:] = self.index[:, :, :-1]
        self.hop[-1][:, :, 0] = self.index[:, :, -1]
        self.hop[-1] = self.hop[-1].flatten()

        self.hop = np.array(self.hop).astype(int)
        self.index = self.index.flatten()
