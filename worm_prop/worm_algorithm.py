from typing import *

import numpy as np
from tqdm.notebook import tqdm

from .lattice import Lattice
from .grid_world import GridWorld, SpaceTimePoint
from .worm import Worm, Direction, TemporalDirection


class WormAlgorithm:
    def __init__(self, lattice: Lattice, mu: float, epsilon: float, beta: float):
        self.mu = mu
        self.epsilon = epsilon
        self.beta = beta
        self.time_length = int(beta / epsilon)
        self.lattice = lattice
        self.grid_world = self.set_init_conf()
        self.Z = 0
        self.sum_nW = 0
        self.sum_eW = 0
        self.W_list = []
        self.Z_list = []
        self.e_list = []
        self.n_list = []
        self.n_shop_list = []
        self.n_thop_list = []

    def set_init_conf(self):
        grid_world = GridWorld(self.time_length, self.lattice)
        while not grid_world.trajectory_initialization():
            grid_world = GridWorld(self.time_length, self.lattice)
            continue
        return grid_world

    def single_run(self) -> Tuple[float, float]:
        self.worm_initialization()
        self.worm.update()
        n_exp, e_exp = self.calculate_expectation()
        return n_exp, e_exp

    def run(self, n_iteration):
        n_list = []
        e_list = []
        for _ in tqdm(range(n_iteration), total=n_iteration):
            n_exp, e_exp = self.single_run()
            n_list.append(n_exp)
            e_list.append(e_exp)
        return np.array(n_list), np.array(e_list)

    def calculate_expectation(self) -> Tuple[float, float]:
        n_spatial_hop: int = 0
        n_temporal_hop: int = 0
        n: int = 0
        for space_index_start in range(self.lattice.N):
            cursor = self.grid_world.grid[0][space_index_start]
            if not cursor.occupied:
                continue
            n += 1
            for i in range(self.time_length - 1):
                next = cursor.tails[0]
                if next.space_index == cursor.space_index:
                    n_temporal_hop += 1
                else:
                    n_spatial_hop += 1
                cursor = next
        W = (
            (self.epsilon) ** n_spatial_hop
            * (1 - 6 * self.epsilon) ** n_temporal_hop
            * (np.exp(self.mu * self.epsilon)) ** (n_spatial_hop + n_temporal_hop)
        )

        k = self.beta / self.epsilon
        e = (-1.0 / (k * self.epsilon)) * n_spatial_hop + (
            6 / (k * (1 - 6 * self.epsilon))
        ) * n_temporal_hop
        self.Z += W
        self.sum_eW += e * W
        self.sum_nW += n * W
        self.Z_list.append(self.Z)
        self.e_list.append(e)
        self.n_list.append(n)
        self.W_list.append(W)
        self.n_shop_list.append(n_spatial_hop)
        self.n_thop_list.append(n_temporal_hop)
        return self.sum_nW / self.Z, self.sum_eW / self.Z

    def worm_initialization(self):
        ## Initialize the worm
        # To begin we pick a lattice site at random
        rand_t, rand_r = (
            np.random.randint(self.time_length),
            np.random.randint(self.lattice.N),
        )
        # If the picked site is empty
        cursor = self.grid_world.grid[rand_t][rand_r]
        if not cursor.occupied:
            if np.random.rand() < min(1, np.exp(self.mu * self.epsilon)):
                # Begin Worm Sector
                self.worm = Worm(
                    self.grid_world,
                    Direction.Outgoing,
                    TemporalDirection.Forward,
                    cursor,
                    self.mu,
                    self.epsilon,
                    self.beta,
                )
            else:
                # Back to partition function sectorz
                return
        # If the picked site has a particle worldline
        else:
            # Begin Worm Sector
            self.worm = Worm(
                self.grid_world,
                Direction.Outgoing,
                TemporalDirection.Backward,
                cursor,
                self.mu,
                self.epsilon,
                self.beta,
            )
