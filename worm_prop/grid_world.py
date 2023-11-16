from typing import *
import numpy as np
import matplotlib.pyplot as plt
from .lattice import Lattice

TSpaceTimePoint = TypeVar("TSpaceTimePoint", bound="SpaceTimePoint")


class SpaceTimePoint:
    def __init__(
        self,
        time_index: int,
        space_index: int,
        heads: List[TSpaceTimePoint],
        tails: List[TSpaceTimePoint],
    ):
        self.time_index = time_index
        self.space_index = space_index
        self.occupied: bool = False
        self.heads = heads
        self.tails = tails
        self.exception_list: List[TSpaceTimePoint] = []

    def __str__(self):
        return f"SpaceTimePoint at ({self.time_index}, {self.space_index})"

    def __repr__(self):
        return f"SpaceTimePoint at ({self.time_index}, {self.space_index})"


class GridWorld:
    def __init__(self, len_time: int, lattice: Lattice):
        self.len_time = len_time
        self.lattice = lattice
        self.len_space = self.lattice.N
        self.grid = [
            [
                SpaceTimePoint(i_time, i_space, [], [])
                for i_space in range(self.len_space)
            ]
            for i_time in range(len_time)
        ]

    # def trajectory_initialization(
    #     self,
    #     n_line: Optional[int] = 1,
    # ):
    #     positions = [i for i in range(self.lattice.N)]
    #     ## Draw lines
    #     line_positions = []
    #     for i in range(n_line):
    #         time_index = 0
    #         space_index = np.random.choice(positions)

    #         positions.remove(space_index)
    #         line_positions.append(space_index)

    #         node = self.grid[time_index][space_index]
    #         node.heads.append(self.grid[-1][space_index])
    #         self.grid[-1][space_index].tails.append(node)
    #         node.occupied = True

    #         for _ in range(self.len_time - 1):
    #             time_index = (time_index + 1) % self.len_time
    #             next_node = self.grid[time_index][space_index]
    #             next_node.occupied = True
    #             next_node.heads.append(node)
    #             node.tails.append(next_node)
    #             node = next_node

    #     return True

    def trajectory_initialization(self):
        time_index = 0
        space_index = np.random.randint(self.len_space)
        first_space_index = space_index
        node = self.grid[time_index][space_index]
        last_node = self.grid[-1][space_index]
        node.occupied = True
        last_node.occupied = True
        node.heads.append(last_node)
        last_node.tails.append(node)

        for _ in range(self.len_time - 3):
            space_index = self.lattice.hop[np.random.randint(1, 7)][space_index]
            time_index = (time_index + 1) % self.len_time
            next_node = self.grid[time_index][space_index]
            next_node.occupied = True
            next_node.heads.append(node)
            node.tails.append(next_node)
            node = next_node

        first_available_hop = set(self.lattice.hop[:, first_space_index])
        last_available_hop = set(self.lattice.hop[:, space_index])
        intersect = list(first_available_hop.intersection(last_available_hop))

        if len(intersect) > 0:
            space_index = np.random.choice(intersect)
            time_index = (time_index + 1) % self.len_time
            next_node = self.grid[time_index][space_index]
            next_node.occupied = True
            next_node.heads.append(node)
            node.tails.append(next_node)
            node = next_node

            time_index = (time_index + 1) % self.len_time
            next_node = self.grid[time_index][first_space_index]
            next_node.occupied = True
            next_node.heads.append(node)
            node.tails.append(next_node)
            node = next_node

            return True
        return False

    def occupation_map(self):
        occupation = np.zeros_like(self.grid).astype(bool)
        for time_index, sub_grid in enumerate(self.grid):
            for space_index, point in enumerate(sub_grid):
                if point.occupied:
                    occupation[time_index, space_index] = True

        return occupation

    def world_map(self):
        fig = plt.figure(figsize=(18, 4))
        ax = fig.add_subplot(111)
        ax.set_xlim([-1, self.len_time])
        ax.set_ylim([-1, self.len_space])
        ax.set_aspect("equal")
        for i in range(self.len_time):
            ax.vlines(x=i, ymin=-0.3, ymax=self.len_space, color="k", lw=3)
        for i in range(self.len_space):
            ax.hlines(y=i, xmin=-0.3, xmax=self.len_time, color="k", lw=3)

        for time_index, sub_grid in enumerate(self.grid):
            for space_index, point in enumerate(sub_grid):
                for tail in point.tails:
                    ax.plot(
                        [point.time_index, point.time_index + 1],
                        [point.space_index, tail.space_index],
                        c="r",
                        lw=5,
                    )

        ax.set_xticks(np.arange(self.len_time))
        ax.set_yticks(np.arange(self.len_space))
        ax.tick_params(axis="both", which="major", labelsize=10)
        fig.tight_layout()
        return fig, ax
