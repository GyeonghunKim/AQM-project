from enum import Enum
from typing import *

import numpy as np
import matplotlib.pyplot as plt

from .grid_world import SpaceTimePoint, GridWorld


class Direction(Enum):
    Incoming = 1
    Outgoing = -1


class TemporalDirection(Enum):
    Forward = 1
    Backward = -1


class Worm:
    def __init__(
        self,
        grid_world: GridWorld,
        direction: Direction,
        temporal_direction: TemporalDirection,
        position: SpaceTimePoint,
        mu: float,
        epsilon: float,
        beta: float,
    ):
        self.grid_world = grid_world
        self.direction = direction
        self.temporal_direction = temporal_direction
        self.position = position
        self.mu = mu
        self.epsilon = epsilon
        self.beta = beta

    def __str__(self):
        return f"Worm at ({self.position.time_index}, {self.position.space_index}) with {self.temporal_direction.name} and {self.direction.name}"

    def __repr__(self):
        return f"Worm at ({self.position.time_index}, {self.position.space_index}) with {self.temporal_direction.name} and {self.direction.name}"

    def update(self, n_iter: Optional[int] = None, show_step: Optional[int] = None):
        if n_iter is None:
            i_iter = 0
            while True:
                # print(self)
                if self.temporal_direction.value == TemporalDirection.Forward.value:
                    keepgo, do_site_update = self.forward_bond_update()
                    if do_site_update:
                        self.site_update()
                    if show_step is not None:
                        if i_iter % show_step == 0:
                            fig, ax = self.show_worm()
                            # print(self.position.heads, self.position.tails)
                    if not keepgo:
                        break
                else:
                    keepgo, do_site_update = self.backward_bond_update()
                    if do_site_update:
                        self.site_update()
                    if show_step is not None:
                        if i_iter % show_step == 0:
                            fig, ax = self.show_worm()
                            # print(self.position.heads, self.position.tails)
                    if not keepgo:
                        break
                i_iter += 1
        else:
            for _ in range(n_iter):
                if self.temporal_direction.value == TemporalDirection.Forward.value:
                    keepgo, do_site_update = self.forward_bond_update()
                    if do_site_update:
                        self.site_update()
                    if show_step is not None:
                        if i_iter % show_step == 0:
                            fig, ax = self.show_worm()
                            # print(self.position.heads, self.position.tails)
                    if not keepgo:
                        break
                else:
                    keepgo, do_site_update = self.backward_bond_update()
                    if do_site_update:
                        self.site_update()
                    if show_step is not None:
                        if i_iter % show_step == 0:
                            fig, ax = self.show_worm()
                            # print(self.position.heads, self.position.tails)
                    if not keepgo:
                        break

    def forward_bond_update(self) -> Tuple[bool, bool]:
        """Implementation of Page 64 and 67"""
        self.position.occupied = True
        # Determine hopping direction
        hopping_direction = int(np.random.rand() // self.epsilon + 1)
        if hopping_direction > 6:
            hopping_direction = 0

        new_position: SpaceTimePoint = self.grid_world.grid[
            (self.position.time_index + 1) % self.grid_world.len_time
        ][self.grid_world.lattice.hop[hopping_direction, self.position.space_index]]

        # Link current and next spacetimepoint
        self.position.tails.append(new_position)
        if len(new_position.heads) > 0:
            new_position.exception_list.append(self.position)
            new_position.heads.append(self.position)
            new_position.occupied = True

            # Update current state
            self.position = new_position
            self.direction = Direction.Outgoing
            self.temporal_direction = TemporalDirection.Backward
            return True, False

        new_position.heads.append(self.position)
        new_position.occupied = True

        # Update current state
        self.position = new_position
        self.direction = Direction.Incoming
        self.temporal_direction = TemporalDirection.Forward
        return True, True

    def backward_bond_update(self) -> Tuple[bool, bool]:
        """Implementation of Page 68"""
        self.position.occupied = False
        if len(self.position.heads) == 0:
            return False, False

        if len(self.position.heads) > 1:
            for pos in self.position.heads:
                if pos not in self.position.exception_list:
                    new_position = pos
                else:
                    self.position.exception_list.remove(pos)
        new_position = self.position.heads[0]
        self.position.heads.remove(new_position)
        new_position.tails.remove(self.position)
        self.direction = Direction.Incoming
        self.position = new_position
        return True, True

    def site_update(self):
        """Implementation of Page 65 66, 69, and 70"""
        self.direction = Direction(self.direction.value * -1)
        if np.random.rand() > min(
            1, np.exp(self.temporal_direction.value * self.mu * self.epsilon)
        ):
            self.temporal_direction = TemporalDirection(
                self.temporal_direction.value * -1
            )

    def show_worm(self):
        fig, ax = self.grid_world.world_map()
        ax.arrow(
            self.position.time_index,
            self.position.space_index,
            0.5 if self.temporal_direction.value == 1 else -0.5,
            0,
            width=0.2,
            head_width=0.5,
            head_length=0.2,
            color="b",
            length_includes_head=True,
            head_starts_at_zero=True if self.direction.value == 1 else False,
        )
        ax.set_xlim([self.position.time_index - 10, self.position.time_index + 10])
        fig.tight_layout()
        plt.show()
        return fig, ax
