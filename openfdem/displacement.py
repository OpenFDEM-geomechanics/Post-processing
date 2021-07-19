import glob
import os
import os.path as path
import re

import pyvista as pv
import time
import concurrent.futures
from multiprocessing import Process
from threading import Thread
import openfdem as fd


def platen_displacement(self, material_id=None, boundary_condition_id=None, location=None):
    # TODO: load based on threshold points (boundary condition)
    #   can not do as there is a -1 boundary condition! And
    #   this should be limited only to the platen boundary condition (e.g. confinement)
    '''
    Checks if Compression Simulations.
    Calculate the strain based on the displacement in the platens.

    :return: strain (%) = displacement / height of sample
    :rtype: list

    Example:
        >>> data = pv.read("../example_outputs/Irazu_UCS")
        >>> axial_disp = data.platen_displacement()
        [0.0, 0.009259259235881877, 0.018518518471763754, 0.02777777770764563, 0.03703703694352751, 0.046296296179409384, 0.05555555541529126, 0.06481481465117314, 0.07407407388705502, 0.08333333312293689, 0.09259259235881877, 0.10185185159470064, 0.11111111083058252, 0.1203703700664644, 0.12962962930234628, 0.13888888853822817, 0.14814814777411003, 0.15740740700999192, 0.16666666624587378, 0.17592592548175565, 0.18518518471763754, 0.1944444439535194, 0.2037037031894013, 0.21296296242528318, 0.22222222166116504, 0.23148148089704693, 0.2407407401329288, 0.24999999936881068, 0.25925925860469257, 0.2685185178405744, 0.27777777707645634]

    '''

    history_strain = []  # List to hold the stress_history
    avg_platen_disp = [0.0, 0.0, 0.0]  # Dummy cell
    axis_of_loading = 1  # Axis of loading in Y direction.

    ## Get rock dimension.
    self.rock_sample_dimensions()
    ## Check UCS Simulation
    if self.simulation_type() != "UCS Simulation":
        print("Simulation appears to be not for compressive strength")

    # Load each timestep
    # backup copy compressed below

    def history_strain_func(openfdem_model_ts):
        platen = (openfdem_model_ts.threshold([self.platen_cells_elem_id, self.platen_cells_elem_id],
                                              self.var_data["mineral_type"]))
        top, bottom = (platen.get_data_range(self.var_data["boundary"]))

        avg_top_platen_disp = self.platen_info(openfdem_model_ts, top, self.var_data["platen_displacement"])
        avg_bottom_platen_disp = self.platen_info(openfdem_model_ts, bottom,
                                                  self.var_data["platen_displacement"])

        for i in range(0, self.number_of_points_per_cell):
            avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

        strain_from_platen = avg_platen_disp[axis_of_loading] / self.sample_height * 100.0
        history_strain.append(strain_from_platen)
        return strain_from_platen

    def main():
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = executor.map(history_strain_func, self)  # is self the list we are iterating over
        return results

    if __name__ == '__main__':
        main()

    return history_strain
