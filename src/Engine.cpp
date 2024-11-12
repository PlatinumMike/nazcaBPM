//
// Created by mike on 9/29/24.
//

#include "Engine.h"
#include "IO/Readers.h"
#include "PML.h"
#include "ModeHandler.h"
#include "AuxiliaryFunctions.h"
#include "IndexMonitor.h"
#include "BpmSolver.h"
#include "ModeSolver.h"
#include "Port.h"

#include <iostream>
#include <ostream>
#include <chrono>

#include "RectangularGrid.h"


Engine::Engine(const std::string &inputFileName) {
    std::cout << "Initializing engine" << std::endl;
    //initialize
    _inputs = Readers::readJSON(inputFileName);

    std::cout << "Lx = " << _inputs.domain_len_x << ", Ly = " << _inputs.domain_len_y << ", Lz = " << _inputs.
            domain_len_z << std::endl;
    std::cout << "Nx = " << _inputs.numx << ", Ny = " << _inputs.numy << ", Nz = " << _inputs.numz << std::endl;
    std::cout << "dx = " << _inputs.dx << ", dy = " << _inputs.dy << ", dz = " << _inputs.dz << std::endl;
}

void Engine::run() const {
    const Parameters inputs = _inputs; //make a const copy, to avoid accidental modification of input variables.

    //TODO: loop over ports in input file, then for each do a mini simulation to get the modes.
    // then use one of those ports as the mode source for the real simulation.

    auto xgrid = AuxiliaryFunctions::linspace(inputs.xmin, inputs.xmax, inputs.numx);
    auto ygrid = AuxiliaryFunctions::linspace(inputs.ymin, inputs.ymax, inputs.numy);
    auto zgrid = AuxiliaryFunctions::linspace(inputs.zmin, inputs.zmax, inputs.numz);

    const RectangularGrid grid(xgrid, ygrid, zgrid);

    // define geometry
    const Geometry geometry(inputs.shapes, inputs.background_index);
    // set perfectly matched layers
    const PML pmly(inputs.pml_thickness, inputs.pml_strength, grid.get_ymin(), grid.get_ymax());
    //adds a PML on both sides
    const PML pmlz(inputs.pml_thickness, inputs.pml_strength, grid.get_zmin(), grid.get_zmax());
    // setup mode source
    const ModeHandler gaussianSource("Gaussian");

    // Define index monitors.
    IndexMonitor monitor_x0(grid.get_ymin(), grid.get_ymax(), grid.get_zmin(), grid.get_zmax(), 'x',
                            0.0,
                            inputs.numy, inputs.numz);
    monitor_x0.populate(geometry);
    monitor_x0.save_data("index_yz_start.h5");

    IndexMonitor monitor_x1(grid.get_ymin(), grid.get_ymax(), grid.get_zmin(), grid.get_zmax(), 'x',
                            xgrid.back(),
                            inputs.numy, inputs.numz);
    monitor_x1.populate(geometry);
    monitor_x1.save_data("index_yz_end.h5");

    IndexMonitor monitor_y(grid.get_xmin(), grid.get_xmax(), grid.get_zmin(), grid.get_zmax(), 'y', 0.0,
                           inputs.numx, inputs.numz);
    monitor_y.populate(geometry);
    monitor_y.save_data("index_xz.h5");

    IndexMonitor monitor_z(grid.get_xmin(), grid.get_xmax(), grid.get_ymin(), grid.get_ymax(), 'z', 0.0,
                           inputs.numx, inputs.numy);
    monitor_z.populate(geometry);
    monitor_z.save_data("index_xy.h5");

    if (inputs.dry_run) {
        std::cout << "Finished dry run." << std::endl;
        return;
    }

    /*
     * Create a 'mini' Solver for each port, that is used to find the modes.
     */
    //todo: loop over all in/out ports, get modes, then get input profile and pass that into the BpmSolver.
    Port inport0("a0", "left", 0.0, 0.0, 0.0, 4.0, 4.0);

    int numy_a0 = static_cast<int>(inport0.get_yspan() * inputs.resolution_y);
    int numz_a0 = static_cast<int>(inport0.get_zspan() * inputs.resolution_z);

    auto ygrid_a0 = AuxiliaryFunctions::linspace(inport0.get_ymin(), inport0.get_ymax(), numy_a0);
    auto zgrid_a0 = AuxiliaryFunctions::linspace(inport0.get_zmin(), inport0.get_zmax(), numz_a0);

    const RectangularGrid grid_a0(xgrid, ygrid_a0, zgrid_a0);

    //todo: some issue with the PMl inside of the mode solver, use metal walls for now, but check later on.
    const PML pmly_a0(inputs.pml_thickness, 0 * inputs.pml_strength, inport0.get_ymin(), inport0.get_ymax());
    const PML pmlz_a0(inputs.pml_thickness, 0 * inputs.pml_strength, inport0.get_zmin(), inport0.get_zmax());


    //todo: the mode solver should use a smaller grid that it gets from the port, not same big BPM grid.
    ModeSolver mode_solver(geometry, pmly_a0, pmlz_a0, gaussianSource, grid_a0, inputs.scheme_parameter, inputs.k0,
                           inputs.reference_index, inport0);
    mode_solver.run();
    auto intial_field = mode_solver.interpolate_field(grid.get_ygrid(), grid.get_zgrid());

    /*
     * Now create the main simulation Solver, and run it.
     */
    BpmSolver bpm_solver(geometry, pmly, pmlz, grid, inputs.scheme_parameter, inputs.k0,
                         inputs.reference_index);
    bpm_solver.run(intial_field);
    //todo: call function solver.save_data() or something to save/extract data after the solve is completed.

    //todo: use bpm_solver.interpolate_field() to get the field at the output pins, and do mode overlap.
}
