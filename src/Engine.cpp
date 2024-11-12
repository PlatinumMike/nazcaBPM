//
// Created by mike on 9/29/24.
//

#include "Engine.h"
#include "IO/Readers.h"
#include "PML.h"
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
    // lot of code duplication here. Maybe grid info of a port should be in the Port class?
    // Also update the input file such that it contains a list of ports that we read in.
    Port inport0("a0", "left", 0.0, 0.0, 0.0, 6.0, 2.0);

    int numy_a0 = static_cast<int>(inport0.get_yspan() * inputs.resolution_y);
    int numz_a0 = static_cast<int>(inport0.get_zspan() * inputs.resolution_z);

    auto ygrid_a0 = AuxiliaryFunctions::linspace(inport0.get_ymin(), inport0.get_ymax(), numy_a0);
    auto zgrid_a0 = AuxiliaryFunctions::linspace(inport0.get_zmin(), inport0.get_zmax(), numz_a0);

    const RectangularGrid grid_a0(xgrid, ygrid_a0, zgrid_a0);

    //todo: some issue with the PMl inside of the mode solver, use metal walls for now, but check later on.
    const PML pmly_a0(inputs.pml_thickness, 0 * inputs.pml_strength, inport0.get_ymin(), inport0.get_ymax());
    const PML pmlz_a0(inputs.pml_thickness, 0 * inputs.pml_strength, inport0.get_zmin(), inport0.get_zmax());


    ModeSolver mode_solver(geometry, pmly_a0, pmlz_a0, grid_a0, inputs.scheme_parameter, inputs.k0,
                           inputs.reference_index, inport0);
    mode_solver.run();
    auto intial_field = mode_solver.interpolate_field(grid.get_ygrid(), grid.get_zgrid());


    //todo: mode solving is not always stable. It does not always converge. Also, a wider WG should have a larger neff, but that is not what the BPM finds...
    // output ports
    Port outport0("b0", "right", xgrid.back(), 0.0, 0.0, 6.0, 2.0);

    int numy_b0 = static_cast<int>(outport0.get_yspan() * inputs.resolution_y);
    int numz_b0 = static_cast<int>(outport0.get_zspan() * inputs.resolution_z);

    auto ygrid_b0 = AuxiliaryFunctions::linspace(outport0.get_ymin(), outport0.get_ymax(), numy_b0);
    auto zgrid_b0 = AuxiliaryFunctions::linspace(outport0.get_zmin(), outport0.get_zmax(), numz_b0);

    const RectangularGrid grid_b0(xgrid, ygrid_b0, zgrid_b0);

    const PML pmly_b0(inputs.pml_thickness, 0 * inputs.pml_strength, outport0.get_ymin(), outport0.get_ymax());
    const PML pmlz_b0(inputs.pml_thickness, 0 * inputs.pml_strength, outport0.get_zmin(), outport0.get_zmax());


    ModeSolver mode_solver_out(geometry, pmly_b0, pmlz_b0, grid_b0, inputs.scheme_parameter, inputs.k0,
                               inputs.reference_index, outport0);
    mode_solver_out.run();

    /*
     * Now create the main simulation Solver, and run it.
     */
    BpmSolver bpm_solver(geometry, pmly, pmlz, grid, inputs.scheme_parameter, inputs.k0,
                         inputs.reference_index);
    bpm_solver.run(intial_field);
    //todo: call function solver.save_data() or something to save/extract data after the solve is completed.

    auto final_field = bpm_solver.interpolate_field(ygrid_b0, zgrid_b0);

    double coef_b0 = mode_solver_out.get_mode_overlap(final_field);
    std::cout << std::format("Coupling coefficient is {}", coef_b0) << std::endl;
}
