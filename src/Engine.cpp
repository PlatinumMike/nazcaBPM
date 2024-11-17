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
#include <format>

#include "RectangularGrid.h"


Engine::Engine(const std::string &inputFileName) {
    std::cout << "Initializing engine" << std::endl;
    //initialize
    _inputs = Readers::readJSON(inputFileName);

    std::cout << std::format("Lx = {}, Ly = {}, Lz = {}\n",
                             _inputs.domain_len_x, _inputs.domain_len_y, _inputs.domain_len_z);
    std::cout << std::format("Nx = {}, Ny = {}, Nz = {}\n", _inputs.numx, _inputs.numy, _inputs.numz);
    std::cout << std::format("dx = {}, dy = {}, dz = {}\n", _inputs.dx, _inputs.dy, _inputs.dz);
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
    //todo: lot of code duplication here. Maybe grid info of a port should be in the Port class?
    //todo: mode solving is not always stable. It does not always converge. Also, a wider WG should have a larger neff, but that is not what the BPM finds...
    // even if I use a straight waveguide it somehow does not find the same mode at the start and end... investigate further.
    // output ports
    std::vector<ModeSolver> input_solvers;
    std::vector<ModeSolver> output_solvers;
    for (const auto &input_port: inputs.input_ports) {
        int numy_input = static_cast<int>(input_port.get_yspan() * inputs.resolution_y);
        int numz_input = static_cast<int>(input_port.get_zspan() * inputs.resolution_z);

        auto ygrid_input = AuxiliaryFunctions::linspace(input_port.get_ymin(), input_port.get_ymax(), numy_input);
        auto zgrid_input = AuxiliaryFunctions::linspace(input_port.get_zmin(), input_port.get_zmax(), numz_input);

        const RectangularGrid grid_input(xgrid, ygrid_input, zgrid_input);

        //todo: some issue with the PMl inside of the mode solver, use metal walls for now, but check later on.
        const PML pmly_input(inputs.pml_thickness, 0 * inputs.pml_strength, input_port.get_ymin(),
                             input_port.get_ymax());
        const PML pmlz_input(inputs.pml_thickness, 0 * inputs.pml_strength, input_port.get_zmin(),
                             input_port.get_zmax());


        ModeSolver mode_solver(geometry, pmly_input, pmlz_input, grid_input, inputs.scheme_parameter, inputs.k0,
                               inputs.reference_index, input_port);
        mode_solver.run();

        input_solvers.push_back(mode_solver);
    }

    for (const auto &output_port: inputs.output_ports) {
        int numy_output = static_cast<int>(output_port.get_yspan() * inputs.resolution_y);
        int numz_output = static_cast<int>(output_port.get_zspan() * inputs.resolution_z);

        auto ygrid_output = AuxiliaryFunctions::linspace(output_port.get_ymin(), output_port.get_ymax(), numy_output);
        auto zgrid_output = AuxiliaryFunctions::linspace(output_port.get_zmin(), output_port.get_zmax(), numz_output);

        const RectangularGrid grid_output(xgrid, ygrid_output, zgrid_output);

        const PML pmly_output(inputs.pml_thickness, 0 * inputs.pml_strength, output_port.get_ymin(),
                              output_port.get_ymax());
        const PML pmlz_output(inputs.pml_thickness, 0 * inputs.pml_strength, output_port.get_zmin(),
                              output_port.get_zmax());


        ModeSolver mode_solver(geometry, pmly_output, pmlz_output, grid_output, inputs.scheme_parameter, inputs.k0,
                                   inputs.reference_index, output_port);
        mode_solver.run();

        output_solvers.push_back(mode_solver);
    }

    /*
     * Now create the main simulation Solver, and run it.
     */
    BpmSolver bpm_solver(geometry, pmly, pmlz, grid, inputs.scheme_parameter, inputs.k0,
                         inputs.reference_index);

    // get initial field from the input mode.
    //todo: for now this only works for one input port at a time. If we have multiple input ports only the first is used.
    // update such that you can use multiple at the same time.
    auto intial_field = input_solvers.front().interpolate_field(grid.get_ygrid(), grid.get_zgrid());
    bpm_solver.run(intial_field);
    //todo: call function solver.save_data() or something to save/extract data after the solve is completed.

    //todo: update the bit below.
    // auto ygrid_output = AuxiliaryFunctions::linspace(output_port.get_ymin(), output_port.get_ymax(), numy_output);
    // auto zgrid_output = AuxiliaryFunctions::linspace(output_port.get_zmin(), output_port.get_zmax(), numz_output);
    // auto final_field = bpm_solver.interpolate_field(ygrid_output, zgrid_output);
    //
    // double coef_b0 = output_solvers.front().get_mode_overlap(final_field);
    // std::cout << std::format("Coupling coefficient is {}", coef_b0) << std::endl;
}
