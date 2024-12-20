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

#include "RectangularGrid3D.h"


Engine::Engine(const std::string &inputFileName) {
    std::cout << "Initializing engine" << std::endl;
    //initialize
    _inputs = Readers::readJSON(inputFileName);
    std::cout << "Simulation data to be stored in " << _inputs.absolute_path_output.string() << std::endl;

    std::cout << std::format("Lx = {}, Ly = {}, Lz = {}",
                             _inputs.domain_len_x, _inputs.domain_len_y, _inputs.domain_len_z) << std::endl;
    std::cout << std::format("Nx = {}, Ny = {}, Nz = {}", _inputs.numx, _inputs.numy, _inputs.numz) << std::endl;
    std::cout << std::format("dx = {}, dy = {}, dz = {}", _inputs.dx, _inputs.dy, _inputs.dz) << std::endl;
}

void Engine::run() const {
    const Parameters inputs = _inputs; //make a const copy, to avoid accidental modification of input variables.

    auto index_path = inputs.absolute_path_output;

    auto xgrid = AuxiliaryFunctions::linspace(inputs.xmin, inputs.xmax, inputs.numx);
    auto ygrid = AuxiliaryFunctions::linspace(inputs.ymin, inputs.ymax, inputs.numy);
    auto zgrid = AuxiliaryFunctions::linspace(inputs.zmin, inputs.zmax, inputs.numz);

    const RectangularGrid3D grid(xgrid, ygrid, zgrid);

    // define geometry
    const Geometry geometry(inputs.shapes, inputs.xs_map);
    // set perfectly matched layers
    const PML pmly(inputs.pml_thickness, inputs.pml_strength, grid.get_ymin(), grid.get_ymax());
    //adds a PML on both sides
    const PML pmlz(inputs.pml_thickness, inputs.pml_strength, grid.get_zmin(), grid.get_zmax());

    // Define index monitors.
    IndexMonitor monitor_x0(grid.get_ymin(), grid.get_ymax(), grid.get_zmin(), grid.get_zmax(), 'x',
                            grid.get_xmin(),
                            inputs.numy, inputs.numz);
    monitor_x0.populate(geometry);
    index_path.append("index_yz_start.h5");
    monitor_x0.save_data(index_path.string());

    IndexMonitor monitor_x1(grid.get_ymin(), grid.get_ymax(), grid.get_zmin(), grid.get_zmax(), 'x',
                            grid.get_xmax(),
                            inputs.numy, inputs.numz);
    monitor_x1.populate(geometry);
    index_path.replace_filename("index_yz_end.h5");
    monitor_x1.save_data(index_path);

    IndexMonitor monitor_y(grid.get_xmin(), grid.get_xmax(), grid.get_zmin(), grid.get_zmax(), 'y', 0.0,
                           inputs.numx, inputs.numz);
    monitor_y.populate(geometry);
    index_path.replace_filename("index_xz.h5");
    monitor_y.save_data(index_path);

    IndexMonitor monitor_z(grid.get_xmin(), grid.get_xmax(), grid.get_ymin(), grid.get_ymax(), 'z', 0.0,
                           inputs.numx, inputs.numy);
    monitor_z.populate(geometry);
    index_path.replace_filename("index_xy.h5");
    monitor_z.save_data(index_path);

    if (inputs.dry_run) {
        std::cout << "Finished dry run." << std::endl;
        return;
    }

    /*
     * Create a 'mini' Solver for each port, that is used to find the modes.
     */
    //todo: mode solving is not always stable.
    std::vector<ModeSolver> input_solvers;
    std::vector<ModeSolver> output_solvers;
    double increment_x = grid.get_dx();
    // input ports
    for (const auto &input_port: inputs.input_ports) {
        // Im-Dis method not stable in combination with PML. So just for the mode solver we will use metal walls.
        const PML pmly_input(inputs.pml_thickness, 0 * inputs.pml_strength, input_port.get_ymin(),
                             input_port.get_ymax());
        const PML pmlz_input(inputs.pml_thickness, 0 * inputs.pml_strength, input_port.get_zmin(),
                             input_port.get_zmax());


        ModeSolver mode_solver(geometry, pmly_input, pmlz_input, input_port, inputs.scheme_parameter, inputs.k0,
                               inputs.reference_index, inputs.absolute_path_output, ERROR);
        mode_solver.run(increment_x);

        input_solvers.push_back(mode_solver);
    }
    // output ports
    for (const auto &output_port: inputs.output_ports) {
        const PML pmly_output(inputs.pml_thickness, 0 * inputs.pml_strength, output_port.get_ymin(),
                              output_port.get_ymax());
        const PML pmlz_output(inputs.pml_thickness, 0 * inputs.pml_strength, output_port.get_zmin(),
                              output_port.get_zmax());


        ModeSolver mode_solver(geometry, pmly_output, pmlz_output, output_port, inputs.scheme_parameter, inputs.k0,
                               inputs.reference_index, inputs.absolute_path_output, ERROR);
        mode_solver.run(increment_x);

        output_solvers.push_back(mode_solver);
    }

    /*
     * Now create the main simulation Solver, and run it.
     */
    BpmSolver bpm_solver(geometry, pmly, pmlz, grid, inputs.scheme_parameter, inputs.k0,
                         inputs.reference_index, inputs.absolute_path_output);

    // get initial field from the input mode.
    //todo: for now this only works for one input port at a time. If we have multiple input ports only the first is used.
    // update such that you can use multiple at the same time.
    auto intial_field = input_solvers.front().interpolate_field(grid.get_ygrid(), grid.get_zgrid());
    bpm_solver.run(intial_field);
    //todo: call function solver.save_data() or something to save/extract data after the solve is completed.

    auto output_port = inputs.output_ports.front();
    auto final_field = bpm_solver.interpolate_field(output_port.get_ygrid(), output_port.get_zgrid());

    double coef_b0 = output_solvers.front().get_mode_overlap(final_field);
    std::cout << std::format("Coupling coefficient is {}", coef_b0) << std::endl;
}
