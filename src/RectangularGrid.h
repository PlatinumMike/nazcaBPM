//
// Created by mike on 10/22/24.
//

#ifndef RECTANGULARGRID_H
#define RECTANGULARGRID_H

#include <vector>

#include <boost/multi_array.hpp>
using boost::multi_array;

/**
 * Class that holds 2D grid information
 */
class RectangularGrid {
public:
    RectangularGrid(std::vector<double> const &ygrid, std::vector<double> const &zgrid);

    double get_ymin() const;

    double get_ymax() const;

    double get_zmin() const;

    double get_zmax() const;

    double get_dy() const;

    double get_dz() const;

    size_t get_numy() const;

    size_t get_numz() const;

    [[nodiscard]] double get_y0() const;

    [[nodiscard]] double get_yspan() const;

    [[nodiscard]] double get_z0() const;

    [[nodiscard]] double get_zspan() const;

    const std::vector<double> &get_ygrid() const;

    const std::vector<double> &get_zgrid() const;

    static multi_array<double, 1> vector_to_multi_array(const std::vector<double> &vec);

private:
    const std::vector<double> m_ygrid;
    const std::vector<double> m_zgrid;
};


#endif //RECTANGULARGRID_H
