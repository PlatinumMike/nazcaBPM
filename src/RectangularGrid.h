//
// Created by mike on 10/22/24.
//

#ifndef RECTANGULARGRID_H
#define RECTANGULARGRID_H

#include <vector>

#include <boost/multi_array.hpp>
using boost::multi_array;


class RectangularGrid {
public:
    RectangularGrid(std::vector<double> const &xgrid, std::vector<double> const &ygrid,
                    std::vector<double> const &zgrid);

    double get_xmin() const;

    double get_ymin() const;

    double get_xmax() const;

    double get_ymax() const;

    double get_zmin() const;

    double get_zmax() const;

    double get_dx() const;

    double get_dy() const;

    double get_dz() const;

    size_t get_numx() const;

    size_t get_numy() const;

    size_t get_numz() const;

    const std::vector<double> &get_xgrid() const;

    const std::vector<double> &get_ygrid() const;

    const std::vector<double> &get_zgrid() const;

    static multi_array<double, 1> vector_to_multi_array(const std::vector<double> &vec);

private:
    const std::vector<double> m_xgrid;
    const std::vector<double> m_ygrid;
    const std::vector<double> m_zgrid;
};


#endif //RECTANGULARGRID_H
