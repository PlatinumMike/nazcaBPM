//
// Created by mike on 10/22/24.
//

#include "RectangularGrid.h"
using boost::extents;


RectangularGrid::RectangularGrid(std::vector<double> const &xgrid, std::vector<double> const &ygrid,
                                 std::vector<double> const &zgrid) : m_xgrid(xgrid), m_ygrid(ygrid), m_zgrid(zgrid) {
}

double RectangularGrid::get_xmin() const {
    return m_xgrid.front();
}

double RectangularGrid::get_ymin() const {
    return m_ygrid.front();
}

double RectangularGrid::get_xmax() const {
    return m_xgrid.back();
}

double RectangularGrid::get_ymax() const {
    return m_ygrid.back();
}

double RectangularGrid::get_zmin() const {
    return m_zgrid.front();
}

double RectangularGrid::get_zmax() const {
    return m_zgrid.back();
}

double RectangularGrid::get_dx() const {
    return m_xgrid[1] - m_xgrid[0];
}

double RectangularGrid::get_dy() const {
    return m_ygrid[1] - m_ygrid[0];
}

double RectangularGrid::get_dz() const {
    return m_zgrid[1] - m_zgrid[0];
}

size_t RectangularGrid::get_numx() const {
    return m_xgrid.size();
}

size_t RectangularGrid::get_numy() const {
    return m_ygrid.size();
}

size_t RectangularGrid::get_numz() const {
    return m_zgrid.size();
}

const std::vector<double> &RectangularGrid::get_xgrid() const {
    return m_xgrid;
}

const std::vector<double> &RectangularGrid::get_ygrid() const {
    return m_ygrid;
}

const std::vector<double> &RectangularGrid::get_zgrid() const {
    return m_zgrid;
}

multi_array<double, 1> RectangularGrid::vector_to_multi_array(const std::vector<double> &vec) {
    const int size = static_cast<int>(vec.size());
    multi_array<double, 1> vec_m(extents[size]);
    for (int i = 0; i < size; i++) {
        vec_m[i] = vec[i];
    }
    return vec_m;
}