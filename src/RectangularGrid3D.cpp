//
// Created by mike on 11/17/24.
//

#include "RectangularGrid3D.h"

RectangularGrid3D::RectangularGrid3D(std::vector<double> const &xgrid, std::vector<double> const &ygrid,
                                     std::vector<double> const &zgrid) : RectangularGrid(ygrid, zgrid), m_xgrid(xgrid) {
}

double RectangularGrid3D::get_xmin() const {
    return m_xgrid.front();
}

double RectangularGrid3D::get_xmax() const {
    return m_xgrid.back();
}

double RectangularGrid3D::get_dx() const {
    return m_xgrid[1] - m_xgrid[0];
}

size_t RectangularGrid3D::get_numx() const {
    return m_xgrid.size();
}

const std::vector<double> &RectangularGrid3D::get_xgrid() const {
    return m_xgrid;
}
