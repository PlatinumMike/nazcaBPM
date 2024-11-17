//
// Created by mike on 11/17/24.
//

#ifndef RECTANGULARGRID3D_H
#define RECTANGULARGRID3D_H

#include "RectangularGrid.h"


/**
 * Class that holds 3D grid information. It simply extends the normal 2D grid by adding a xgrid.
 * This split is done because not every class needs all three axes. It is best to pass only the information that is used
 * so we do not accidentally use the xgrid when it is not allowed.
 */
class RectangularGrid3D : public RectangularGrid {
public:
    RectangularGrid3D(std::vector<double> const &xgrid, std::vector<double> const &ygrid,
                      std::vector<double> const &zgrid);

    double get_xmin() const;


    double get_xmax() const;


    double get_dx() const;


    size_t get_numx() const;


    const std::vector<double> &get_xgrid() const;

private:
    const std::vector<double> m_xgrid;
};


#endif //RECTANGULARGRID3D_H
