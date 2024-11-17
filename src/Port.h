//
// Created by mike on 11/11/24.
//

#ifndef PORT_H
#define PORT_H
#include <string>
#include "RectangularGrid.h"


class Port : public RectangularGrid {
public:
    Port(const std::string &name, const std::string &placement,
         double x0, std::vector<double> const &ygrid, std::vector<double> const &zgrid);

    [[nodiscard]] std::string get_name() const;

    [[nodiscard]] double get_x0() const;

private:
    const std::string name;
    const double x0;
    //todo: placement not used currently, remove?
    const std::string placement; //left (x=xmin) or right (x=xmax)

    //todo: add own grid?
};


#endif //PORT_H
