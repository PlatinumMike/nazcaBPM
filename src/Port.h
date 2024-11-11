//
// Created by mike on 11/11/24.
//

#ifndef PORT_H
#define PORT_H
#include <string>


class Port {
public:
    Port(const std::string &name, const std::string &placement, double x0, double y0, double z0, double yspan,
         double zspan);

    [[nodiscard]] std::string get_name() const;

    [[nodiscard]] double get_x0() const;

    [[nodiscard]] double get_y0() const;

    [[nodiscard]] double get_yspan() const;

    [[nodiscard]] double get_z0() const;

    [[nodiscard]] double get_zspan() const;

private:
    const std::string name;
    const double x0;
    const double y0;
    const double z0;
    const double yspan;
    const double zspan;
    //todo: placement not used currently, remove?
    const std::string placement; //left (x=xmin) or right (x=xmax)

    //todo: add own grid?
};


#endif //PORT_H
