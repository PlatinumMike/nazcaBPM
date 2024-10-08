//
// Created by mike on 10/7/24.
//

#ifndef GRIDINTERPOLATOR_H
#define GRIDINTERPOLATOR_H
#include <utility>
#include <vector>
#include <boost/multi_array.hpp>


/**simple bilinear interpolation on a 2D grid.
 * the grid is assumed to be rectangular.
 * Note, this is not optimized at all. Probably much faster searching is possible but the current implementation is fast enough.
*/

template<typename T>

class GridInterpolator {
public:
    GridInterpolator(std::vector<double> gridCoordinate1, std::vector<double> gridCoordinate2,
                     boost::multi_array<T, 2> data, T defaultValue);

    T get_value(double coordinate1, double coordinate2) const;

private:
    //can be the x direction, y, z, angle, or whatever.
    // the precise meaning of "coordinate1" is irrelevant, we are just given a grid, and we will interpolate it.
    // the only thing that matters is that coordinate1 matches the first dimension of "data", and coordinate2 matches the second.
    std::vector<double> m_gridCoordinate1;
    std::vector<double> m_gridCoordinate2;
    boost::multi_array<T, 2> m_data;
    T m_defaultValue;

    int get_grid_index_lower(double x, const std::vector<double> &grid) const;
};


template<typename T>
GridInterpolator<T>::GridInterpolator(std::vector<double> gridCoordinate1, std::vector<double> gridCoordinate2,
                                      boost::multi_array<T, 2> data, T defaultValue) {
    assert(("array size does not match grid size, dimension 0", data.shape()[0]==gridCoordinate1.size()));
    assert(("array size does not match grid size, dimension 1", data.shape()[1]==gridCoordinate2.size()));
    m_gridCoordinate1 = gridCoordinate1;
    m_gridCoordinate2 = gridCoordinate2;
    m_defaultValue = defaultValue;
    // make deep copy
    const int num1 = data.shape()[0];
    const int num2 = data.shape()[1];
    m_data.resize(boost::extents[num1][num2]);
    m_data = data;
}

template<typename T>
T GridInterpolator<T>::get_value(const double coordinate1, const double coordinate2) const {
    //Note, this will be useful for evaluating the field at arbitrary points, as well as interpolating the mode source on the grid.

    // check if outside of grid
    if (coordinate1 < m_gridCoordinate1.front() || coordinate1 > m_gridCoordinate1.back()) {
        return m_defaultValue;
    }
    if (coordinate2 < m_gridCoordinate2.front() || coordinate2 > m_gridCoordinate2.back()) {
        return m_defaultValue;
    }

    //todo: handle edge case when x=xgrid[10] for instance, so there is no lower index because it is exactly on a grid point.

    const int idx_lower = get_grid_index_lower(coordinate1, m_gridCoordinate1);
    const int idx_upper = idx_lower + 1;
    const int idx_lower2 = get_grid_index_lower(coordinate2, m_gridCoordinate2);
    const int idx_upper2 = idx_lower2 + 1;

    const double dx1 = m_gridCoordinate1[idx_upper] - m_gridCoordinate1[idx_lower];
    const double dx2 = m_gridCoordinate2[idx_upper2] - m_gridCoordinate2[idx_lower2];
    double coef_sw = (m_gridCoordinate1[idx_upper] - coordinate1) * (m_gridCoordinate2[idx_upper2] - coordinate2);
    double coef_se = (coordinate1 - m_gridCoordinate1[idx_lower]) * (m_gridCoordinate2[idx_upper2] - coordinate2);
    double coef_nw = (m_gridCoordinate1[idx_upper] - coordinate1) * (coordinate2 - m_gridCoordinate2[idx_lower2]);
    double coef_ne = (coordinate1 - m_gridCoordinate1[idx_lower]) * (coordinate2 - m_gridCoordinate2[idx_lower2]);
    return (m_data[idx_lower][idx_lower2] * coef_sw + m_data[idx_upper][idx_lower2] * coef_se + m_data[idx_lower][
                idx_upper2] * coef_nw + m_data[idx_upper][idx_upper2] * coef_ne) / (dx1 * dx2);
}

template<typename T>
int GridInterpolator<T>::get_grid_index_lower(const double x, const std::vector<double> &grid) const {
    for (int i = 0; i < grid.size(); i++) {
        if (grid[i] > x) {
            return i - 1;
        }
    }
    // last element is N-1, so return N-2 as the lower index.
    return static_cast<int>(grid.size() - 2);
}
#endif //GRIDINTERPOLATOR_H
