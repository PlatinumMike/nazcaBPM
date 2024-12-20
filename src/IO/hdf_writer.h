//
// Copied this from https://stackoverflow.com/a/15221077
//

#ifndef HDF_WRITER_H
#define HDF_WRITER_H

#include <vector>
#include <cstdint>
#include <H5Cpp.h>
#include <boost/multi_array.hpp>

//!_______________________________________________________________________________________
//!
//!     map types to HDF5 types
//!
//!
//!     \author lg (04 March 2013)
//!_______________________________________________________________________________________

template<typename T> struct get_hdf5_data_type
{   static H5::PredType type()
    {
        //static_assert(false, "Unknown HDF5 data type");
        return H5::PredType::NATIVE_DOUBLE;
    }
};
template<> struct get_hdf5_data_type<char>                  {   H5::IntType type    {   H5::PredType::NATIVE_CHAR       };  };
//template<> struct get_hdf5_data_type<unsigned char>       {   H5::IntType type    {   H5::PredType::NATIVE_UCHAR      };  };
//template<> struct get_hdf5_data_type<short>               {   H5::IntType type    {   H5::PredType::NATIVE_SHORT      };  };
//template<> struct get_hdf5_data_type<unsigned short>      {   H5::IntType type    {   H5::PredType::NATIVE_USHORT     };  };
//template<> struct get_hdf5_data_type<int>                 {   H5::IntType type    {   H5::PredType::NATIVE_INT        };  };
//template<> struct get_hdf5_data_type<unsigned int>        {   H5::IntType type    {   H5::PredType::NATIVE_UINT       };  };
//template<> struct get_hdf5_data_type<long>                {   H5::IntType type    {   H5::PredType::NATIVE_LONG       };  };
//template<> struct get_hdf5_data_type<unsigned long>       {   H5::IntType type    {   H5::PredType::NATIVE_ULONG      };  };
template<> struct get_hdf5_data_type<long long>             {   H5::IntType type    {   H5::PredType::NATIVE_LLONG      };  };
template<> struct get_hdf5_data_type<unsigned long long>    {   H5::IntType type    {   H5::PredType::NATIVE_ULLONG     };  };
template<> struct get_hdf5_data_type<int8_t>                {   H5::IntType type    {   H5::PredType::NATIVE_INT8       };  };
template<> struct get_hdf5_data_type<uint8_t>               {   H5::IntType type    {   H5::PredType::NATIVE_UINT8      };  };
template<> struct get_hdf5_data_type<int16_t>               {   H5::IntType type    {   H5::PredType::NATIVE_INT16      };  };
template<> struct get_hdf5_data_type<uint16_t>              {   H5::IntType type    {   H5::PredType::NATIVE_UINT16     };  };
template<> struct get_hdf5_data_type<int32_t>               {   H5::IntType type    {   H5::PredType::NATIVE_INT32      };  };
template<> struct get_hdf5_data_type<uint32_t>              {   H5::IntType type    {   H5::PredType::NATIVE_UINT32     };  };
template<> struct get_hdf5_data_type<int64_t>               {   H5::IntType type    {   H5::PredType::NATIVE_INT64      };  };
template<> struct get_hdf5_data_type<uint64_t>              {   H5::IntType type    {   H5::PredType::NATIVE_UINT64     };  };
template<> struct get_hdf5_data_type<float>                 {   H5::FloatType type  {   H5::PredType::NATIVE_FLOAT      };  };
template<> struct get_hdf5_data_type<double>                {   H5::FloatType type  {   H5::PredType::NATIVE_DOUBLE     };  };
template<> struct get_hdf5_data_type<long double>           {   H5::FloatType type  {   H5::PredType::NATIVE_LDOUBLE    };  };

//!_______________________________________________________________________________________
//!
//!     write_hdf5 multi_array
//!
//!     \author leo Goodstadt (04 March 2013)
//!
//!_______________________________________________________________________________________
template<typename T, std::size_t DIMENSIONS, typename hdf5_data_type>
void do_write_hdf5(H5::H5File file, const std::string& data_set_name, const boost::multi_array<T, DIMENSIONS>& data, hdf5_data_type& datatype)
{
    // Little endian for x86
    //FloatType datatype(get_hdf5_data_type<T>::type());
    datatype.setOrder(H5T_ORDER_LE);

    std::vector<hsize_t> dimensions(data.shape(), data.shape() + DIMENSIONS);
    H5::DataSpace dataspace(DIMENSIONS, dimensions.data());

    H5::DataSet dataset = file.createDataSet(data_set_name, datatype, dataspace);

    dataset.write(data.data(), datatype);
}

template<typename T, std::size_t DIMENSIONS>
void write_hdf5(H5::H5File file, const std::string& data_set_name, const boost::multi_array<T, DIMENSIONS>& data )
{

    get_hdf5_data_type<T> hdf_data_type;
    do_write_hdf5(file, data_set_name, data, hdf_data_type.type);
}

inline void write_cmplx_hdf5(const std::string &filename,
                             const boost::multi_array<std::complex<double>, 2> &data, const multi_array<double,1>& grid1,
                             const multi_array<double,1>& grid2, const std::string &label1, const std::string &label2) {
    const auto num1 = static_cast<int>(data.shape()[0]);
    const auto num2 = static_cast<int>(data.shape()[1]);

    multi_array<double, 2> data_real(boost::extents[num1][num2]);
    multi_array<double, 2> data_imag(boost::extents[num1][num2]);
    for (int i = 0; i < num1; i++) {
        for (int j = 0; j < num2; j++) {
            data_real[i][j] = data[i][j].real();
            data_imag[i][j] = data[i][j].imag();
        }
    }

    H5::H5File file(filename, H5F_ACC_TRUNC);
    write_hdf5(file, "real", data_real);
    write_hdf5(file, "imag", data_imag);
    write_hdf5(file, label1, grid1);
    write_hdf5(file, label2, grid2);
    file.close();
}

#endif //HDF_WRITER_H
