#ifndef __MDN_constants__
#define __MDN_constants__

#include <string>

namespace mdn{

    namespace cs
    {
        constexpr const int mpi_root = 0;
    }
    namespace files
    {
        constexpr const std::string_view data("data.json");
        constexpr const std::string_view output_datafile_basename("data.bp");
        constexpr const std::string_view distribution_cache("distribution.cache");
        constexpr const std::string_view storages_cache("storages.cache");
        constexpr const std::string_view matrix_csv("matrix.csv");
        constexpr const std::string_view temps_csv("temps.csv");
    }
    namespace folders
    {
        constexpr const std::string_view log("logs/MDNCPP");
        constexpr const std::string_view post_process_folder("post_processing");
        constexpr const std::string_view cache("cache");
    }
    namespace fields
    {
        constexpr const std::string_view storages("storages");
        constexpr const std::string_view Natoms("N_atoms");
        constexpr const std::string_view Volume("Volume");
        constexpr const std::string_view Dims("dimensions");
        constexpr const std::string_view data_processing_folder("post_process_folder");
    }
    namespace lcf
    {
        constexpr const std::string_view timestep("ntimestep");
        constexpr const std::string_view atoms("atoms");
        constexpr const std::string_view natoms("natoms");
        constexpr const std::string_view boxxhi("boxxhi");
        constexpr const std::string_view boxyhi("boxyhi");
        constexpr const std::string_view boxzhi("boxzhi");
        constexpr const std::string_view boxxlo("boxxlo");
        constexpr const std::string_view boxylo("boxylo");
        constexpr const std::string_view boxzlo("boxzlo");
    }
} // namespace mdn

#endif // !__MDN_constants__
