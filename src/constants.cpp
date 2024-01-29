#ifndef __MDN_constants__
#define __MDN_constants__

#include <string>

namespace mdn::cs
{
    constexpr const int mpi_root = 0;
}
namespace mdn::files
{
    constexpr const std::string_view data("data.json");
}
namespace mdn::fields
{
    constexpr const std::string_view storages("storages");
    constexpr const std::string_view Natoms("N_atoms");
    constexpr const std::string_view Volume("Volume");
    constexpr const std::string_view Dims("dimensions");
    constexpr const std::string_view data_processing_folder("post_process_folder");
}
namespace mdn::lcf
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
#endif // !__MDN_constants__