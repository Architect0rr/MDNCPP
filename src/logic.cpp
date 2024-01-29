#ifndef __MDN_LOGIC__
#define __MDN_LOGIC__

#include "adios2.h"
#include <xtensor/xarray.hpp>
#include <xtensor/xadapt.hpp>

#include "mpi.h"

#include "utils.cpp"
#include "constants.cpp"

#include <filesystem>

namespace mdn::logic
{
    namespace fs = std::filesystem;
    RETURN_CODES run(utils::MC &sts, const fs::path &ws, const std::map<fs::path, std::pair<int, int>> &storages, const uint64_t _Natoms)
    {
        adios2::ADIOS adios = adios2::ADIOS(MPI_COMM_SELF);
        adios2::IO dataio = adios.DeclareIO("DATAWRITER");
        adios2::IO lmpio = adios.DeclareIO("LAMMPSReader");
        dataio.SetEngine("BP5");
        lmpio.SetEngine("BP5");

        adios2::Engine writer = dataio.Open(ws, adios2::Mode::Write, MPI_COMM_SELF);

        adios2::Variable<uint64_t> varFloats = dataio.DefineVariable<uint64_t>("ntimestep");
        adios2::Variable<double> varTTemp = dataio.DefineVariable<double>("total_temp", {_Natoms}, {0}, {_Natoms}, adios2::ConstantDims);
        adios2::Variable<double> varSizes = dataio.DefineVariable<double>("sizes", {_Natoms}, {0}, {_Natoms}, adios2::ConstantDims);
        adios2::Variable<double> varDist = dataio.DefineVariable<double>("dist", {_Natoms}, {0}, {_Natoms}, adios2::ConstantDims);
        adios2::Variable<double> varSCnt = dataio.DefineVariable<double>("scnt", {_Natoms}, {0}, {_Natoms}, adios2::ConstantDims);
        adios2::Variable<double> varTemps = dataio.DefineVariable<double>("temps", {_Natoms}, {0}, {_Natoms}, adios2::ConstantDims);

        uint64_t timestep{};
        uint64_t Natoms{};
        double boxxhi{};
        double boxyhi{};
        double boxzhi{};
        double boxxlo{};
        double boxylo{};
        double boxzlo{};
        double Volume{};
        double rho{};

        size_t ndim = 3;
        size_t nprops = 9;
        std::shared_ptr<double[]> Atoms_buf(new double[_Natoms * nprops]);
        xt::xarray<double>::shape_type lineshape = {_Natoms * nprops};
        xt::xarray<double>::shape_type rightshape = {_Natoms, nprops};

        for (const auto &[storage, steps] : storages)
        {
            adios2::Engine reader = lmpio.Open(ws, adios2::Mode::Read, MPI_COMM_SELF);
            uint64_t currentStep = 0;
            while (steps.first != currentStep && reader.BeginStep() != adios2::StepStatus::EndOfStream)
            {
                currentStep = reader.CurrentStep();
                reader.EndStep();
            }
            adios2::Variable<uint64_t> varNstep = lmpio.InquireVariable<uint64_t>(std::string(lcf::timestep));
            adios2::Variable<uint64_t> varNatoms = lmpio.InquireVariable<uint64_t>(std::string(lcf::natoms));
            adios2::Variable<double> varBoxxhi = lmpio.InquireVariable<double>(std::string(lcf::boxxhi));
            adios2::Variable<double> varBoxyhi = lmpio.InquireVariable<double>(std::string(lcf::boxyhi));
            adios2::Variable<double> varBoxzhi = lmpio.InquireVariable<double>(std::string(lcf::boxzhi));
            adios2::Variable<double> varBoxxlo = lmpio.InquireVariable<double>(std::string(lcf::boxxlo));
            adios2::Variable<double> varBoxylo = lmpio.InquireVariable<double>(std::string(lcf::boxylo));
            adios2::Variable<double> varBoxzlo = lmpio.InquireVariable<double>(std::string(lcf::boxzlo));
            // id c_clusters mass vx vy vz x y z
            adios2::Variable<double> varAtoms = lmpio.InquireVariable<double>(std::string(lcf::atoms));

            reader.Get(varNatoms, Natoms);
            reader.Get(varBoxxhi, boxxhi);
            reader.Get(varBoxyhi, boxyhi);
            reader.Get(varBoxzhi, boxzhi);
            reader.Get(varBoxxlo, boxxlo);
            reader.Get(varBoxylo, boxylo);
            reader.Get(varBoxzlo, boxzlo);

            Volume = abs((boxxhi - boxxlo) * (boxyhi - boxylo) * (boxzhi - boxzlo));
            rho = Natoms / Volume;

            auto Atoms = xt::adapt_smart_ptr(Atoms_buf, rightshape);

            reader.Close();
        }

        return RETURN_CODES::OK;
    }
}

#endif // !__MDN_LOGIC__