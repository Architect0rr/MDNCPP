#ifndef __MDN_LOGIC__
#define __MDN_LOGIC__

#include <span>
#include <version>
#include <filesystem>
#include <unordered_set>
#include <chrono>

#include "adios2.h"
#include "mpi.h"
#include "utils.cpp"
#include "constants.cpp"
#include "config.hpp"

#ifdef __cpp_lib_ranges_zip
    #include <ranges>
#else
    #include "zip.cpp"
#endif // !__cpp_lib_ranges_zip

#ifdef __MDN_BUILD_DEBUG_CODE__
    #define __MDN_PROFILING__
#endif // !__MDN_BUILD_DEBUG_CODE__

namespace mdn{

    struct timer{
        uint64_t counter = 0;
        std::chrono::system_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        std::chrono::system_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        void start(){
            t1 = std::chrono::high_resolution_clock::now();
            t2 = t1;
            counter = 0;
        }
        void cutoff(){
            counter += 1;
        }
        void check(){
            t2 = std::chrono::high_resolution_clock::now();
        }
        long double current(){
            check();
            if (counter > 0)
                return std::chrono::duration<long double, std::milli>(t2 - t1).count() / counter;
            else
                return std::chrono::duration<long double, std::milli>(0).count();
        }
        long double count(){
            return std::chrono::duration<long double, std::milli>(t2 - t1).count();
        }
    };

    namespace fs = std::filesystem;
    RETURN_CODES MDN::run(const fs::path &ws, const std::map<fs::path, std::pair<int, int>> &storages, const uint64_t _Natoms){
        adios2::ADIOS adios = adios2::ADIOS(MPI_COMM_SELF);
        adios2::IO dataio   = adios.DeclareIO("DATAWRITER");
        adios2::IO lmpsio   = adios.DeclareIO("LAMMPSReader");
        dataio.SetEngine("BP4");
        lmpsio.SetEngine("BP4");


        adios2::Engine writer = dataio.Open(ws, adios2::Mode::Write, MPI_COMM_SELF);

        adios2::Variable<uint64_t> varFloats = dataio.DefineVariable<uint64_t>("ntimestep");
        adios2::Variable<double>   varTTemp  = dataio.DefineVariable<double>("total_temp", {_Natoms}, {0}, {_Natoms}, adios2::ConstantDims);
        adios2::Variable<double>   varSizes  = dataio.DefineVariable<double>("sizes",      {_Natoms}, {0}, {_Natoms}, adios2::ConstantDims);
        adios2::Variable<double>   varDist   = dataio.DefineVariable<double>("dist",       {_Natoms}, {0}, {_Natoms}, adios2::ConstantDims);
        adios2::Variable<double>   varSCnt   = dataio.DefineVariable<double>("scnt",       {_Natoms}, {0}, {_Natoms}, adios2::ConstantDims);
        adios2::Variable<double>   varTemps  = dataio.DefineVariable<double>("temps",      {_Natoms}, {0}, {_Natoms}, adios2::ConstantDims);

        uint64_t timestep{}, Natoms{};
        double boxxhi{}, boxyhi{}, boxzhi{}, boxxlo{}, boxylo{}, boxzlo{}, Volume{}, rho{};

        int ndim = 3;
        int nprops = 9;

        std::unique_ptr<double[]> Atoms_buf(new double[_Natoms * nprops]);
        std::span<const double> particle_ids(Atoms_buf.get() + 0 * _Natoms, 1 * _Natoms);
        std::span<const double> cluster_ids (Atoms_buf.get() + 1 * _Natoms, 1 * _Natoms);
        std::span<const double> masses      (Atoms_buf.get() + 2 * _Natoms, 1 * _Natoms);
        std::span<const double> velocities  (Atoms_buf.get() + 3 * _Natoms, 3 * _Natoms);
        std::span<const double> xs          (Atoms_buf.get() + 6 * _Natoms, 1 * _Natoms);
        std::span<const double> ys          (Atoms_buf.get() + 7 * _Natoms, 1 * _Natoms);
        std::span<const double> zs          (Atoms_buf.get() + 8 * _Natoms, 1 * _Natoms);

        std::unordered_set<uint64_t> unique_cluster_ids;
        unique_cluster_ids.reserve(_Natoms);
        uint64_t Nclusters{}, size{};
        std::map<uint64_t, std::vector<uint64_t>> particles_by_cluster_id, cluster_ids_by_size, particles_by_size;
        std::vector<uint64_t> sizes_counts(_Natoms + 1, 0UL);
        std::vector<double> kes;
        kes.reserve(_Natoms);
        double total_temp{}, op{};
        std::map<uint64_t, double> temps_by_size;
        // std::map<ul, double> kes_by_size, particle_count_by_size, ndofs_by_size;

        adios2::Variable<uint64_t> varNstep;
        adios2::Variable<uint64_t> varNatoms;
        adios2::Variable<double>   varBoxxhi;
        adios2::Variable<double>   varBoxyhi;
        adios2::Variable<double>   varBoxzhi;
        adios2::Variable<double>   varBoxxlo;
        adios2::Variable<double>   varBoxylo;
        adios2::Variable<double>   varBoxzlo;
        adios2::Variable<double>   varAtoms;
        // id c_clusters mass vx vy vz x y z

        #ifdef __MDN_PROFILING__
            if (rank == 0)
                std::cout << "Profiling enabled" << std::endl;
            timer global_timer;
            timer per_storage_timer;
            timer per_step_timer;
            global_timer.start();
            per_storage_timer.start();
        #endif // !__MDN_PROFILING__

        for (const auto &[storage, steps] : storages){

            try{
                adios2::Engine reader = lmpsio.Open(ws, adios2::Mode::Read, MPI_COMM_SELF);
                // if (rank == 0)
                //     std::cout << "Storage open, scrolling" << std::endl;

                // varNstep  = lmpsio.InquireVariable<uint64_t>(std::string(lcf::timestep));
                // varNatoms = lmpsio.InquireVariable<uint64_t>(std::string(lcf::natoms  ));
                // varBoxxhi = lmpsio.InquireVariable<double>  (std::string(lcf::boxxhi  ));
                // varBoxyhi = lmpsio.InquireVariable<double>  (std::string(lcf::boxyhi  ));
                // varBoxzhi = lmpsio.InquireVariable<double>  (std::string(lcf::boxzhi  ));
                // varBoxxlo = lmpsio.InquireVariable<double>  (std::string(lcf::boxxlo  ));
                // varBoxylo = lmpsio.InquireVariable<double>  (std::string(lcf::boxylo  ));
                // varBoxzlo = lmpsio.InquireVariable<double>  (std::string(lcf::boxzlo  ));
                // varAtoms  = lmpsio.InquireVariable<double>  (std::string(lcf::atoms   ));

                uint64_t currentStep = 0;
                while (currentStep != steps.first)
                {
                    // if (rank == 0)
                    //     std::cout << "Enter loop" << std::endl;
                    if (reader.BeginStep() == adios2::StepStatus::EndOfStream)
                    {
                        logger.error("End of stream happened while scrolling to begin step");
                        logger.error("Storage: {}, begin: {}, end: {}", storage.string(), steps.first, steps.second);
                        throw std::logic_error("End of stream happened while scrolling to begin step");
                    }
                    currentStep = reader.CurrentStep();
                    reader.EndStep();
                    // if (rank == 0)
                    //     std::cout << "Current step: " << currentStep << std::endl;
                }
                // if (rank == 0)
                //     std::cout << "End loop" << std::endl;
                #ifdef __MDN_PROFILING__
                    // if (rank == 0)
                    //     std::cout << "Start timer" << std::endl;
                    per_step_timer.start();
                    // if (rank == 0)
                    //     std::cout << "Timer started" << std::endl;
                #endif // !__MDN_PROFILING__
                if (rank == 0)
                    std::cout << "Starting calculations" << std::endl;
                while (currentStep != steps.first + steps.second)
                {
                    if (rank == 0)
                        std::cout << "In loop" << std::endl;
                    if (reader.BeginStep() == adios2::StepStatus::EndOfStream)
                        if (rank == 0)
                            std::cout << "EOS" << std::endl;
                        break;
                    currentStep = reader.CurrentStep();
                    if (rank == 0)
                        std::cout << "Current step: " << currentStep << std::endl;

                    varNstep  = lmpsio.InquireVariable<uint64_t>(std::string(lcf::timestep));
                    varNatoms = lmpsio.InquireVariable<uint64_t>(std::string(lcf::natoms  ));
                    varBoxxhi = lmpsio.InquireVariable<double>  (std::string(lcf::boxxhi  ));
                    varBoxyhi = lmpsio.InquireVariable<double>  (std::string(lcf::boxyhi  ));
                    varBoxzhi = lmpsio.InquireVariable<double>  (std::string(lcf::boxzhi  ));
                    varBoxxlo = lmpsio.InquireVariable<double>  (std::string(lcf::boxxlo  ));
                    varBoxylo = lmpsio.InquireVariable<double>  (std::string(lcf::boxylo  ));
                    varBoxzlo = lmpsio.InquireVariable<double>  (std::string(lcf::boxzlo  ));
                    varAtoms  = lmpsio.InquireVariable<double>  (std::string(lcf::atoms   ));
                    if (!(varNstep && varNatoms && varBoxxhi && varBoxyhi && varBoxzhi && varBoxxlo && varBoxylo && varBoxzlo && varAtoms)){

                    }
                    if (rank == 0)
                        std::cout << "Inquired variables" << std::endl;
                    reader.Get(varNstep, timestep);
                    reader.Get(varNatoms, Natoms);
                    reader.Get(varBoxxhi, boxxhi);
                    reader.Get(varBoxyhi, boxyhi);
                    reader.Get(varBoxzhi, boxzhi);
                    reader.Get(varBoxxlo, boxxlo);
                    reader.Get(varBoxylo, boxylo);
                    reader.Get(varBoxzlo, boxzlo);
                    if (rank == 0)
                        std::cout << "Got variables" << std::endl;
                    Volume = abs((boxxhi - boxxlo) * (boxyhi - boxylo) * (boxzhi - boxzlo));
                    rho = Natoms / Volume;

                    for (const uint64_t i : cluster_ids)
                    {
                        unique_cluster_ids.insert(i);
                    }
                    Nclusters = unique_cluster_ids.size();

                    size = 0;
                    for (const uint64_t &i : unique_cluster_ids)
                    {
                        for (uint64_t j = 0; j < Natoms; ++j)
                        {
                            if (cluster_ids[j] == i)
                                particles_by_cluster_id[i].push_back(j);
                        }
                        size = particles_by_cluster_id[i].size();
                        cluster_ids_by_size[size].push_back(i);
                        particles_by_size[size].insert(particles_by_size[size].end(), particles_by_cluster_id[i].begin(), particles_by_cluster_id[i].end());
                    }

                    for (const auto &[k, v] : cluster_ids_by_size)
                    {
                        sizes_counts[k] = v.size();
                    }

                    total_temp = 0;

                    #ifndef __cpp_lib_ranges_zip
                        for (const auto &[vel, mass] : zip<const std::span<const double> &, const std::span<const double> &>(velocities, masses))
                    #else
                        for (const auto &[vel, mass] : std::ranges::views::zip(velocities, masses))
                    #endif // !__cpp_lib_ranges_zip
                    {
                            kes.emplace_back(mass * vel / 2);
                            total_temp += kes.back();
                    }
                    total_temp = std::pow(total_temp / ((Natoms - 1) * ndim), 2);

                    for (const auto &[k, v] : particles_by_size)
                    {
                        op = 0;
                        for (const uint64_t& i : v)
                        {
                            op += kes[i];
                        }
                        // kes_by_size[k] = op;
                        // particle_count_by_size = v.size();
                        // ndofs_by_size[k] = (v.size() - 1) * ndim;
                        temps_by_size.emplace(k, std::pow(op / ((v.size() - 1) * ndim), 2));
                    }


                    reader.EndStep();

                    unique_cluster_ids.clear();
                    particles_by_cluster_id.clear();
                    cluster_ids_by_size.clear();
                    particles_by_size.clear();
                    // sizes_counts.clear();
                    std::fill(sizes_counts.begin(), sizes_counts.end(), 0);
                    kes.clear();
                    temps_by_size.clear();

                    #ifdef __MDN_PROFILING__
                        per_step_timer.cutoff();
                        // if (global_timer.count() > 1000*60)
                        std::cout << "Profiling: " << per_step_timer.current() << " ms per step, " << per_storage_timer.current() << " ms per storage" << std::endl;
                        logger.info("Profiling: {} ms per step, {} ms per storage", per_step_timer.current(), per_storage_timer.current());
                    #endif // !__MDN_PROFILING__

                }
                if (rank == 0)
                    std::cout << "End storage" << std::endl;

                reader.Close();
            }catch (std::logic_error& e){
                logger.error("Some std::logic_error happened on storage {}, steps: ()", storage.string(), steps.first, steps.second);
                logger.error(e.what());
                logger.error("Assuming we can go to the next storage");
            }catch(std::exception& e){
                logger.error("Some std::exception happened on storage {}, steps: ()", storage.string(), steps.first, steps.second);
                logger.error(e.what());
                logger.error("Assuming we can go to the next storage");
            }catch (...){
                logger.error("Some error happened on storage {}, steps: ()", storage.string(), steps.first, steps.second);
                logger.error("Error is not inherits std::exception, so exiting... Probably the whole program may be aborted");
                throw;
            }
            #ifdef __MDN_PROFILING__
                per_storage_timer.cutoff();
            #endif // !__MDN_PROFILING__
        }
        #ifdef __MDN_PROFILING__
            per_step_timer.check();
            per_storage_timer.check();
            global_timer.check();
            global_timer.cutoff();
            std::cout << "Profiling: " << per_step_timer.current() << " ms per step, " << per_storage_timer.current() << " ms per storage" << std::endl;
            std::cout << "Profiling: " << global_timer.current() << " ms total wall time" << std::endl;
            logger.info("Final profiling results: {} ms per step, {} ms per storage", per_step_timer.current(), per_storage_timer.current());
            logger.info("Profiling: {} ms total wall time", global_timer.current());
        #endif // !__MDN_PROFILING__

        return RETURN_CODES::OK;
    }
} // namespace

#endif // !__MDN_LOGIC__
