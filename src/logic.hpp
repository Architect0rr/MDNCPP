#ifndef __MDN_LOGIC__
#define __MDN_LOGIC__

#include <span>
#include <version>
#include <filesystem>
#include <unordered_set>
#include <chrono>
#include <numbers>

#include "adios2.h"
#include "mpi.h"
#include "enums.hpp"
#include "constants.hpp"
#include "config.hpp"

#ifdef __cpp_lib_ranges_zip
    #include <ranges>
#else
    #include "zip.hpp"
#endif // !__cpp_lib_ranges_zip

#ifdef __MDN_BUILD_DEBUG_CODE__
    #define __MDN_PROFILING__
#endif // !__MDN_BUILD_DEBUG_CODE__

#define __CALC_ENTHROPY__
#define __KE_PE_PRESENT__
#ifdef __CALC_ENTHROPY__
    #ifndef __KE_PE_PRESENT__
        #define __CALC_PE__
    #endif // !__KE_PE_PRESENT__
#endif // __CALC_ENTHROPY__


namespace mdn{

    template <typename T, uint64_t N>
    constexpr std::array<T, N> _linspace(T start, T stop){
        std::array<T, N> a;
        T dx = (stop - start) / N;
        for (auto& el : a){
            el = start;
            start += dx;
        }
        return a;
    }

    template <typename T, uint64_t N>
    constexpr std::array<T, N-1> _centers(const std::array<T, N>& arr){
        std::array<T, N-1> a;
        for (size_t i = 0; i < N-1; ++i) a[i] = (arr[i] + arr[i+1])/2;
        return a;
    }

    template <typename T, uint64_t N>
    constexpr std::array<T, N-1> _sl_volumes(const std::array<T, N>& arr, const uint64_t Np_bins, const uint64_t Nt_bins){
        std::array<T, N-1> lv;
        for (size_t i = 0; i < N-1; ++i) lv[i] = 4 * std::numbers::pi * (std::pow(arr[i+1], 3) - std::pow(arr[i], 3)) / 3 / Np_bins / Nt_bins;
        return lv;
    }

    namespace fs = std::filesystem;
    RETURN_CODES MDN::run(const fs::path &ws, const std::map<fs::path, std::pair<int, int>> &storages, const uint64_t _Natoms, uint64_t& max_cluster_size){
        adios2::ADIOS adios = adios2::ADIOS("adios2_BP4_config.xml", MPI_COMM_SELF);
        adios2::IO dataio   = adios.DeclareIO("DataWriter");
        adios2::IO lmpsio   = adios.DeclareIO("LAMMPSReader");
        // dataio.SetEngine("BP4");
        // lmpsio.SetEngine("BP4");

        adios2::Engine writer = dataio.Open(ws, adios2::Mode::Write, MPI_COMM_SELF);

        adios2::Variable<uint64_t> WvarNstep  = dataio.DefineVariable<uint64_t>("ntimestep");
        adios2::Variable<uint64_t> WvarN      = dataio.DefineVariable<uint64_t>("Natoms");
        adios2::Variable<uint64_t> WvarDist   = dataio.DefineVariable<uint64_t>("dist",     {_Natoms + 1}, {0}, {_Natoms + 1}, adios2::ConstantDims);
        adios2::Variable<double>   WvarTemps  = dataio.DefineVariable<double>  ("temps",    {_Natoms + 1}, {0}, {_Natoms + 1}, adios2::ConstantDims);
        adios2::Variable<double>   WvarVol    = dataio.DefineVariable<double>  ("volume");
        adios2::Variable<double>   WvarTTemp  = dataio.DefineVariable<double>  ("total_temp");
        #ifdef __CALC_ENTHROPY__
        adios2::Variable<double>   WvarEnthr  = dataio.DefineVariable<double>  ("enthr",    {_Natoms + 1}, {0}, {_Natoms + 1}, adios2::ConstantDims);
        adios2::Variable<double>   WvarEnerg  = dataio.DefineVariable<double>  ("energ",    {_Natoms + 1}, {0}, {_Natoms + 1}, adios2::ConstantDims);
        #endif // __CALC_ENTHROPY__

        adios2::Variable<uint64_t> varNstep;
        adios2::Variable<uint64_t> varNatoms;
        adios2::Variable<double>   varBoxxhi;
        adios2::Variable<double>   varBoxyhi;
        adios2::Variable<double>   varBoxzhi;
        adios2::Variable<double>   varBoxxlo;
        adios2::Variable<double>   varBoxylo;
        adios2::Variable<double>   varBoxzlo;
        adios2::Variable<double>   varAtoms;

        constexpr int ndim = 3;
        #ifdef __CALC_ENTHROPY__
            #ifdef __KE_PE_PRESENT__
                constexpr int nprops = 11;
            #else
                constexpr int nprops = 9;
            #endif // __KE_PE_PRESENT__
        #else
            constexpr int nprops = 6;
        #endif // __CALC_ENTHROPY__
        constexpr double rcut = 2.5;

        #ifdef __CALC_ENTHROPY__
            constexpr size_t Nr_bins = 10;
            constexpr size_t Np_bins = 36;
            constexpr size_t Nt_bins = 18;
            constexpr double rho_normalize = 4 * std::numbers::pi * std::pow(rcut, 3) / 3;
            constexpr double bin_r_w = (rcut - 0) / Nr_bins;
            constexpr double bin_p_w = (2 * std::numbers::pi) / Np_bins;
            constexpr double bin_t_w = (std::numbers::pi - 0) / Nt_bins;

            constexpr std::array<double, Nr_bins + 1> r_mesh = _linspace<double, Nr_bins + 1>(0.0, rcut);
            constexpr std::array<double, Nr_bins> rbs = _centers(r_mesh);
            constexpr std::array<double, Nr_bins> lv  = _sl_volumes(r_mesh, Np_bins, Nt_bins);
        #endif // __CALC_ENTHROPY__

        std::unique_ptr<double[]> Atoms_buf(new double[_Natoms * nprops]);
        std::span<const double> particle_ids(Atoms_buf.get() + 0 * _Natoms, 1 * _Natoms);
        std::span<const double> cluster_ids (Atoms_buf.get() + 1 * _Natoms, 1 * _Natoms);
        std::span<const double> masses      (Atoms_buf.get() + 2 * _Natoms, 1 * _Natoms);
        std::span<const double> velocities  (Atoms_buf.get() + 3 * _Natoms, 3 * _Natoms);
        #ifdef __CALC_ENTHROPY__
            std::span<const double> xs          (Atoms_buf.get() + 6 * _Natoms, 1 * _Natoms);
            std::span<const double> ys          (Atoms_buf.get() + 7 * _Natoms, 1 * _Natoms);
            std::span<const double> zs          (Atoms_buf.get() + 8 * _Natoms, 1 * _Natoms);
        #endif // __CALC_ENTHROPY__
        #ifdef __KE_PE_PRESENT__
            std::span<const double> kes         (Atoms_buf.get() + 9 * _Natoms, 1 * _Natoms);
            std::span<const double> pes         (Atoms_buf.get() + 10 * _Natoms, 1 * _Natoms);
        #endif // __KE_PE_PRESENT__

        // Main algorithm containers
        uint64_t timestep{}, Natoms{}, Nclusters{}, size{};
        double boxxhi{}, boxyhi{}, boxzhi{}, boxxlo{}, boxylo{}, boxzlo{}, Volume{}, total_temp{}, ke{}, pe{};
        std::unordered_set<uint64_t> unique_cluster_ids;
        unique_cluster_ids.reserve(_Natoms);
        std::map<uint64_t, std::vector<uint64_t>> particles_by_cluster_id, cluster_ids_by_size, particles_by_size;
        std::vector<uint64_t> sizes_counts(_Natoms + 1, 0UL);

        #ifndef __KE_PE_PRESENT__
            std::vector<double> kes;
            kes.reserve(_Natoms);
            #ifdef __CALC_PE__
                std::vector<double> pes(_Natoms, 0.0);
                double dist{}, en{};
            #endif // !__CALC_PE__
        #endif // !__KE_PE_PRESENT__

        // std::map<uint64_t, double> temps_by_size;
        std::vector<double> temps_by_size(_Natoms + 1, 0.0);
        // std::map<uint64_t, double> particle_count_by_size, ndofs_by_size;

        // Enthropy containers
        #ifdef __CALC_ENTHROPY__
            std::vector<std::array<double, ndim>> pts; // cleared
            pts.reserve(_Natoms);
            std::vector<std::array<double, ndim>> _pts; // cleared
            _pts.reserve(_Natoms);
            std::vector<double> rho_a; // not required
            rho_a.resize(_Natoms, 0.0);

            std::vector<double> rvi(_Natoms, 0.0);   // cleaning not required
            size_t r_bin_n{}, p_bin_n{}, t_bin_n{};
            double rho{}, s{}, phi{}, theta{};
            std::array<std::vector<int>, Nr_bins+1> r_bin; // cleared
            std::array<std::array<std::array<std::vector<int>, Np_bins+1>, Nt_bins+1>, Nr_bins+1> ang_bin; // cleared
            for (size_t i = 0; i < Nr_bins+1; ++i){
                r_bin[i].reserve(_Natoms);
                for (size_t j = 0; j < Nt_bins+1; ++j)
                    for (size_t k = 0; k < Np_bins+1; ++k)
                        ang_bin[i][j][k].reserve(_Natoms);
            }
            std::array<double, Nr_bins> rs;  // cleaning not required
            std::array<std::array<std::array<double, Np_bins>, Nt_bins>, Nr_bins> ang_s;  // cleaning not required
            std::array<double, Nr_bins> enth_ang({0.0}); // cleaning not required
            double r_sum{}, t_sum{}, enth{};
            std::vector<double> enth_by_size(_Natoms + 1, 0.0); // cleared
            // std::map<uint64_t, double> kes_by_size; // cleared
            #ifdef __KE_PE_PRESENT__
                // std::vector<double> f_energy(_Natoms + 1, 0.0);
                // std::vector<double> kes_by_size(_Natoms + 1, 0.0); // cleared
                // std::vector<double> pes_by_size(_Natoms + 1, 0.0); // cleared
                std::vector<double> eng_by_size(_Natoms + 1, 0.0); // cleared
            #endif // __KE_PE_PRESENT__

            // std::map<size_t, std::vector<std::pair<std::array<double, 3>, std::array<double, Nr_bins>>>> stats;
            // std::array<double, Nr_bins> ang_int;
            // double rr_integral=0;
        #endif // __CALC_ENTHROPY__

        // id c_clusters mass vx vy vz [x y z [ke pe]]

        #ifdef __MDN_PROFILING__
            timer global_timer;
            timer per_storage_timer;
            timer per_step_timer;

            global_timer.start();
            per_storage_timer.start();
        #endif // !__MDN_PROFILING__

        for (const auto &[storage, steps] : storages){
            try{
                adios2::Engine reader = lmpsio.Open(storage, adios2::Mode::Read, MPI_COMM_SELF);

                uint64_t currentStep = 0;
                while (currentStep != steps.first)
                {
                    if (reader.BeginStep() == adios2::StepStatus::EndOfStream)
                    {
                        logger.error("End of stream happened while scrolling to begin step");
                        logger.error("Storage: {}, begin: {}, end: {}", storage.string(), steps.first, steps.second);
                        throw std::logic_error("End of stream happened while scrolling to begin step");
                    }
                    currentStep = reader.CurrentStep();
                    reader.EndStep();
                }

                #ifdef __MDN_PROFILING__
                    per_step_timer.start();
                #endif // !__MDN_PROFILING__

                while (currentStep != steps.first + steps.second)
                {
                    if (reader.BeginStep() == adios2::StepStatus::EndOfStream)
                        break;

                    currentStep = reader.CurrentStep();

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
                        reader.EndStep();
                        logger.error("Error on read variables at step {}", currentStep);
                        continue;
                    }

                    reader.Get(varNstep, timestep);
                    reader.Get(varNatoms, Natoms);
                    reader.Get(varBoxxhi, boxxhi);
                    reader.Get(varBoxyhi, boxyhi);
                    reader.Get(varBoxzhi, boxzhi);
                    reader.Get(varBoxxlo, boxxlo);
                    reader.Get(varBoxylo, boxylo);
                    reader.Get(varBoxzlo, boxzlo);
                    reader.Get(varAtoms, Atoms_buf.get());
                    reader.PerformGets();
                    reader.EndStep();

                    Volume = abs((boxxhi - boxxlo) * (boxyhi - boxylo) * (boxzhi - boxzlo));

                    for (uint64_t j = 0; j < Natoms; ++j) particles_by_cluster_id[cluster_ids[j]].emplace_back(j);

                    for (const uint64_t& i : cluster_ids) unique_cluster_ids.emplace(i);

                    Nclusters = unique_cluster_ids.size();


                    size = 0;
                    for (const uint64_t &i : unique_cluster_ids)
                    {
                        size = particles_by_cluster_id[i].size();
                        cluster_ids_by_size[size].push_back(i);
                        particles_by_size[size].insert(particles_by_size[size].end(), particles_by_cluster_id[i].begin(), particles_by_cluster_id[i].end());
                    }

                    for (const auto &[k, v] : cluster_ids_by_size)
                    {
                        sizes_counts[k] = v.size();
                        if (k > max_cluster_size) max_cluster_size = k;
                    }

                    #ifdef __CALC_PE__
                        for (size_t i = 0; i < Natoms; ++i)
                            for (size_t j = i; j < Natoms; ++j)
                            {
                                dist = std::sqrt(std::pow(xs[i] - xs[j], 2) + std::pow(ys[i] - ys[j], 2) + std::pow(zs[i] - zs[j], 2));
                                if (dist < rcut){
                                    en = 4*(std::pow(dist, -12) - std::pow(dist, -6));
                                    pes[i] += en;
                                    pes[j] += en;
                                }
                            }
                    #endif // __CALC_PE__

                    total_temp = 0;
                    #ifndef __KE_PE_PRESENT__
                        #ifndef __cpp_lib_ranges_zip
                            for (const auto &[vel, mass] : zip<const std::span<const double> &, const std::span<const double> &>(velocities, masses))
                        #else
                            for (const auto &[vel, mass] : std::ranges::views::zip(velocities, masses))
                        #endif // !__cpp_lib_ranges_zip
                        {
                                kes.emplace_back(mass * vel / 2);
                                total_temp += kes.back();
                        }
                    #else
                        for (const double& ke: kes) total_temp += ke;
                    #endif // !__KE_PE_PRESENT__
                    total_temp = 2 * total_temp / ((Natoms - 1) * ndim);

                    for (const auto &[size, v] : particles_by_size)
                    {
                        ke = 0;
                        #ifdef __KE_PE_PRESENT__
                            pe = 0;
                        #endif // __KE_PE_PRESENT__
                        for (const uint64_t& i : v)
                        {
                            ke += kes[i];
                            #ifdef __CALC_ENTHROPY__
                                pe += pes[i];
                            #endif // __CALC_ENTHROPY__
                        }
                        #ifdef __CALC_ENTHROPY__
                            eng_by_size[size] = (ke + pe) / v.size();
                            // kes_by_size[size] = ke / v.size();
                            // pes_by_size[size] = pe / v.size();
                        #endif // __CALC_ENTHROPY__

                        // particle_count_by_size = v.size();
                        // ndofs_by_size[k] = (v.size() - 1) * ndim;
                        temps_by_size[size] = 2 * ke / ((v.size() - 1) * ndim);
                        // temps_by_size.emplace(k, std::pow(ke / ((v.size() - 1) * ndim), 2));
                    }

                    #ifdef __CALC_ENTHROPY__
                        for (const auto&[size, cl_ids] : cluster_ids_by_size){
                            for (const uint64_t& cl_id: cl_ids){

                                pts.clear();
                                for (const uint64_t& i : particles_by_cluster_id[cl_id]){
                                    pts[i] = {xs[i], ys[i], zs[i]};
                                }

                                for (size_t i = 0; i < Nr_bins+1; ++i){
                                    r_bin[i].clear();
                                    r_bin[i].resize(size, 0);
                                    for (size_t j = 0; j < Nt_bins+1; ++j)
                                        for (size_t k = 0; k < Np_bins+1; ++k){
                                            ang_bin[i][j][k].clear();
                                            ang_bin[i][j][k].resize(size, 0);
                                        }
                                }

                                size_t i = 0;
                                for (size_t id = 0; id < size; ++id){
                                    _pts.clear();
                                    i = 0;

                                    for (size_t _id = 0; _id < size; ++_id){
                                        if (_id == id) continue;

                                        rvi[i] = std::sqrt(std::pow(pts[_id][0] - pts[id][0], 2) + std::pow(pts[_id][1] - pts[id][1], 2) + std::pow(pts[_id][2] - pts[id][2], 2));

                                        if (rvi[i] < rcut){
                                            _pts.push_back(pts[_id]);
                                            ++i;
                                        }
                                    }

                                    rho_a[id] = static_cast<double>(_pts.size());

                                    for (size_t _id = 0; _id < _pts.size(); ++_id){
                                        r_bin_n = static_cast<size_t>(rvi[_id] / bin_r_w);

                                        phi = std::atan2(_pts[_id][1] - pts[id][1], _pts[_id][0] - pts[id][0]) + std::numbers::pi;
                                        // if (phi < 0) phi += 2 * std::numbers::pi;
                                        p_bin_n = static_cast<size_t>(phi / bin_p_w);

                                        theta = std::atan((_pts[_id][2] - pts[id][2])/(_pts[_id][1] - pts[id][1]));
                                        // if (theta < - std::numbers::pi / 2) theta = - (std::numbers::pi + theta);
                                        // if (theta >  std::numbers::pi / 2) theta = std::numbers::pi - theta;
                                        t_bin_n = static_cast<size_t>((theta + std::numbers::pi / 2) / bin_t_w);

                                        ++ang_bin[r_bin_n][t_bin_n][p_bin_n][id];
                                        ++r_bin[r_bin_n][id];
                                    }
                                }

                                for (size_t i = 0; i < size; ++i)
                                {
                                    r_bin[Nr_bins-1][i] += r_bin[Nr_bins][i];
                                }
                                for (size_t i = 0; i < Nr_bins; ++i)
                                    for (size_t j = 0; j < Nt_bins; ++j)
                                        for (size_t k = 0; k < size; ++k)
                                            ang_bin[i][j][Np_bins - 1][k] += ang_bin[i][j][Np_bins][k];
                                for (size_t i = 0; i < Nr_bins; ++i)
                                    for (size_t j = 0; j < Np_bins; ++j)
                                        for (size_t k = 0; k < size; ++k)
                                            ang_bin[i][Nt_bins - 1][j][k] += ang_bin[i][Nt_bins][j][k];
                                for (size_t i = 0; i < Nt_bins; ++i)
                                    for (size_t j = 0; j < Np_bins; ++j)
                                        for (size_t k = 0; k < size; ++k)
                                            ang_bin[Nr_bins - 1][i][j][k] += ang_bin[Nr_bins][i][j][k];

                                rho = 0;
                                for (size_t i = 0; i < size; ++i)
                                    rho += rho_a[i];
                                rho /= size * rho_normalize;

                                s = 0;
                                for (size_t i = 0; i < Nr_bins; ++i)
                                {
                                    s = 0;
                                    for (size_t j = 0; j < size; ++j) s += r_bin[i][j];
                                    rs[i] = s / (size * lv[i] * rho);
                                    // rs[i] = s / (size * 4 * std::numbers::pi * std::pow(rbs[i], 2) * bin_r_w * rho);
                                }

                                for (size_t i = 0; i < Nr_bins; ++i){
                                    for (size_t j = 0; j < Nt_bins; ++j)
                                        for (size_t k = 0; k < Np_bins; ++k){
                                            s = 0;
                                            for (size_t l = 0; l < size; ++l) s += ang_bin[i][j][k][l];
                                            ang_s[i][j][k] = s / (size * lv[i] * rho);
                                            // ang_s[i][j][k] = s / (size * 4 * std::numbers::pi * std::pow(rbs[i], 2) * bin_r_w * rho);
                                        }
                                }

                                for (size_t i = 0; i < Nr_bins; ++i){
                                    r_sum = 0;
                                    for (size_t j = 0; j < Nt_bins; ++j){
                                        t_sum = 0;
                                        // ternary avoid nan in 0*ln(0)
                                        for (size_t k = 0; k < Np_bins; ++k) t_sum += (ang_s[i][j][k]>0) ? ang_s[i][j][k] * std::log(ang_s[i][j][k]) * bin_p_w : 0;
                                        r_sum += t_sum * bin_t_w;
                                    }
                                    enth_ang[i] = -r_sum / (4 * std::numbers::pi);
                                }

                                enth = 0;
                                for (size_t i = 0; i < Nr_bins; ++i)
                                    enth += rs[i] * enth_ang[i] * std::pow(rbs[i], 2);
                                enth *= 4 * std::numbers::pi * bin_r_w;
                                enth_by_size[size] += enth;
                            }
                            enth_by_size[size] /= cl_ids.size();
                            // f_energy[size] = kes_by_size[size] + pes_by_size[size] - enth_by_size[size]*temps_by_size[size];
                            // kes_by_size[size] + pes_by_size[size] - enth_by_size[size]*temps_by_size[size];
                        }
                    #endif // __CALC_ENTHROPY__

                    writer.BeginStep();
                    writer.Put(WvarNstep, timestep);
                    writer.Put(WvarDist,   sizes_counts.data());
                    writer.Put(WvarVol,   Volume);
                    writer.Put(WvarN,     Natoms);
                    writer.Put(WvarTTemp, total_temp);
                    writer.Put(WvarTemps,  temps_by_size.data());
                    #ifdef __CALC_ENTHROPY__
                        writer.Put(WvarEnthr,  enth_by_size.data());
                        writer.Put(WvarEnerg,  eng_by_size.data());
                    #endif // __CALC_ENTHROPY__
                    writer.PerformPuts();
                    // writer.PerformDataWrite();
                    writer.EndStep();

                    unique_cluster_ids.clear();
                    particles_by_cluster_id.clear();
                    cluster_ids_by_size.clear();
                    particles_by_size.clear();
                    std::fill(sizes_counts.begin(), sizes_counts.end(), 0);
                    std::fill(temps_by_size.begin(), temps_by_size.end(), 0.0);

                    #ifdef __CALC_ENTHROPY__
                        std::fill(enth_by_size.begin(), enth_by_size.end(), 0.0);
                        std::fill(eng_by_size.begin(), eng_by_size.end(), 0.0);
                        // std::fill(kes_by_size.begin(), kes_by_size.end(), 0.0);
                        // std::fill(pes_by_size.begin(), pes_by_size.end(), 0.0);
                        #ifndef __KE_PE_PRESENT__
                            kes.clear();
                            pes.clear();
                        // #else
                        //     std::fill(f_energy.begin(), f_energy.end(), 0.0);
                        #endif // !__KE_PE_PRESENT__
                    #endif // __CALC_ENTHROPY__

                    #ifdef __MDN_PROFILING__
                        per_step_timer.cutoff();
                        if (global_timer.count() > 1000*60)
                            logger.info("Profiling: {} ms per step, {} ms per storage", per_step_timer.current(), per_storage_timer.current());
                        logger.info("Profiling: {} ms per step, {} ms per storage based on {} steps", per_step_timer.current(), per_storage_timer.current(), per_step_timer.counter);
                    #endif // !__MDN_PROFILING__
                    logger.flush();
                }

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

            logger.info("End storage: {}", storage.string());
        }

        writer.Close();
        logger.debug("Closed writer");

        #ifdef __MDN_PROFILING__
            per_step_timer.check();
            per_storage_timer.check();
            global_timer.check();
            global_timer.cutoff();
            logger.info("Final profiling results: {} ms per step, {} ms per storage", per_step_timer.current(), per_storage_timer.current());
            logger.info("Profiling: {} ms total wall time", global_timer.current());
        #endif // !__MDN_PROFILING__

        logger.info("End of processing");
        return RETURN_CODES::OK;
    }
} // namespace

#endif // !__MDN_LOGIC__
