#ifndef __MDN_SSTAGE__
#define __MDN_SSTAGE__

#include <vector>
#include <numeric>
#include <tuple>
#include <numbers>
#include <filesystem>
#include <fstream>
#include <span>

#include "mpi.h"
#include "adios2.h"

#include "MDN.hpp"
#include "utils.hpp"
#include "config.hpp"
#include "constants.hpp"
#include "csv.hpp"

#define __MDN_PROFILING__

namespace mdn{

    namespace props{

        inline constexpr double nl(const double& T) noexcept(true) {
            return 0.9367410354674542 * std::exp(-0.46391125909214476 * std::pow(T, 2.791206046910478));
        }
        inline constexpr double sigma(const double& T) noexcept(true) {
            return -1.8111682291065432 * T + 1.8524737887189553;
        }
        inline constexpr double nvs(const double& T) noexcept(true) {
            return 5.85002406256173*std::exp(-(-1.9395800010433437*T + 6.089581679273376)/T);
        }
        inline constexpr double n1s(const double& T) noexcept(true) {
            return 2.05457139*std::exp(-(-1.10143417*T + 4.91381273)/T);
        }
        inline constexpr double nvs_reverse(const double& nv) noexcept(true) {
            return -6.08958/std::log(0.02457499600439634*nv);
        }
    }

    double mean_size(const std::vector<uint64_t>& sizes, const std::span<const uint64_t>& counts) noexcept(true) {
        uint64_t c_sum = 0, p_sum = 0;
        for (uint64_t i = 0; i < sizes.size(); ++i) {
            c_sum += counts[i];
            p_sum += sizes[i]*counts[i];
        }
        return static_cast<double>(p_sum) / c_sum;
    }

    constexpr double condensation_degree(const std::vector<uint64_t>& sizes, const std::span<const uint64_t>& counts, const uint64_t& N, const uint64_t& km) noexcept(true) {
        uint64_t p_sum = 0;
        for (uint64_t i = 0; i < km; ++i) p_sum += sizes[i]*counts[i];
        return 1 - static_cast<double>(p_sum) / N;
    }

    uint64_t max_size(std::span<uint64_t>& counts) noexcept(true) {
        uint64_t i = counts.size() - 1;
        while (counts[i] == 0) --i;
        return i;
    }

    constexpr double nv(const std::vector<uint64_t>& sizes, const std::span<const uint64_t>& counts, const double& volume, const uint64_t& km) noexcept(true) {
        uint64_t p_sum = 0;
        for (uint64_t i = 0; i < km; ++i) p_sum += sizes[i]*counts[i];
        return static_cast<double>(p_sum) / volume;
    }

    constexpr double nd(const std::span<const uint64_t>& counts, const double& volume, const uint64_t& km) noexcept(true) {
        uint64_t p_sum = 0;
        for (uint64_t i = 0; i < km; ++i) p_sum += counts[i];
        return static_cast<double>(p_sum) / volume;
    }

    double nvs(const std::vector<uint64_t>& sizes, const std::span<const uint64_t>& counts, const double& volume, const uint64_t& km, const double& T) noexcept(true) {
        double rl = rl = std::pow(3 / (4 * std::numbers::pi * props::nl(T)), 1.0/3);
        double eco = 2 * props::sigma(T) / (props::nl(T) * T * rl);
        double num = 0, denum = 0, buf = 0;
        for (uint64_t i = km; i < counts.size(); ++i) {
            buf  = std::pow(sizes[i], 2.0/3)*counts[i];
            num += buf;
            denum += buf * std::exp(eco / std::pow(sizes[i], 1.0/3));
        }
        num = num*counts[1] / volume / volume;
        return num / denum;
    }

    constexpr double Srh_props(const double& nv, const double& T) noexcept(true) {return nv / props::nvs(T);}

    constexpr double Srh_nv(const double& _nv, const double& _nvs) noexcept(true) {
        double _Srh = _nv/_nvs;

        if (std::isnan(_Srh) || std::isinf(_Srh)) return 0;
        else return _Srh;
    }

    double S1(const std::span<const uint64_t>& counts, const double& T, const double& volume) noexcept(true) {
        double _S1 = counts[1] / volume / props::n1s(T);

        if (std::isnan(_S1) || std::isinf(_S1)) return 0;
        else return _S1;
    }

    //         0-step    1-r.time  2-x     3-nv    4-T     5-Srh   6-Srh_p 7-nd    8-nd    9-S
    std::tuple<uint64_t, uint64_t, double, double, double, double, double, double, double, double>
    MDN::get_row(const uint64_t& step, const std::vector<uint64_t>& sizes, const std::span<const uint64_t>& counts, const double& volume, const uint64_t& km, const double& T, const uint64_t& N_atoms){
        double _nv = nv(sizes, counts, volume, km);
        double _nvs = nvs(sizes, counts, volume, km, T);

        return std::make_tuple(
            step,
            static_cast<uint64_t>(std::round(time_step*step*distance)),
            condensation_degree(sizes, counts, N_atoms, km),
            _nv,
            T,
            _nvs,
            Srh_nv(_nv, _nvs),
            Srh_props(_nv, T),
            nd(counts, volume, km),
            S1(counts, T, volume)
        );
    }

    void MDN::sstage(){
        adios2::ADIOS adios = adios2::ADIOS(MPI_COMM_SELF);
        logger.debug("ADIOS2 initialized");

        adios2::IO dataio   = adios.DeclareIO("DataReader");
        dataio.SetEngine("BP4");
        adios2::Engine reader = dataio.Open(args.outfile, adios2::Mode::Read, MPI_COMM_SELF);
        logger.debug("ADIOS2 IO initialized");

        AK::CSV_t<std::ofstream&, 10, uint64_t, uint64_t, double, double, double, double, double, double, double, double>
                                      data_csv(args.data_csv, std::ios::out);
        AK::CSV<std::ofstream&, true> dist_csv(args.dist_csv, std::ios::out);
        AK::CSV<std::ofstream&, true> temp_csv(args.temp_csv, std::ios::out);
        uint64_t to_cut = 0;
        if (args.cut > 0){
            dist_csv.fix_columns(args.cut + 1);
            temp_csv.fix_columns(args.cut + 2);
            to_cut = args.cut;
        } else {
            dist_csv.fix_columns(max_cluster_size + 1);
            temp_csv.fix_columns(max_cluster_size + 2);
            to_cut = max_cluster_size;
        }

        adios2::Variable<uint64_t> WvarNstep;
        adios2::Variable<uint64_t> WvarN;
        adios2::Variable<uint64_t> WvarDist;
        adios2::Variable<double>   WvarTemps;
        adios2::Variable<double>   WvarVol;
        adios2::Variable<double>   WvarTTemp;
        #ifdef __CALC_ENTHROPY__
        adios2::Variable<double>   WvarEnthr;
        adios2::Variable<double>   WvarEnerg;
        #endif // __CALC_ENTHROPY__

        uint64_t ntimestep = 0, natoms = 0;
        double volume = 0, total_temp = 0;

        std::unique_ptr<uint64_t[]> dist_buff(new uint64_t[_Natoms + 1]);
        std::unique_ptr<double[]>   temp_buff(new double  [_Natoms + 1]);

        const std::span<const uint64_t> dist_view(dist_buff.get(), _Natoms + 1);
        const std::span<const double> temp_view(temp_buff.get(), _Natoms + 1);

        std::vector<uint64_t> sizes(_Natoms + 1);
        for (size_t i = 0; i < _Natoms + 1; i++) sizes[i] = i;

        timer Tmisc;
        Tmisc.b();
        #ifdef __MDN_PROFILING__
            timer Tglobal;
            timer TPstep;

            timer TADIOS_get_data;
            timer TCSVWrite;
            timer TCSVandGETROW;

            Tglobal.start();
        #endif // __MDN_PROFILING__

        uint64_t currentStep = 0;
        logger.debug("Starting main loop");
        try{
            while (currentStep != total_steps) {
                #ifdef __MDN_PROFILING__
                    TPstep.start();
                #endif // __MDN_PROFILING__

                if (reader.BeginStep() == adios2::StepStatus::EndOfStream){
                    #ifdef __MDN_TRACE_OUT__
                        logger.trace("EOS reached");
                    #endif // __MDN_TRACE_OUT__
                    break;
                }

                currentStep = reader.CurrentStep();

                WvarNstep = dataio.InquireVariable<uint64_t>(std::string("ntimestep" ));
                WvarN     = dataio.InquireVariable<uint64_t>(std::string("natoms"    ));
                WvarDist  = dataio.InquireVariable<uint64_t>(std::string("dist"      ));
                WvarTemps = dataio.InquireVariable<double>  (std::string("temps"     ));
                WvarVol   = dataio.InquireVariable<double>  (std::string("volume"    ));
                WvarTTemp = dataio.InquireVariable<double>  (std::string("total_temp"));
                #ifdef __CALC_ENTHROPY__
                WvarEnthr = dataio.InquireVariable<double>  (std::string("enthr"));
                WvarEnerg = dataio.InquireVariable<double>  (std::string("energ"));
                #endif // __CALC_ENTHROPY__

                if (!(WvarNstep && WvarN && WvarDist && WvarTemps && WvarVol && WvarTTemp)){
                    reader.EndStep();
                    logger.error("Error on read variables at step {}, continuing to next step", currentStep);
                    continue;
                }

                #ifdef __MDN_PROFILING__
                    TADIOS_get_data.b();
                #endif // __MDN_PROFILING__
                reader.Get(WvarNstep, &ntimestep);
                reader.Get(WvarN, &natoms);
                reader.Get(WvarVol, &volume);
                reader.Get(WvarTTemp, &total_temp);
                reader.Get(WvarDist, dist_buff.get());
                reader.Get(WvarTemps, temp_buff.get());
                reader.PerformGets();
                reader.EndStep();
                #ifdef __MDN_PROFILING__
                    TADIOS_get_data.e();
                #endif // __MDN_PROFILING__

                #ifdef __MDN_PROFILING__
                    TCSVWrite.b();
                #endif // __MDN_PROFILING__
                dist_csv << ntimestep;
                dist_csv.write(dist_view);
                temp_csv << ntimestep;
                temp_csv << total_temp;
                temp_csv.write(temp_view);
                #ifdef __MDN_PROFILING__
                    TCSVWrite.e();
                #endif // __MDN_PROFILING__

                #ifdef __MDN_PROFILING__
                    TCSVandGETROW.b();
                #endif // __MDN_PROFILING__
                data_csv.write(get_row(ntimestep, sizes, dist_view, volume, 9, total_temp, natoms));
                #ifdef __MDN_PROFILING__
                    TCSVandGETROW.e();
                #endif // __MDN_PROFILING__

                #ifdef __MDN_PROFILING__
                    TPstep.cutoff();
                #endif // __MDN_PROFILING__

                Tmisc.check();
                if (Tmisc.count() > 60*1000){
                    Tmisc.b();
                    logger.debug("Progress:     {}%", (static_cast<long double>(currentStep) / total_steps)*100.0);

                    #ifdef __MDN_PROFILING__
                        logger.info("###################");
                        logger.info("Profiling:");
                        logger.info("└─Per step:    {} ms", TPstep.current());
                        logger.info("  ├─Input:     {} ms", TADIOS_get_data.sum());
                        logger.info("  ├─CSVWrite:  {} ms", TCSVWrite.sum());
                        logger.info("  ├─CSVandROW: {} ms", TCSVandGETROW.sum());
                        long double auxtime = TPstep.current() - TADIOS_get_data.sum() - TCSVWrite.sum() - TCSVandGETROW.sum();
                        logger.info("  └─Not stated   {} ms", auxtime);
                        logger.info("###################");
                    #endif // !__MDN_PROFILING__

                }
            }
        }catch (std::logic_error& e){
            logger.error("Some std::logic_error happened");
            logger.error(e.what());
            logger.error("Failing");
        }catch(std::exception& e){
            logger.error("Some std::exception happened");
            logger.error(e.what());
            logger.error("Failing");
        }catch (...){
            logger.error("Some error happened");
            logger.error("Error is not inherits std::exception, so exiting... Probably the whole program may be aborted");
            throw;
        }

    }

} // namespace

#endif // !__MDN_SSTAGE__
