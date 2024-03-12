#ifndef __MDN_ROOT__
#define __MDN_ROOT__

#include "nlohmann/json.hpp"
#include "adios2.h"
// #include "csv.hpp"

#include "mpi.h"

#include "constants.cpp"
#include "utils.cpp"
#include "simple.cpp"
#include "barriers.cpp"
#include "logic.cpp"
#include "MDN.hpp"

#include <fstream>
#include <memory>
#include <functional>

namespace mdn {

    using json = nlohmann::json;
    namespace fs = std::filesystem;

    const int MDN_root::get_total_steps(adios2::IO &io, const fs::path &storage){
        adios2::Engine reader = io.Open(storage, adios2::Mode::Read);
        int steps = reader.Steps();
        reader.Close();
        return steps;
    }

    const bool MDN_root::save_storages(const std::map<fs::path, int>& _storages, const unsigned long hash){
        std::map<std::string, int> storages;
        for (const auto& [k, v]: _storages){
            storages.emplace(k.string(), v);
        }
        logger.info("Trying to cache storages");
        fs::path cache_folder = cwd / folders::cache;
        fs::path cache_file = cwd / folders::cache / files::storages_cache;

        if (create_d_if_not(cache_folder, "Cache") != RETURN_CODES::OK){
            logger.warn("Not caching storages");
            return false;
        }
        std::error_code ecc;
        if (!fs::remove(cache_file, ecc)){
            if (ecc){
                logger.error("Cannot remove previous cache file: {}, so cache cannot be saved", cache_file.string());
                return false;
            }
        }
        logger.info("Trying to write storages to cache");
        TRY
            auto odo = YAS_OBJECT_NVP("storages", ("storages", storages));
            auto odon = YAS_OBJECT_NVP("hash", ("hash", hash));
            yas::file_ostream osf(cache_file.string().c_str());
            // yas::binary_oarchive<yas::file_ostream> oaf(osf);
            // oaf.serialize(oo);
            // oaf.serialize(oon);
            yas::save<yas::file|yas::json>(osf, odo);
            yas::save<yas::file|yas::json>(osf, odon);
            osf.flush();
        CATCH_NOTHROW_RET("Error while writing cache (yas error)", false)
        logger.info("Storages cache written to cache file: {}", cache_file.string());
        return true;
    }

    const bool MDN_root::load_storages(std::map<fs::path, int>& _storages, const unsigned long hash){
        std::map<std::string, int> storages;
        logger.info("Trying to load cached storages");
        fs::path cache_folder = cwd / folders::cache;
        fs::path cache_file = cwd / folders::cache / files::storages_cache;
        std::error_code ecc;
        if (fs::exists(cache_folder, ecc)){
                if(fs::exists(cache_file, ecc)){
                    unsigned long _hash{};
                    yas::file_istream isf(cache_file.string().c_str());
                    yas::load<yas::file|yas::json>(isf, YAS_OBJECT_NVP(
                        "storages", ("storages", storages)
                        ));
                    yas::load<yas::file|yas::json>(isf, YAS_OBJECT_NVP(
                        "hash", ("hash", _hash)
                        ));
                    if (_hash == hash){
                        for (const auto& [k, v]: storages){
                            _storages.emplace(fs::path(k), v);
                        }
                        return true;
                    }else{
                        logger.warn("Cached storages mismatch hash, so storages will be computed and cached");
                        return false;
                    }
                }else{
                    if (!ecc){
                        logger.warn("Cache file {} does not exists, so cache is not loaded", cache_file.string());
                    }else{
                        logger.warn("Cannot check cache file {} existense due to error {}, so cache is not loaded", cache_file.string(), ecc.message());
                    }
                    return false;
                }
        }else{
            if (!ecc){
                logger.warn("Cache directory {} does not exists, so cache is not loaded", cache_folder.string());
            }else{
                logger.warn("Cannot check cache directory {} existense due to error {}, so cache is not loaded", cache_folder.string(), ecc.message());
            }
            return false;
        }
        return false;
    }

    const std::map<fs::path, int> MDN_root::storage_rsolve(const json &_storages){
        adios2::ADIOS adios = adios2::ADIOS(MPI_COMM_SELF);
        adios2::IO io = adios.DeclareIO("Steps_reader");
        logger.debug("ADIOS2.IO initialized");
        std::map<fs::path, int> storages;
        int steps = 0;
        fs::path rsolved_storage;
        for (const auto &storage : _storages.items()){
            TRY
            rsolved_storage = fs::absolute(storage.key());
            logger.trace("Storage '{}', resolved path: {}", storage.key(), rsolved_storage.string());
            steps = get_total_steps(io, rsolved_storage);
            storages.emplace(rsolved_storage, steps);
            logger.trace("Storage '{}' contains {} steps", rsolved_storage.string(), steps);
            CATCH("Error while resolving storage: " + storage.key())
        }
        return storages;
    }

    const uint64_t MDN_root::bearbeit(const fs::path &storage, double *dims){
        adios2::ADIOS adios = adios2::ADIOS();
        adios2::IO io = adios.DeclareIO("PReader");
        adios2::Engine reader;
        TRY
            reader = io.Open(storage, adios2::Mode::Read);
        CATCH("Error while opening storage" + storage.string())

        reader.BeginStep();

        adios2::Variable<uint64_t> varNatoms;
        adios2::Variable<double>   varboxxhi;
        adios2::Variable<double>   varboxyhi;
        adios2::Variable<double>   varboxzhi;

        TRY
            varNatoms = io.InquireVariable<uint64_t>(std::string(lcf::natoms));
            varboxxhi = io.InquireVariable<double>(std::string(lcf::boxxhi));
            varboxyhi = io.InquireVariable<double>(std::string(lcf::boxyhi));
            varboxzhi = io.InquireVariable<double>(std::string(lcf::boxzhi));
        CATCH("Error while inquiring variables")

        uint64_t Natoms = 0;

        TRY
            reader.Get(varNatoms, Natoms);
            reader.Get(varboxxhi, dims + 0);
            reader.Get(varboxyhi, dims + 1);
            reader.Get(varboxzhi, dims + 2);
        CATCH("Error while getting variables")

        reader.EndStep();
        reader.Close();
        return Natoms;
    }

    uint64_t MDN_root::dns(std::map<fs::path, std::pair<int, int>> &_distrib){
        logger.debug("Reading datafile");
        fs::path datafile_path = cwd / files::data;
        if (!fs::exists(datafile_path)){
            logger.error("Datafile (file: '" + datafile_path.string() + "') does not exists");
            throw std::runtime_error("Datafile (file: '" + datafile_path.string() + "') does not exists");
        }
        std::ifstream datafile_in(datafile_path);
        json data;
        TRY
            data = json::parse(datafile_in);
        CATCH("Error while parsing json from datafile")
        logger.debug("Successfully read datafile");

        logger.debug("Resolving storages");
        std::map<fs::path, int> storages;

        const unsigned long shash = std::hash<json>{}(data.at(fields::storages));
        if (!load_storages(storages, shash)){
            TRY
                storages = storage_rsolve(data.at(fields::storages));
            CATCH("Error while resolving storages")
        }
        save_storages(storages, shash);

        fs::path sto_check = storages.begin()->first;
        logger.debug("Simulation parameters will be obtained from storage: {}", sto_check.string());

        std::unique_ptr<double[]> dims(new double[3]);
        uint64_t Natoms{};
        TRY
            Natoms = bearbeit(sto_check, dims.get());
        CATCH("Error while getting simulation parameters from one of storages")

        const std::array<double, 3> _dims = {*dims.get(), *(dims.get() + 1), *(dims.get() + 2)};
        const double Volume = dims[0] * dims[1] * dims[2];
        logger.debug("Number of atoms: {}", Natoms);
        logger.debug("Dimensions: [{}, {}, {}]", dims[0], dims[1], dims[2]);
        logger.debug("System volume: {}", Volume);

        TRY
            data[fields::Natoms] = Natoms;
            data[fields::Dims] = _dims;
            data[fields::Volume] = Volume;
            logger.debug("Updating datafile");
            std::ofstream datafile_out(datafile_path);
            datafile_out << std::setw(4) << data << std::endl;
        CATCH("Error while writing datafile")

        fs::path data_processing_folder = cwd / folders::post_process_folder;
        logger.debug("Creating post-processing directory: " + data_processing_folder.string());

        std::error_code ec;
        if (!fs::exists(data_processing_folder, ec)){
            if (!ec){
                if (!(fs::create_directory(data_processing_folder, ec) && ec)){
                    logger.error("Can not create directory '{}' due to error with code: {}. Message", data_processing_folder.string(), ec.value());
                    logger.error(ec.message());
                    throw std::runtime_error("Can not create post-processing directory: " + data_processing_folder.string());
                }
            }else{
                logger.error("Can not check directory '{}' existense due to error with code: {}. Message", data_processing_folder.string(), ec.value());
                logger.error(ec.message());
                throw std::runtime_error("Can not check post-processing directory existense: " + data_processing_folder.string());
            }
        }

        logger.debug("Distributing storages");
        TRY
            distribute(storages, _distrib);
        CATCH("Error while distributing storages")

        logger.debug("Bcasting Natoms");
        MPI_Bcast(&Natoms, 1, MPI_UINT64_T, cs::mpi_root, wcomm);
        logger.debug("Sent Natoms");
        logger.info("Distribution complete");

        return Natoms;
    }

    RETURN_CODES MDN_root::sanity(){
        int rnd = std::rand();
        std::cout << "      Bcasting: " << rnd << std::endl;
        MPI_Bcast(&rnd, 1, MPI_INT, cs::mpi_root, wcomm);
        // int responce = 0;
        int *responces = new int[size];
        std::cout << "      Gathering..." << std::endl;
        MPI_Gather(&rnd, 1, MPI_INT, responces, 1, MPI_INT, cs::mpi_root, wcomm);

        if (responces[0] != rnd){
            std::cerr << "Root sanity doesn't passed" << std::endl;
            throw std::runtime_error("Root sanity doesn't passed");
        }

        // logger.info(std::string("Random int: " + std::to_string(rnd)));
        // logger.info(std::string("Received arrray:"));
        for (int i = 0; i < size; ++i)
        {
            // logger.info(std::to_string(responces[i]));
            if (responces[i] - i != rnd)
            {
                std::cerr << "Root sanity doesn't passed" << std::endl;
                throw std::runtime_error("Root sanity doesn't passed");
            }
        }
        std::cout << "      Root sanity OK" << std::endl;
        std::cout << "      Releasing sanity barrier" << std::endl;
        SANITY_BARRIER(wcomm);
        std::cout << "      Sanity barrier released" << std::endl;

        delete[] responces;
        return RETURN_CODES::OK;
    }

    void MDN_root::gather_storages(std::map<int, std::string>& storages){
        std::set<int> ncw;
        for (int i = 1; i < size; ++i) ncw.emplace(i);
        MPI_Status status;
        int flag = 0;
        int count{};
        while (!ncw.empty()){
            for (const int i : ncw){
                MPI_Iprobe(i, MPI_TAGS::STOR, wcomm, &flag, &status);
                if (flag){
                    MPI_Get_count(&status, MPI_CHAR, &count);
                    char* buf = new char[count];
                    MPI_Recv(buf, count, MPI_CHAR, i, MPI_TAGS::STOR, wcomm, &status);
                    storages.emplace(i, std::string(buf, count));
                    delete buf;
                    flag = 0;
                    status = MPI_Status();
                    ncw.erase(i);
                }
            }
        }
    }

    void MDN_root::gen_matrix(std::map<int, std::string>& storages, const uint64_t Natoms, const u_int64_t max_cluster_size){
        adios2::ADIOS adios = adios2::ADIOS(MPI_COMM_SELF);
        adios2::IO dataio   = adios.DeclareIO("DATAREADER");
        dataio.SetEngine("BP4");

        adios2::Variable<uint64_t> WvarNstep;
        adios2::Variable<uint64_t> WvarN;
        adios2::Variable<uint64_t> WvarDist;
        adios2::Variable<double>   WvarTemps;
        adios2::Variable<double>   WvarVol;
        adios2::Variable<double>   WvarTTemp;

        std::unique_ptr<uint64_t[]> _sizes_counts(new uint64_t[Natoms + 1]);
        std::unique_ptr<double[]> _temps_by_size(new double[Natoms + 1]);
        std::span<const uint64_t> sizes_counts(_sizes_counts.get(), max_cluster_size);
        std::span<const double> temps_by_size(_temps_by_size.get(), max_cluster_size);

        uint64_t timestep{}, _Natoms{};
        double Volume{}, total_temp{};
        std::cout << std::endl;
        // std::vector<uint64_t> sizes_counts(Natoms + 1, 0UL);
        // std::vector<double> temps_by_size(Natoms + 1, 0.0);

        // std::ofstream matrix_csv_file(args.outfile.parent_path() / files::matrix_csv, std::ios_base::out);
        // std::ofstream temps_csv_file(args.outfile.parent_path() / files::temps_csv, std::ios_base::out);
        // auto matrix_writer = csv::make_csv_writer(matrix_csv_file);
        // auto temps_writer = csv::make_csv_writer(temps_csv_file);
        csvWriter matrix_writer(args.outfile.parent_path() / files::matrix_csv);
        csvWriter temps_writer(args.outfile.parent_path() / files::temps_csv);

        for (const auto&[rank, storage] : storages){
            adios2::Engine reader = dataio.Open(storage, adios2::Mode::Read, MPI_COMM_SELF);
            uint64_t currentStep = 0;
            while ((reader.BeginStep() == adios2::StepStatus::EndOfStream)){
                currentStep = reader.CurrentStep();

                WvarNstep  = dataio.InquireVariable<uint64_t>("ntimestep");
                WvarN      = dataio.InquireVariable<uint64_t>("Natoms");
                WvarDist   = dataio.InquireVariable<uint64_t>("dist");
                WvarTemps  = dataio.InquireVariable<double>("temps");
                WvarVol    = dataio.InquireVariable<double>("volume");
                WvarTTemp  = dataio.InquireVariable<double>("total_temp");

                if (!(WvarNstep && WvarN && WvarDist && WvarTemps && WvarVol && WvarTTemp))
                    continue;

                reader.Put(WvarNstep, timestep);
                reader.Put(WvarDist,   _sizes_counts.get());
                reader.Put(WvarVol,   Volume);
                reader.Put(WvarN,     _Natoms);
                reader.Put(WvarTTemp, total_temp);
                reader.Put(WvarTemps,  _temps_by_size.get());

                matrix_writer << timestep << _Natoms << sizes_counts << "\n";
                temps_writer << timestep << total_temp << temps_by_size << "\n";

                reader.EndStep();
            }
            matrix_writer.flush();
            temps_writer.flush();
        }
    }

    RETURN_CODES MDN_root::setup(){
        std::cout << "Sanity:" << std::endl;
        TRY
            sanity();
        CATCH_NOLOGGER("Error while checking sanity")

        std::error_code ec;
        fs::path log_folder = cwd / folders::log;
        if (!fs::exists(log_folder)){
            if (!fs::create_directories(log_folder, ec)){
                std::cerr << "Cannot create non-existent directories: " << log_folder.string() << std::endl;
                std::cerr << "Error code: " << ec.value() << std::endl;
                std::cerr << "Message: " << ec.message() << std::endl;
                throw std::runtime_error("Cannot create non-existent directory: " + log_folder.string());
            }else{
                std::cout << "Created logging directory: " << log_folder << std::endl;
            }
        }else{
            std::cout << "Logging folder exists" << std::endl;
        }

        std::cout << "  Releasing setup barrier" << std::endl;
        SETUP_BARRIER(wcomm);
        std::cout << "  Setup barrier released" << std::endl;

        TRY
            setup_logger(rank);
        CATCH("Error while setting up logger")

        std::cout << "  Root setup OK" << std::endl;
        return RETURN_CODES::OK;
    }

    RETURN_CODES MDN_root::entry(){
        std::cout << "Setup..." << std::endl;
        TRY
            setup();
        CATCH_NOLOGGER("Error wahile setting up")
        logger.info("Initialized");

        std::map<fs::path, std::pair<int, int>> distribution;
        uint64_t summ = 0;
        uint64_t Natoms{};
        TRY
            Natoms = dns(distribution);
        CATCH("Error while getting distribution")
        for (const auto &[k,v]: distribution){
            logger.info("{}: {}({}) - {}", k.string(), v.first, v.second, v.first + v.second);
            summ += v.second;
        }
        logger.info("Total: {}", summ);

        fs::path ss = args.outfile.parent_path() / (args.outfile.filename().string() + "." + std::to_string(rank));
        logger.info("Output data will be written to {}", ss.string());

        uint64_t max_cluster_size = 0;
        logger.info("Starting calculations...");
        TRY
            run(ss, distribution, Natoms, max_cluster_size);
        CATCH("Error while doing caculations")
        logger.info("Calculations ended, gathering result info");

        uint64_t *max_sizes = new uint64_t[size];
        MPI_Gather(&max_cluster_size, 1, MPI_UINT64_T, max_sizes, 1, MPI_UINT64_T, cs::mpi_root, wcomm);
        max_cluster_size = *std::max_element(max_sizes, max_sizes + size);
        delete max_sizes;
        logger.debug("Gathered max cluster size");

        std::map<int, std::string> storages;
        gather_storages(storages);
        storages.emplace(rank, ss.string());
        logger.debug("Gathered storages");

        logger.info("Generating output matrices");
        gen_matrix(storages, Natoms, max_cluster_size);

        logger.info("Exiting entry point. NO RETURN");
        return RETURN_CODES::OK;
    }
} // namespace mdn

#endif // !__MDN_ROOT__
