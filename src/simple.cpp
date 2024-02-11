#ifndef __MDN_SIMPLE__
#define __MDN_SIMPLE__

#include "adios2.h"
#include "nlohmann/json.hpp"

#include <yas/serialize.hpp>
#include <yas/object.hpp>
#include <yas/std_types.hpp>
#include <yas/mem_streams.hpp>
#include <yas/binary_iarchive.hpp>
#include <yas/binary_oarchive.hpp>

#include "mpi.h"

#include "constants.cpp"
#include "utils.cpp"
#include "barriers.cpp"

#include <fstream>

namespace mdn {
    namespace fs = std::filesystem;

    std::map<int, std::map<fs::path, std::pair<int, int>>> MDN_root::make_distribution(const std::map<fs::path, int> &storages){
        int total_step_count = 0;
        for (const auto &[key, val] : storages)
            total_step_count += val;
        logger.info("Distributing storages: {} steps over {} workers", total_step_count, size);

        int steps_per_worker = std::ceil(static_cast<double>(total_step_count) / size);
        logger.debug("Steps per worker: {}", steps_per_worker);
        int remaining_steps = total_step_count % size;
        logger.debug("Remaining steps: {}", remaining_steps);

        std::map<fs::path, int>::const_iterator it = storages.begin();
        int worker_remain = 0;
        int storage_count = it->second;
        int storage_begin = 0;
        bool advance = true;
        std::map<int, std::map<fs::path, std::pair<int, int>>> distribution;
        for (size_t i = 0; i <= size; ++i){
            logger.trace("Processing {} worker", i);
            if (advance){
                logger.trace("Advancing, {} steps remaining for this worker", steps_per_worker);
                worker_remain = steps_per_worker;
            }else{
                logger.trace("Not advancing, current worker: {}, {} steps remaining", i - 1, worker_remain);
                --i;
                advance = true;
            }
            logger.trace("Storage: {}, {} steps remaining starting from {}", it->first.string(), storage_count, storage_begin);
            if (storage_count < worker_remain){
                logger.trace("Step count in storage less than worker remaining steps");
                logger.trace("Assigning {}, [{}, {}]] to worker", it->first.string(), storage_begin, storage_count);
                distribution[i].insert(std::make_pair(it->first, std::make_pair(storage_begin, storage_count)));
                worker_remain -= storage_count;
                logger.trace("Worker remaining: {}", worker_remain);
                if (++it == storages.end())
                {
                    break;
                }
                storage_count = it->second;
                storage_begin = 0;
                advance = false;
                logger.trace("New storage: {}, [0, {}]", it->first.string(), storage_count);
                logger.trace("Not advancing to next worker");
            }else if (storage_count > worker_remain){
                logger.trace("Remaining steps in storage greater than worker remaining steps, splitting storage");
                logger.trace("Assigning {} [{}, {}]", it->first.string(), storage_begin, storage_count);
                distribution[i].insert(std::make_pair(it->first, std::make_pair(storage_begin, worker_remain)));
                storage_begin += worker_remain;
                storage_count -= worker_remain;
                logger.trace("Storage {} [{}, {}], advancing to the next worker", it->first.string(), storage_begin, storage_count);
            }else if (storage_count == worker_remain){
                logger.trace("Remaining steps in storage equal to worker remaining steps");
                logger.trace("Assigning {} [{}, {}]", it->first.string(), storage_begin, storage_count);
                distribution[i].insert(std::make_pair(it->first, std::make_pair(storage_begin, worker_remain)));
                logger.trace("Advancing to new worker");
                if (++it == storages.end())
                    break;
                storage_count = it->second;
                storage_begin = 0;
                worker_remain = steps_per_worker;
                logger.trace("New storage: {} [0, {}]", it->first.string(), storage_count);
            }
        }
        return distribution;
    }



    RETURN_CODES MDN_root::distribute(const std::map<fs::path, int> &storages, std::map<fs::path, std::pair<int, int>> &_distrib){

        std::map<int, std::map<std::string, std::pair<int, int>>> _distribution;
        std::map<int, std::map<fs::path, std::pair<int, int>>> distribution;

        bool cache_loaded = false;
        if (args.cache){
            logger.info("Trying to load cached distribution");
            std::error_code ecc;
            if (fs::exists(cwd / folders::cache, ecc)){
                if (!ecc)
                    if(fs::exists(cwd / folders::cache / files::cache, ecc)){
                        if (!ecc){
                            int Nworker{};
                            yas::file_istream isf((cwd / folders::cache / files::cache).string().c_str());
                            yas::load<yas::file|yas::json>(isf, YAS_OBJECT_NVP(
                                "distribution", ("distribution", _distribution)
                                ));
                            yas::load<yas::file|yas::json>(isf, YAS_OBJECT_NVP(
                                "NW", ("NW", Nworker)
                                ));
                            if (Nworker == size){
                                ds2p(distribution, _distribution);
                                cache_loaded = true;
                            }else{
                                logger.info("Cached distribution mismatch number of worker, so new distribution will be cumputed and cached");
                            }
                        }
                    }else{
                        logger.info("Cache file {} does not exists, so cached distribution is not loaded", (cwd / folders::cache / files::cache).string());
                    }
            }else{
                logger.info("Cache directory {} does not exists, so cached distribution is not loaded", (cwd / folders::cache).string());
            }
        }
        if (!cache_loaded)
            distribution = make_distribution(storages);

        dp2s(distribution, _distribution);

        logger.debug("Completed distribution:");
        for (const auto &[worker, stors] : distribution){
            logger.debug("Worker: {}", worker);
            for (const auto &[stor, bounds] : stors){
                logger.debug("Storage: '{}' [{}, {}]", stor.string(), bounds.first, bounds.second);
            }
        }

        yas::mem_ostream os;
        yas::binary_oarchive<yas::mem_ostream> oa(os);

        auto oo = YAS_OBJECT("distribution", ("i", _distribution));
        // auto oon = YAS_OBJECT("NW", ("NWM", size));

        if (args.cache){
            logger.info("Cache enabled");
            if (create_d_if_not(cwd / folders::cache, "Cache") != RETURN_CODES::OK){
                logger.warn("Not caching distribution");
            }else{
                logger.info("Trying to write distribution to cache");
                TRY
                    auto odo = YAS_OBJECT_NVP("distribution", ("distribution", _distribution));
                    auto odon = YAS_OBJECT_NVP("NW", ("NW", size));
                    yas::file_ostream osf((cwd / folders::cache / files::cache).string().c_str());
                    // yas::binary_oarchive<yas::file_ostream> oaf(osf);
                    // oaf.serialize(oo);
                    // oaf.serialize(oon);
                    yas::save<yas::file|yas::json>(osf, odo);
                    yas::save<yas::file|yas::json>(osf, odon);
                    osf.flush();
                CATCH_NOTHROW("Error while writing cache (yas error)")
                logger.info("Distribution for {} workers written to cache", size);
            }
        }

        oa.serialize(oo);
        auto buf = os.get_intrusive_buffer();

        logger.debug("Broadcasting distribution buffer size: {}", buf.size);
        MPI_Bcast(const_cast<std::size_t *>(&buf.size), 1, MPI_UNSIGNED_LONG, cs::mpi_root, wcomm);
        // DISTRIBUTION_BARRIER(MPI_COMM_WORLD);
        logger.debug("Broadcasting distribution");
        MPI_Bcast(const_cast<char *>(buf.data), buf.size, MPI_BYTE, cs::mpi_root, wcomm);
        logger.debug("Sent distribution");
        // MPI_Send(buf.data, buf.size, MPI_BYTE, 1, 123, MPI_COMM_WORLD);
        // DISTRIBUTION_BARRIER(MPI_COMM_WORLD);
        // distribution.at(rank);
        for (const auto &[k,v]: distribution.at(rank)){
            _distrib.emplace(k ,v);
        }
        #ifdef __MDN_BUILD_DEBUG_CODE__
        logger.debug("Root distribution emplaced by ptr");
        #endif // __MDN_BUILD_DEBUG_CODE__

        return RETURN_CODES::OK;
    }
}  // namespace mdn

#endif // !__MDN_SIMPLE__
