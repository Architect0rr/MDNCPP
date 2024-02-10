#ifndef __MDN_SIMPLE__
#define __MDN_SIMPLE__

#include "adios2.h"
#include "nlohmann/json.hpp"
#include <yas/serialize.hpp>
#include <yas/std_types.hpp>
#include <yas/mem_streams.hpp>
#include <yas/binary_iarchive.hpp>
#include <yas/binary_oarchive.hpp>

#include "mpi.h"

#include "constants.cpp"
#include "utils.cpp"
#include "barriers.cpp"

#include <fstream>

// namespace mdn {
namespace fs = std::filesystem;

std::map<int, std::map<fs::path, std::pair<int, int>>> MDN_root::make_distribution(const std::map<fs::path, int> &storages)
{
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
    for (size_t i = 0; i <= size; ++i)
    {
        logger.trace("Processing {} worker", i);
        if (advance)
        {
            logger.trace("Advancing, {} steps remaining for this worker", steps_per_worker);
            worker_remain = steps_per_worker;
        }
        else
        {
            logger.trace("Not advancing, current worker: {}, {} steps remaining", i - 1, worker_remain);
            --i;
            advance = true;
        }
        logger.trace("Storage: {}, {} steps remaining starting from {}", it->first.string(), storage_count, storage_begin);
        if (storage_count < worker_remain)
        {
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
        }
        else if (storage_count > worker_remain)
        {
            logger.trace("Remaining steps in storage greater than worker remaining steps, splitting storage");
            logger.trace("Assigning {} [{}, {}]", it->first.string(), storage_begin, storage_count);
            distribution[i].insert(std::make_pair(it->first, std::make_pair(storage_begin, worker_remain)));
            storage_begin += worker_remain;
            storage_count -= worker_remain;
            logger.trace("Storage {} [{}, {}], advancing to the next worker", it->first.string(), storage_begin, storage_count);
        }
        else if (storage_count == worker_remain)
        {
            logger.trace("Remaining steps in storage equal to worker remaining steps");
            logger.trace("Assigning {} [{}, {}]", it->first.string(), storage_begin, storage_count);
            distribution[i].insert(std::make_pair(it->first, std::make_pair(storage_begin, worker_remain)));
            logger.trace("Advancing to new worker");
            if (++it == storages.end())
            {
                break;
            }
            storage_count = it->second;
            storage_begin = 0;
            worker_remain = steps_per_worker;
            logger.trace("New storage: {} [0, {}]", it->first.string(), storage_count);
        }
    }
    return distribution;
}

RETURN_CODES MDN_root::distribute(const std::map<fs::path, int> &storages, std::map<fs::path, std::pair<int, int>> *_distrib)
{
    std::map<int, std::map<fs::path, std::pair<int, int>>> distribution = make_distribution(storages);
    std::map<int, std::map<std::string, std::pair<int, int>>> _distribution;

    logger.debug("Completed distribution:");
    for (const auto &[worker, stors] : distribution)
    {
        logger.debug("Worker: {}", worker);
        for (const auto &[stor, bounds] : stors)
        {
            _distribution[worker].insert(std::make_pair(stor.string(), bounds));
            logger.debug("Storage: '{}' [{}, {}]", stor.string(), bounds.first, bounds.second);
        }
    }

    yas::mem_ostream os;
    yas::binary_oarchive<yas::mem_ostream> oa(os);

    auto oo = YAS_OBJECT("distribution", _distribution);
    oa.serialize(oo);
    auto buf = os.get_intrusive_buffer();

    logger.debug("Broadcasting distribution buffer size: {}", buf.size);
    MPI_Bcast(const_cast<std::size_t *>(&buf.size), 1, MPI_UNSIGNED_LONG, cs::mpi_root, wcomm);
    // DISTRIBUTION_BARRIER(MPI_COMM_WORLD);
    logger.debug("Broadcasting distribution");
    MPI_Bcast(const_cast<char *>(buf.data), buf.size, MPI_BYTE, cs::mpi_root, wcomm);
    // MPI_Send(buf.data, buf.size, MPI_BYTE, 1, 123, MPI_COMM_WORLD);
    // DISTRIBUTION_BARRIER(MPI_COMM_WORLD);
    *_distrib = distribution.at(rank);
    return RETURN_CODES::OK;
}
// }  // namespace mdn
#endif // !__MDN_SIMPLE__