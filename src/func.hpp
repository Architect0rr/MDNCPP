#ifndef __MDN_COMPUTE__
#define __MDN_COMPUTE__

#include <vector>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <span>

#include "memory.hpp"
#include "enums.hpp"

namespace mdn{

    class ComputeModule{
    public:
        ComputeModule(Memory* memptr, uint64_t _Natoms, int ndim): _Natoms(_Natoms), ndim(ndim), memptr(memptr) {
            ucids.reserve(_Natoms);
            sizes_counts.resize(_Natoms + 1, 0UL);
            temps_by_size.resize(_Natoms + 1, 0.0);
            sqvels.resize(_Natoms, 0.0);
        }

        void c_pid_by_cid(const std::span<const double> &cids){
            pid_by_cid.clear();
            for (uint64_t i = 0; i < Natoms; ++i) pid_by_cid[cids[i]].emplace_back(i);
        }

        void c_ucids(const std::span<const double> &cids){
            ucids.clear();
            for (size_t i = 0; i < Natoms; ++i) ucids.emplace(cids[i]);
        }

        void c_by_size(){
            cluster_ids_by_size.clear();
            particles_by_size.clear();
            uint64_t c_size = 0;
            for (const uint64_t &i : ucids) {
                c_size = pid_by_cid[i].size();
                cluster_ids_by_size[c_size].push_back(i);
                particles_by_size[c_size].insert(particles_by_size[c_size].end(), pid_by_cid[i].begin(), pid_by_cid[i].end());
            }
        }

        void c_sqvels(){
            std::fill(sqvels.begin(), sqvels.end(), 0);
            for (size_t i = 0; i < _Natoms; ++i) {
                sqvels[i] += std::pow(velX[i], 2) + std::pow(velY[i], 2) + std::pow(velZ[i], 2);
            }
        }

        void c_ke(){
            std::fill(kes.begin(), kes.end(), 0);
            for (size_t i = 0; i < _Natoms; ++i) kes[i] = masses[i]*sqvels[i] / 2;
        }

        void c_total_temp(){
            total_temp = 0;
            for (const double& ke: kes) total_temp += ke;
            total_temp = 2 * total_temp / ((Natoms - 1) * ndim);
        }

        inline uint64_t Nclusters(){
            return ucids.size();
        }

    private:
        int ndim = 3;
        uint64_t _Natoms = 0;
        uint64_t timestep = 0, Natoms = 0, Nclusters = 0;
        double boxxhi = 0, boxyhi = 0, boxzhi = 0, boxxlo = 0, boxylo = 0, boxzlo = 0, Volume = 0;
        double total_temp = 0, ke = 0, pe = 0;
        std::unordered_set<uint64_t> ucids;
        std::map<uint64_t, std::vector<uint64_t>> pid_by_cid, cluster_ids_by_size, particles_by_size;
        std::vector<uint64_t> sizes_counts;
        std::vector<double> temps_by_size;
        std::vector<double> sqvels;
        // std::map<uint64_t, double> particle_count_by_size, ndofs_by_size;

        Memory * memptr;
    };

} // namespace

#endif // !__MDN_COMPUTE__
