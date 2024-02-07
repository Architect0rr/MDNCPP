#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include <xtensor/xadapt.hpp>
#include <xtensor/xsort.hpp>
#include "xtensor/xmath.hpp"
#include "xtensor/xrandom.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xbroadcast.hpp"
#include "xtensor/xmasked_view.hpp"

#include <iostream>
#include <vector>
#include <chrono>
#include <ranges>
#include <unordered_set>

typedef unsigned long ul;

struct runner1
{
    int ndim = 3;
    ul Nsamples{}, Nclusters{}, size{};
    xt::xarray<ul> unique_cluster_ids{}, sizes_counts, unique_sizes_row, particle_counts_by_size, ndofs_by_size;
    std::map<ul, xt::xarray<bool>> particle_mask_by_cluster_id, particle_mask_by_size;
    std::map<ul, std::vector<ul>> cluster_ids_by_size;
    xt::xarray<double> kes, buf_by_size, temps_by_size;
    double total_temp{};

    double run(const xt::xarray<double> &velocities, const xt::xarray<double> &masses, const xt::xarray<ul> &particle_ids, const xt::xarray<ul> &cluster_ids)
    {
        Nsamples = particle_ids.shape()[0];

        unique_cluster_ids = xt::unique(cluster_ids);
        Nclusters = unique_cluster_ids.shape()[0];

        size = 0;
        for (ul i : unique_cluster_ids)
        {
            particle_mask_by_cluster_id.emplace(i, xt::equal(cluster_ids, i));
            size = xt::sum<ul>(particle_mask_by_cluster_id[i])();
            cluster_ids_by_size[size].emplace_back(i);
            particle_mask_by_size[size] += particle_mask_by_cluster_id[i];
        }
        sizes_counts = xt::zeros<ul>({Nsamples + 1});
        for (const auto &[k, v] : cluster_ids_by_size)
        {
            sizes_counts(k) = v.size();
        }
        kes = masses * velocities / 2;
        buf_by_size = xt::zeros<double>({Nsamples + 1});
        for (auto &[k, v] : particle_mask_by_size)
        {
            buf_by_size(k) = xt::sum<double>(xt::index_view(kes, xt::argwhere(v)))();
        }
        unique_sizes_row = xt::arange<ul>(0, Nsamples + 1);
        particle_counts_by_size = unique_sizes_row * sizes_counts;
        ndofs_by_size = (particle_counts_by_size - 1) * ndim;
        temps_by_size = xt::pow(buf_by_size / ndofs_by_size, 2);
        total_temp = std::pow(xt::sum<double>(buf_by_size)() / ((Nsamples - 1) * ndim), 2);
        return total_temp;
    }
    void clear()
    {
        particle_mask_by_cluster_id.clear();
        cluster_ids_by_size.clear();
        particle_mask_by_size.clear();
    }
};

template <typename Iter>
using select_access_type_for = std::conditional_t<
    std::is_same_v<Iter, std::vector<bool>::iterator> ||
        std::is_same_v<Iter, std::vector<bool>::const_iterator>,
    typename Iter::value_type,
    typename Iter::reference>;

template <typename Iter1, typename Iter2>
class zip_iterator
{
public:
    using value_type = std::pair<
        select_access_type_for<Iter1>,
        select_access_type_for<Iter2>>;

    zip_iterator() = delete;

    zip_iterator(Iter1 iter_1_begin, Iter2 iter_2_begin)
        : m_iter_1_begin{iter_1_begin}, m_iter_2_begin{iter_2_begin}
    {
    }

    auto operator++() -> zip_iterator &
    {
        ++m_iter_1_begin;
        ++m_iter_2_begin;
        return *this;
    }

    auto operator++(int) -> zip_iterator
    {
        auto tmp = *this;
        ++*this;
        return tmp;
    }

    auto operator!=(zip_iterator const &other)
    {
        return !(*this == other);
    }

    auto operator==(zip_iterator const &other)
    {
        return m_iter_1_begin == other.m_iter_1_begin ||
               m_iter_2_begin == other.m_iter_2_begin;
    }

    auto operator*() -> value_type
    {
        return value_type{*m_iter_1_begin, *m_iter_2_begin};
    }

private:
    Iter1 m_iter_1_begin;
    Iter2 m_iter_2_begin;
};

template <typename T>
using select_iterator_for = std::conditional_t<
    std::is_const_v<std::remove_reference_t<T>>,
    typename std::decay_t<T>::const_iterator,
    typename std::decay_t<T>::iterator>;

template <typename T, typename U>
class zipper
{
public:
    using Iter1 = select_iterator_for<T>;
    using Iter2 = select_iterator_for<U>;

    using zip_type = zip_iterator<Iter1, Iter2>;

    template <typename V, typename W>
    zipper(V &&a, W &&b)
        : m_a{a}, m_b{b}
    {
    }

    auto begin() -> zip_type
    {
        return zip_type{std::begin(m_a), std::begin(m_b)};
    }
    auto end() -> zip_type
    {
        return zip_type{std::end(m_a), std::end(m_b)};
    }

private:
    T m_a;
    U m_b;
};

template <typename T, typename U>
auto zip(T &&t, U &&u)
{
    return zipper<T, U>{std::forward<T>(t), std::forward<U>(u)};
}

struct runcer{
    int ndim = 3;
    ul Nsamples{}, Nclusters{}, size{};
    std::unordered_set<ul> unique_cluster_ids;
    std::map<ul, std::vector<ul>> particles_by_cluster_id, cluster_ids_by_size, particles_by_size;
    std::vector<ul> sizes_counts;
    std::vector<double> kes;
    double total_temp{}, op{};
    std::map<ul, double> temps_by_size;

    runcer(const ul &Nsamples) : Nsamples(Nsamples), sizes_counts(Nsamples + 1, 0UL), kes(Nsamples, 0.0)
    {
        unique_cluster_ids.reserve(Nsamples);
    }

    double run(const std::vector<double> &velocities, const std::vector<double> &masses, const std::vector<ul> &particle_ids, const std::vector<ul> &cluster_ids)
    {
        Nsamples = particle_ids.size();
        for (const ul &i : cluster_ids)
        {
            unique_cluster_ids.insert(i);
        }
        Nclusters = unique_cluster_ids.size();



        size = 0;
        for (const ul &i : unique_cluster_ids)
        {
            for (ul j = 0; j < Nsamples; ++j)
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
        for (const auto &[vel, mass] : zip(velocities, masses))
        {
            kes.emplace_back(mass * vel / 2);
            total_temp += kes.back();
        }
        total_temp = std::pow(total_temp / ((Nsamples - 1) * ndim), 2);
        // std::map<ul, double> kes_by_size, particle_count_by_size, ndofs_by_size;

        op = 0;
        for (const auto &[k, v] : particles_by_size)
        {
            op = 0;
            for (ul i : v)
            {
                op += kes[i];
            }
            // kes_by_size[k] = op;
            // particle_count_by_size = v.size();
            // ndofs_by_size[k] = (v.size() - 1) * ndim;
            temps_by_size.emplace(k, std::pow(op / ((v.size() - 1) * ndim), 2));
        }
        // for (const auto &[k, v]: temps_by_size){
        //     std::cout << k << ":" << v << std::endl;
        // }
        return total_temp;
    }
    void clear(){
        unique_cluster_ids.clear();
        particles_by_cluster_id.clear();
        cluster_ids_by_size.clear();
        particles_by_size.clear();
        sizes_counts.clear();
        kes.clear();
        temps_by_size.clear();
    }
};

double runc(const std::vector<double> &velocities, const std::vector<double> &masses, const std::vector<ul> &particle_ids, const std::vector<ul> &cluster_ids)
{
    int ndim = 3;
    ul Nsamples = particle_ids.size();
    std::unordered_set<ul> unique_cluster_ids;
    unique_cluster_ids.reserve(Nsamples);
    for (const ul &i : cluster_ids)
    {
        unique_cluster_ids.insert(i);
    }
    ul Nclusters = unique_cluster_ids.size();

    std::map<ul, std::vector<ul>> particles_by_cluster_id;
    std::map<ul, std::vector<ul>> cluster_ids_by_size;
    std::map<ul, std::vector<ul>> particles_by_size;

    ul size = 0;
    for (const ul &i : unique_cluster_ids)
    {
        for (ul j = 0; j < Nsamples; ++j)
        {
            if (cluster_ids[j] == i)
                particles_by_cluster_id[i].push_back(j);
        }
        size = particles_by_cluster_id[i].size();
        cluster_ids_by_size[size].push_back(i);
        particles_by_size[size].insert(particles_by_size[size].end(), particles_by_cluster_id[i].begin(), particles_by_cluster_id[i].end());
    }

    std::vector<ul> sizes_counts(Nsamples + 1, 0UL);
    for (const auto &[k, v] : cluster_ids_by_size)
    {
        sizes_counts[k] = v.size();
    }
    std::vector<double> kes;
    kes.reserve(Nsamples);
    double total_temp{};
    for (const auto &[vel, mass] : zip(velocities, masses))
    {
        kes.emplace_back(mass * vel / 2);
        total_temp += kes.back();
    }
    total_temp = std::pow(total_temp / ((Nsamples - 1) * ndim), 2);
    // std::map<ul, double> kes_by_size, particle_count_by_size, ndofs_by_size;
    std::map<ul, double> temps_by_size;
    double op{};
    for (const auto &[k, v] : particles_by_size)
    {
        op = 0;
        for (ul i : v)
        {
            op += kes[i];
        }
        // kes_by_size[k] = op;
        // particle_count_by_size = v.size();
        // ndofs_by_size[k] = (v.size() - 1) * ndim;
        temps_by_size.emplace(k, std::pow(op / ((v.size() - 1) * ndim), 2));
    }
    // for (const auto &[k, v]: temps_by_size){
    //     std::cout << k << ":" << v << std::endl;
    // }
    return total_temp;
}

double run(const xt::xtensor<double, 1> &velocities, const xt::xtensor<double, 1> &masses, const xt::xtensor<ul, 1> &particle_ids, const xt::xtensor<ul, 1> &cluster_ids)
{
    int ndim = 3;
    ul Nsamples = particle_ids.shape()[0];
    xt::xtensor<ul, 1> unique_cluster_ids(xt::unique(cluster_ids));
    ul Nclusters = unique_cluster_ids.shape()[0];

    std::map<ul, xt::xtensor<bool, 1>> particle_mask_by_cluster_id;
    std::map<ul, std::vector<ul>> cluster_ids_by_size;
    std::map<ul, xt::xtensor<bool, 1>> particle_mask_by_size;

    ul size = 0;
    for (const ul i : unique_cluster_ids)
    {
        particle_mask_by_cluster_id.emplace(i, xt::equal(cluster_ids, i));
        size = xt::sum<ul>(particle_mask_by_cluster_id[i])();
        cluster_ids_by_size[size].emplace_back(i);
        if (!particle_mask_by_size.contains(size))
        {
            particle_mask_by_size.emplace(size, particle_mask_by_cluster_id[i]);
        }
        else
        {
            particle_mask_by_size[size] += particle_mask_by_cluster_id[i];
        }
    }
    xt::xtensor<ul, 1> sizes_counts = xt::zeros<ul>({Nsamples + 1});
    for (const auto &[k, v] : cluster_ids_by_size)
    {
        sizes_counts(k) = v.size();
    }
    xt::xtensor<double, 1> kes = masses * velocities / 2;
    xt::xtensor<double, 1> buf_by_size = xt::zeros<double>({Nsamples + 1});
    for (auto &[k, v] : particle_mask_by_size)
    {
        buf_by_size(k) = xt::sum<double>(xt::index_view(kes, xt::argwhere(v)))();
    }
    xt::xtensor<ul, 1> unique_sizes_row = xt::arange<ul>(0, Nsamples + 1);
    xt::xtensor<ul, 1> particle_counts_by_size = unique_sizes_row * sizes_counts;
    xt::xtensor<ul, 1> ndofs_by_size = (particle_counts_by_size - 1) * ndim;
    xt::xtensor<double, 1> temps_by_size = xt::pow(buf_by_size / ndofs_by_size, 2);
    double total_temp = xt::sum<double>(buf_by_size)() / ((Nsamples - 1) * ndim);
    total_temp = std::pow(total_temp, 2);
    // for (ul ii = 0; ii < temps_by_size.shape()[0]; ++ii)
    // {
    //     if (temps_by_size[ii]!=0){
    //         std::cout << ii << ":" << temps_by_size[ii] << std::endl;
    //     }
    // }
    return total_temp;
}

struct runner2
{
    int ndim = 3;
    ul Nsamples{}, Nclusters{}, size{}, ind{};
    xt::xarray<ul> unique_cluster_ids;
    std::map<ul, xt::xarray<bool>> particle_mask_by_cluster_id;
    xt::xarray<bool> cluster_mask_by_size, particle_mask_by_size;
    xt::xarray<double> kes, buf_by_size, temps_by_size;
    xt::xarray<ul> sizes_counts, unique_sizes_row, particle_counts_by_size, ndofs_by_size;
    double total_temp{};

    double run(const xt::xarray<double> &velocities, const xt::xarray<double> &masses, const xt::xarray<ul> &particle_ids, const xt::xarray<ul> &cluster_ids)
    {
        Nsamples = particle_ids.shape()[0];

        unique_cluster_ids = xt::unique(cluster_ids);
        Nclusters = unique_cluster_ids.shape()[0];

        cluster_mask_by_size = xt::xarray<bool>(xt::xarray<bool>::shape_type({Nsamples + 1, Nclusters}), false);
        particle_mask_by_size = xt::xarray<bool>(xt::xarray<bool>::shape_type({Nsamples + 1, Nsamples}), false);

        size = 0;
        ind = 0;
        for (ul i : unique_cluster_ids)
        {
            particle_mask_by_cluster_id.emplace(i, xt::equal(cluster_ids, i));
            size = xt::sum<ul>(particle_mask_by_cluster_id[i])();
            cluster_mask_by_size(size, ind) = true;
            xt::view(particle_mask_by_size, size, xt::all()) += particle_mask_by_cluster_id[i];
            ++ind;
        }
        sizes_counts = xt::sum<ul>(cluster_mask_by_size, {1});
        kes = masses * velocities / 2;
        buf_by_size = xt::sum<double>(xt::broadcast(kes, {Nsamples + 1, Nsamples}) * particle_mask_by_size, {1});
        unique_sizes_row = xt::arange<ul>(0, Nsamples + 1);
        particle_counts_by_size = unique_sizes_row * sizes_counts;
        ndofs_by_size = (particle_counts_by_size - 1) * ndim;
        temps_by_size = xt::pow(buf_by_size / ndofs_by_size, 2);
        total_temp = std::pow(xt::sum<double>(buf_by_size)() / ((Nsamples - 1) * ndim), 2);
        return total_temp;
    }
    void clear()
    {
        particle_mask_by_cluster_id.clear();
    }
};

double run2(const xt::xarray<double> &velocities, const xt::xarray<double> &masses, const xt::xarray<ul> &particle_ids, const xt::xarray<ul> &cluster_ids)
{
    int ndim = 3;
    ul Nsamples = particle_ids.shape()[0];

    xt::xarray<ul> unique_cluster_ids(xt::unique(cluster_ids));
    ul Nclusters = unique_cluster_ids.shape()[0];

    std::map<ul, xt::xarray<bool>> particle_mask_by_cluster_id;
    xt::xarray<bool> cluster_mask_by_size(xt::xarray<bool>::shape_type({Nsamples + 1, Nclusters}), false);
    xt::xarray<bool> particle_mask_by_size(xt::xarray<bool>::shape_type({Nsamples + 1, Nsamples}), false);

    ul size = 0;
    ul ind = 0;
    for (const ul i : unique_cluster_ids)
    {
        particle_mask_by_cluster_id.emplace(i, xt::equal(cluster_ids, i));
        size = xt::sum<ul>(particle_mask_by_cluster_id[i])();
        cluster_mask_by_size(size, ind) = true;
        xt::view(particle_mask_by_size, size, xt::all()) += particle_mask_by_cluster_id[i];
        ++ind;
    }
    xt::xarray<ul> sizes_counts = xt::sum<ul>(cluster_mask_by_size, {1});
    xt::xarray<double> kes = masses * velocities / 2;
    xt::xarray<double> buf_by_size(xt::sum<double>(xt::broadcast(kes, {Nsamples + 1, Nsamples}) * particle_mask_by_size, {1}));
    xt::xarray<ul> unique_sizes_row = xt::arange<ul>(0, Nsamples + 1);
    xt::xarray<ul> particle_counts_by_size = unique_sizes_row * sizes_counts;
    xt::xarray<ul> ndofs_by_size = (particle_counts_by_size - 1) * ndim;
    xt::xarray<double> temps_by_size = xt::pow(buf_by_size / ndofs_by_size, 2);
    double total_temp = xt::sum<double>(buf_by_size)() / ((Nsamples - 1) * ndim);
    total_temp = std::pow(total_temp, 2);
    return total_temp;
    return 0;
}

template <typename T>
struct measurer
{
    T worker;
    std::chrono::system_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::chrono::system_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    xt::xarray<long double> vc1;
    ul Ncycles, Nloop;

    measurer(T &worker, ul Nloop, ul Ncycles) : worker(worker), Ncycles(Ncycles), Nloop(Nloop), vc1({Nloop}, 0.0) {}

    template <typename... Args>
    void measure(const Args &...);
};

template <typename T>
template <typename... Args>
void measurer<T>::measure(const Args &... args)
{
    for (int i = 0; i < Nloop; ++i)
    {
        t1 = std::chrono::high_resolution_clock::now();
        for (int j = 0; j < Ncycles; ++j)
        {
            worker.run(args...);
            worker.clear();
        }
        t2 = std::chrono::high_resolution_clock::now();
        vc1(i) = std::chrono::duration<long double, std::milli>(t2 - t1).count() / Ncycles;
    }
    auto vc1_mean = xt::mean(vc1);
    auto vc1_dev = xt::sqrt(xt::sum(xt::pow(vc1 - vc1_mean, 2)) / ((Nloop - 1) * Nloop));
    std::cout << "Mean: " << vc1_mean << " ms, dev: " << vc1_dev << " ms" << std::endl;
}

template <typename RFT, typename... Args>
struct func_measurer_t
{
    std::chrono::system_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::chrono::system_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    xt::xarray<long double> vc1;
    ul Ncycles, Nloop;

    RFT(*func)
    (const Args &...);

    func_measurer_t(RFT (*func)(const Args &...), ul Nloop, ul Ncycles) : func(func), Ncycles(Ncycles), Nloop(Nloop), vc1({Nloop}, 0.0) {}

    void measure(const Args &...args)
    {
        for (int i = 0; i < Nloop; ++i)
        {
            t1 = std::chrono::high_resolution_clock::now();
            for (int j = 0; j < Ncycles; ++j)
            {
                func(args...);
            }
            t2 = std::chrono::high_resolution_clock::now();
            vc1(i) = std::chrono::duration<long double, std::milli>(t2 - t1).count() / Ncycles;
        }
        auto vc1_mean = xt::mean(vc1);
        auto vc1_dev = xt::sqrt(xt::sum(xt::pow(vc1 - vc1_mean, 2)) / ((Nloop - 1) * Nloop));
        std::cout << "Mean: " << vc1_mean << " ms, dev: " << vc1_dev << " ms" << std::endl;
    }
};

int main(int argc, const char **argv)
{
    using std::chrono::duration;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    ul Nsamples = 1000000;
    const xt::xtensor<double, 1> velocities = xt::random::randn<double>({Nsamples}, 3.0, 1);
    const xt::xtensor<double, 1> masses({Nsamples}, 1.);
    const xt::xtensor<ul, 1> particle_ids = xt::random::randint<ul>({Nsamples}, 0, 15);
    const xt::xtensor<double, 1> cl_id_buf = xt::random::randn<double>({Nsamples}, 0.2 * Nsamples, 2);
    // const xt::xarray<ul> cluster_ids = xt::random::randn<ul>({Nsamples}, static_cast<ul>(Nsamples * 0.2));
    const xt::xtensor<ul, 1> cluster_ids = xt::round(cl_id_buf - xt::amin(cl_id_buf)());
    // xt::xarray<ul> cluster_ids = {1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7};

    func_measurer_t measurer_f1(&run, 10, 5);
    // func_measurer_t measurer_f2(&run2, 70, 70);
    double a = run(velocities, masses, particle_ids, cluster_ids);
    std::cout << a << std::endl;
    measurer_f1.measure(velocities, masses, particle_ids, cluster_ids);
    // measurer_f2.measure(velocities, masses, particle_ids, cluster_ids);

    const std::vector<double> _velocities(velocities.begin(), velocities.end());
    const std::vector<double> _masses(masses.begin(), masses.end());
    const std::vector<ul> _particle_ids(particle_ids.begin(), particle_ids.end());
    const std::vector<ul> _cluster_ids(cluster_ids.begin(), cluster_ids.end());
    double b = runc(_velocities, _masses, _particle_ids, _cluster_ids);
    std::cout << b << std::endl;
    func_measurer_t m(&runc, 10, 5);
    m.measure(_velocities, _masses, _particle_ids, _cluster_ids);

    runner1 run_1;
    runcer runcc(Nsamples);
    // runner2 run_2;
    measurer meas1(run_1, 10, 5);
    measurer measc(runcc, 10, 5);
    // measurer meas2(run_2, 70, 7);
    meas1.measure(velocities, masses, particle_ids, cluster_ids);
    measc.measure(_velocities, _masses, _particle_ids, _cluster_ids);
    // meas2.measure(velocities, masses, particle_ids, cluster_ids);

    return 0;
}