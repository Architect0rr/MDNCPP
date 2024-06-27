#ifndef __MDN_MEMORY__
#define __MDN_MEMORY__

#include <unordered_set>
#include <unordered_map>
#include <span>

#include "nlohmann/json.hpp"

// #include "enums.hpp"
#include "ct_map.hpp"
#include "stride_view.hpp"

namespace mdn{
    extern spdlog::logger logger;
    using json = nlohmann::json;

    class Memory{
    public:
        static Memory construct(const std::unordered_map<int, CT>& _table){
            std::unordered_set<CT> avails;
            std::unordered_map<CT, int> reverse_table;
            logger.debug("---Parsed properties:");
            for (const auto& [k, v] : _table) {
                avails.emplace(v);
                reverse_table.emplace(v, k);
                logger.debug("    {}->{}. Reversed: {}->{}", k, static_cast<int>(v), static_cast<int>(v), reverse_table[v]);
            }
            return Memory(_table, reverse_table, avails);
        }

        static std::unordered_map<int, CT> parse_props(const json& j) {
            std::unordered_map<int, CT> result;

            if (!j.is_object()) throw std::runtime_error("PROPS JSON is not an object");

            for (const auto& [key, value] : j.items())
                result.emplace(std::stoi(key), CT(value.get<int>()));
            return result;
        }

        Memory(): valid(false) {}

        Memory(Memory&& o):
            valid(o.valid), _nprops(o._nprops), nlines_allocated(o.nlines_allocated),
            table(std::move(o.table)), reverse_table(std::move(o.reverse_table)), avails(std::move(o.avails)), _data(std::move(o._data))
        {
            o.valid = false;
            o._nprops = 0;
            o.nlines_allocated = 0;
        }

        Memory& operator=(Memory&& o){
            valid = o.valid;
            o.valid = false;
            _nprops = o._nprops;
            o._nprops = 0;
            nlines_allocated = o.nlines_allocated;
            o.nlines_allocated = 0;
            table = std::move(o.table);
            reverse_table = std::move(o.reverse_table);
            avails = std::move(o.avails);
            _data = std::move(o._data);
            return *this;
        }

        inline double* get() {return _data.get();}
        inline const double* const get() const {return _data.get();}

        inline double* data() {
            if (valid && nlines_allocated) return _data.get();
            else throw std::runtime_error("Invalid object state");
        }

        inline const double* const data() const {
            if (valid && nlines_allocated) return _data.get();
            else throw std::runtime_error("Invalid object state");
        }

        void allocate(const uint64_t& count){
            _data = std::unique_ptr<double[]>(new double[count*_nprops]);
            _span = std::span<const double>(_data.get(), count*_nprops);
            nlines_allocated = count;
        }

        inline const int nprops() {return _nprops;}

        auto prop(const CT& _prop){
            if (!reverse_table.contains(_prop)) throw std::runtime_error("Invalid property");
            // return _span | std::views::transform([offset = this->reverse_table[_prop], nprops = this->_nprops](size_t i) { return offset + i * nprops; });
            // return _span | std::views::filter([offset = this->reverse_table[_prop], nprops = this->_nprops](size_t i)
            //     { return static_cast<size_t>(static_cast<size_t>(i - offset) % nprops) == 0; });

            return stride_view(std::ranges::drop_view(_span, reverse_table[_prop]), _nprops);
        }

        inline bool contains(const CT& _prop) {return reverse_table.contains(_prop); }

    private:
        Memory(const std::unordered_map<int, CT>& _table, const std::unordered_map<CT, int>& _reverse_table, const std::unordered_set<CT>& _avails):
            table(_table), reverse_table(_reverse_table), _nprops(_table.size()), avails(_avails), valid(true) { }
        Memory(const Memory&) = delete;
        Memory& operator=(const Memory&) = delete;

        bool valid = false;
        int _nprops = 0, nlines_allocated = 0;
        std::unordered_map<int, CT> table;
        std::unordered_map<CT, int> reverse_table;
        std::unordered_set<CT> avails;

        std::unique_ptr<double[]> _data;
        std::span<const double> _span;
    };

}  // namepace mdn

#endif // !__MDN_MEMORY__
