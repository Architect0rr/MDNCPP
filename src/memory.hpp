#ifndef __MDN_MEMORY__
#define __MDN_MEMORY__

#include <unordered_set>
#include <unordered_map>
#include <span>

#include "nlohmann/json.hpp"

#include "enums.hpp"

namespace mdn{
    extern spdlog::logger logger;
    using json = nlohmann::json;

    // naive stride-view implementation
    template <std::ranges::view V>
    class stride_view : public std::ranges::view_interface<stride_view<V>> {
    private:
        V base_;
        std::ranges::range_difference_t<V> stride_;

    public:
        // stride_view() = default;

        constexpr stride_view(V base, std::ranges::range_difference_t<V> stride): base_(std::move(base)), stride_(stride) {
            if (stride_ <= 0) { throw std::invalid_argument("Stride must be positive"); }
        }

        inline constexpr V base() const & { return base_; }
        inline constexpr V base() && { return std::move(base_); }

        inline constexpr auto begin() const {
            return iterator{std::ranges::begin(base_), stride_};
        }

        inline constexpr auto end() const {
            auto end_offset = std::distance(std::ranges::begin(base_), std::ranges::end(base_)) % stride_ - stride_;
            return iterator{std::ranges::end(base_) - end_offset, stride_};
        }

        inline constexpr auto size() const requires std::ranges::sized_range<V> {
            auto n = std::ranges::size(base_);
            return n / stride_ + (n % stride_ != 0);
        }

        // element access
        const auto& operator[](std::ranges::range_difference_t<V> n) const {
            auto pos = std::ranges::begin(base_) + n * stride_;
            // if (pos < std::ranges::begin(base_) || pos >= std::ranges::end(base_)) {
            //     throw std::out_of_range("stride_view index out of range");
            // }
            return *pos;
        }

        auto at(std::ranges::range_difference_t<V> n) const {
            auto pos = std::ranges::begin(base_) + n * stride_;
            if (pos < std::ranges::begin(base_) || pos >= std::ranges::end(base_)) {
                throw std::out_of_range("stride_view index out of range");
            }
            // if (n < 0 || n >= size()) {  // size() should be available
            //     throw std::out_of_range("stride_view index out of range");
            // }
            // return (*this)[n]; // Delegate to operator[]
            return *pos;
        }

        class iterator {
        public:
            using iterator_category = std::random_access_iterator_tag;
            using value_type        = std::ranges::range_value_t<V>;
            using difference_type   = std::ranges::range_difference_t<V>;
            using pointer           = value_type*;
            using reference         = value_type&;

            // Constructors
            iterator() = default;
            constexpr iterator(std::ranges::iterator_t<V> current, difference_type stride) : current_(current), stride_(stride) {}

            // Copyable and movable
            iterator(const iterator&) = default;
            iterator& operator=(const iterator&) = default;
            iterator(iterator&&) = default;
            iterator& operator=(iterator&&) = default;

            // Dereference operators
            inline constexpr reference operator*() const { return *current_; }
            inline constexpr pointer operator->() const { return &*current_; }

            // Increment/Decrement operators
            inline constexpr iterator& operator++() {
                std::advance(current_, stride_);
                return *this;
            }
            inline constexpr iterator operator++(int) {
                auto temp = *this;
                ++(*this);
                return temp;
            }
            inline constexpr iterator& operator--() {
                std::advance(current_, -stride_);
                return *this;
            }
            inline constexpr iterator operator--(int) {
                auto temp = *this;
                --(*this);
                return temp;
            }

            // Arithmetic operators
            inline constexpr iterator operator+(difference_type n) const {
                return iterator(current_ + n * stride_, stride_);
            }
            inline constexpr iterator& operator+=(difference_type n) {
                current_ += n * stride_;
                return *this;
            }
            inline constexpr iterator operator-(difference_type n) const {
                return iterator(current_ - n * stride_, stride_);
            }
            inline constexpr iterator& operator-=(difference_type n) {
                current_ -= n * stride_;
                return *this;
            }
            inline friend constexpr iterator operator+(difference_type n, const iterator& it) {
                return it + n * it.stride_;
            }
            inline friend constexpr difference_type operator-(const iterator& lhs, const iterator& rhs) {
                return (lhs.current_ - rhs.current_) / lhs.stride_;
            }

            // Comparison operators
            inline friend constexpr bool operator==(const iterator& lhs, const iterator& rhs) {
                return lhs.current_ == rhs.current_;
            }

        private:
            std::ranges::iterator_t<V> current_;
            difference_type stride_;
        };  // iterator
    };  // stride_view

    // Helper function to create stride_view
    template<std::ranges::viewable_range R>
    stride_view(R&&, std::ranges::range_difference_t<std::views::all_t<R>>) -> stride_view<std::views::all_t<R>>;



    class Memory{
    public:
        static Memory construct(const std::unordered_map<int, PROPS>& _table){
            std::unordered_set<PROPS> avails;
            std::unordered_map<PROPS, int> reverse_table;
            logger.debug("---Parsed properties:");
            for (const auto& [k, v] : _table) {
                avails.emplace(v);
                reverse_table.emplace(v, k);
                logger.debug("    {}->{}. Reversed: {}->{}", k, static_cast<int>(v), static_cast<int>(v), reverse_table[v]);
            }
            return Memory(_table, reverse_table, avails);
        }

        static std::unordered_map<int, PROPS> parse_props(const json& j) {
            std::unordered_map<int, PROPS> result;

            if (!j.is_object()) throw std::runtime_error("PROPS JSON is not an object");

            for (const auto& [key, value] : j.items())
                result.emplace(std::stoi(key), PROPS(value.get<int>()));
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

        auto prop(const PROPS& _prop){
            if (!reverse_table.contains(_prop)) throw std::runtime_error("Invalid property");
            // return _span | std::views::transform([offset = this->reverse_table[_prop], nprops = this->_nprops](size_t i) { return offset + i * nprops; });
            // return _span | std::views::filter([offset = this->reverse_table[_prop], nprops = this->_nprops](size_t i)
            //     { return static_cast<size_t>(static_cast<size_t>(i - offset) % nprops) == 0; });

            return stride_view(std::ranges::drop_view(_span, reverse_table[_prop]), _nprops);
        }

    private:
        Memory(const std::unordered_map<int, PROPS>& _table, const std::unordered_map<PROPS, int>& _reverse_table, const std::unordered_set<PROPS>& _avails):
            table(_table), reverse_table(_reverse_table), _nprops(_table.size()), avails(_avails), valid(true) { }
        Memory(const Memory&) = delete;
        Memory& operator=(const Memory&) = delete;

        bool valid = false;
        int _nprops = 0, nlines_allocated = 0;
        std::unordered_map<int, PROPS> table;
        std::unordered_map<PROPS, int> reverse_table;
        std::unordered_set<PROPS> avails;

        std::unique_ptr<double[]> _data;
        std::span<const double> _span;
    };

}  // namepace mdn

#endif // !__MDN_MEMORY__
