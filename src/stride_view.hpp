#ifndef __MDN_STRIDE_VIEW__
#define __MDN_STRIDE_VIEW__

#include <ranges>

namespace mdn{
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

} // namespace mdn

#endif // !__MDN_STRIDE_VIEW__