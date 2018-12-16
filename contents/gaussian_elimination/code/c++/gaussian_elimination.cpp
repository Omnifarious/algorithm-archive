#include <array>
#include <cstdlib>
#include <fmt/format.h>
#include <gmpxx.h>
#include <utility>

namespace fmt {

template <>
struct formatter<mpq_class> {
   template <typename ParseContext>
   constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

   template <typename FormatContext>
   auto format(const mpq_class &rat, FormatContext &ctx) {
      mpz_class const &num = rat.get_num();
      mpz_class const &denom = rat.get_den();
      if (denom == 1) {
         return format_to(ctx.begin(), "{:s}", num.get_str());
      } else {
         return format_to(ctx.begin(), "{:s}/{:s}", num.get_str(), denom.get_str());
      }
   }
};

} // namespace fmt

using ary_size_t = ::std::array<int, 1>::size_type;
template <typename T, ary_size_t size, ary_size_t augments> using GJMatrix =
      ::std::array< ::std::array<T, size + augments>, size>;


template <class T, ary_size_t rows, ary_size_t cols>
constexpr ::std::pair<ary_size_t, ary_size_t>
get_array_rows_cols(::std::array< ::std::array<T, cols>, rows> *)
{
   return {rows, cols};
}


template <typename T>
void guassian_elimination(T &A)
{
    using row_t = typename T::value_type;
    using val_t = typename T::value_type::value_type;
    constexpr auto size = get_array_rows_cols((T *)(nullptr)).first;
    constexpr auto total_width = get_array_rows_cols((T *)(nullptr)).second;
    constexpr auto augments = total_width - size;


    auto max_col_row = [&A, &size](ary_size_t col, ary_size_t rowstart) {
        using ::std::abs;
        ary_size_t tmp_max = rowstart;
        val_t tmp_max_val = abs(A[tmp_max][col]);
        ++rowstart;
        while (rowstart < size) {
            if (tmp_max_val < abs(A[rowstart][col])) {
                tmp_max_val = abs(A[rowstart][col]);
                tmp_max = rowstart;
            }
            ++rowstart;
        }
        return tmp_max;
    };

    for (ary_size_t pivot_col = 0; pivot_col < size; ++pivot_col) {
        {
            ary_size_t const max_row = max_col_row(pivot_col, pivot_col);
            if (max_row != pivot_col) {
                A[max_row].swap(A[pivot_col]);
            }
        }
        // Now pivot_col = pivot_row.
        val_t const denom = A[pivot_col][pivot_col];
        for (ary_size_t elim_row = pivot_col + 1; elim_row < size; ++elim_row) {
            val_t const numerator = A[elim_row][pivot_col];
            for (ary_size_t fix_col = pivot_col + 1; fix_col < total_width; ++fix_col) {
                A[elim_row][fix_col] -= (A[pivot_col][fix_col] * numerator) / denom;
            }
            A[elim_row][pivot_col] = 0;
        }
        for (ary_size_t fix_col = pivot_col + 1; fix_col < total_width; ++fix_col) {
            A[pivot_col][fix_col] /= denom;
        }
        A[pivot_col][pivot_col] = 1;
    }
}

int main()
{
    GJMatrix<mpq_class, 3, 1> tmp{ ::std::array<mpq_class, 4>{ 2, 3, 4, 6 },
                                   ::std::array<mpq_class, 4>{ 1, 2, 3, 4 },
                                   ::std::array<mpq_class, 4>{ 3, -4, 0, 10} };
    guassian_elimination(tmp);
    for (auto const &row: tmp) {
       for (auto const &val: row) {
          fmt::print("{} ", val);
       }
       fmt::print("\n");
    }
    return 0;
}
