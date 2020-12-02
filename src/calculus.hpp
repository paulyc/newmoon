/**
 * calculus.hpp
 *
 * Copyright (C) 2020 Paul Ciarlo <paul.ciarlo@gmail.com>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
 **/

#ifndef PAULYC_CALCULUS_HPP
#define PAULYC_CALCULUS_HPP

#include <functional>
#include <optional>
#include <cstdint>
#include <cmath>

#include "quadmath.h"
#include "lalgebra.hpp"

namespace github {
namespace paulyc {

typedef std::function<__float128(__float128)> fun_1d_t;
typedef std::function<__float128(__float128, __float128)> fun_2d_t;
typedef std::function<__float128(__float128, __float128, __float128)> fun_3d_t;
typedef std::function<__float128(__float128, __float128, __float128, __float128)> fun_4d_t;

// identity
fun_1d_t id_op(fun_1d_t fun, const __float128 delta=FLT128_EPSILON) {
    return [fun, delta](__float128 x) -> __float128 {
        return fun(x);
    };
}

//derivative
fun_1d_t d_op(fun_1d_t fun, const __float128 delta=FLT128_EPSILON) {
    const __float128 _2_delta_m1 = 0.5q / delta;
    return [fun, delta, _2_delta_m1](__float128 x) -> __float128 {
        return _2_delta_m1 * (fun(x+delta) - fun(x-delta));
    };
}

//2nd derivative
fun_1d_t d2_op(fun_1d_t fun, const __float128 delta=FLT128_EPSILON) {
    const __float128 _2_delta_m1 = 0.5q / delta;
    return d_op(d_op(fun, delta), delta);
}

std::optional<__float128> min_x(fun_1d_t fun, const __float128 range_min, const __float128 range_max, const __float128 delta=FLT128_EPSILON) {
    auto d_fun = d_op(fun, delta);
    auto d2_fun = d_op(fun, delta);
    // find d_fun zero crossing
    int d_fun_sign = signbitq(d_fun(range_min));
    for (__float128 x = range_min; x < range_max; x += delta) {
        if (d_fun_sign != signbitq(d_fun(x)) && signbitq(d2_fun(x)) > 0) {
            return x;
        }
    }
    return std::nullopt;
}

} /* namespace paulyc */
} /* namespace github */

#endif /* PAULYC_CALCULUS_HPP */
