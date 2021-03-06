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

typedef std::function<long double(long double)> fun_1d_t;
typedef std::function<long double(long double, long double)> fun_2d_t;
typedef std::function<long double(long double, long double, long double)> fun_3d_t;
typedef std::function<long double(long double, long double, long double, long double)> fun_4d_t;

// identity
fun_1d_t id_op(fun_1d_t fun, const long double delta);

//derivative
fun_1d_t d_op(fun_1d_t fun, const long double delta);

//2nd derivative
fun_1d_t d2_op(fun_1d_t fun, const long double delta);

std::optional<long double> min_x(fun_1d_t fun, const long double range_min, const long double range_max, const long double delta);

} /* namespace paulyc */
} /* namespace github */

#endif /* PAULYC_CALCULUS_HPP */
