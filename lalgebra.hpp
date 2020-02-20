/*
* lalgebra.hpp - part of newmoon, moon phase calculator
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
*/

#ifndef _PAULYC_LALGEBRA_HPP_
#define _PAULYC_LALGEBRA_HPP_

#include <array>
#include <cstdint>
#include <cmath>

typedef double vec3d_t[3];
typedef vec3d_t pv_t[2];
typedef std::array<long double, 2> vec2q_t;
typedef std::array<long double, 3> vec3q_t;

static inline vec2q_t drop_z(const vec3d_t &v) {
    return {static_cast<long double>(v[0]), static_cast<long double>(v[1])};
}

static inline long double dotP(const vec2q_t &u, const vec2q_t &v)
{
    return u[0] * v[0] + u[1] * v[1];
}

static inline long double magn(const vec2q_t &u)
{
    return sqrtl(dotP(u,u));
}

static inline long double dotP(const vec3q_t &u, const vec3q_t &v)
{
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

static inline long double magn(const vec3q_t &u)
{
    return sqrtl(dotP(u,u));
}

static inline vec2q_t normalize(const vec2q_t &v)
{
    const long double mag = magn(v);
    return {v[0]/mag,v[1]/mag};
}

static inline vec3q_t normalize(const vec3d_t &v)
{
    const long double x = static_cast<long double>(v[0]);
    const long double y = static_cast<long double>(v[1]);
    const long double z = static_cast<long double>(v[2]);
    const long double mag = sqrtl(x*x+y*y);
    return {x/mag,y/mag,z};
}

static inline long double arg(const vec2q_t &u, const vec2q_t &v)
{
    return acosl(dotP(u,v)/(magn(u) * magn(v)));
}

#endif /* _PAULYC_LALGEBRA_HPP_ */
