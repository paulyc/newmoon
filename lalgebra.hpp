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

#include <cstdint>
#include <cmath>

//typedef double vec3d_t[3];
//typedef std::array<long double, 2> vec2q_t;
//typedef std::array<long double, 3> vec3q_t;

template <typename T=long double>
struct vec2
{
	T x;
	T y;
	long double dotP(const vec2<T> &v) const {
		return x * v.x + y * v.y;
	}
	long double mag() const {
		return sqrtl(dotP(*this));
	}
	long double phase() const {
		return atanl(y/x);
	}
	vec2<T> sum(const vec2<T> &v) const {
		return {x+v.x, y+v.y};
	}
	vec2<long double> normalize() const {
		const long double mag = this->mag();
		return {x/mag, y/mag};
	}
	vec2<long double> magphase() const {
		return {this->mag(), this->phase()};
	}
	long double angle(const vec2<T> &v) const {
		return acosl(dotP(v)/(mag() * v.mag()));
	}
};
typedef vec2<double> vec2d_t;
typedef vec2<long double> vec2q_t;

template <typename T=long double>
struct vec3
{
	T x;
	T y;
	T z;

	vec2<T> drop_z() const {
		return {x, y};
	}
	long double dotP(const vec3<T> &v) const {
		return x * v.x + y * v.y + z * v.z;
	}
	long double mag() const {
		return sqrtl(dotP(*this));
	}
	long double phase() const {
		return atanl(y/x);
	}
	vec3<T> normalize() const {
		const long double mag = this->mag();
		return {x/mag, y/mag, z/mag};
	}
};
typedef vec3<double> vec3d_t;
typedef vec3<long double> vec3q_t;
typedef vec3d_t pv_t[2];

template <typename T=long double>
struct mat3x3
{
	T elems[3][3];
};

#endif /* _PAULYC_LALGEBRA_HPP_ */
