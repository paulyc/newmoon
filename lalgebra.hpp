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

#include <functional>
#include <string>
#include <cstdint>
#include <cmath>

///static constexpr binary128 MMMMMM_PI = 0x3.243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89p0;
static constexpr long double MMM_PI = 0x3.243f6a8885a308d313198a2e03707344p0;
//typedef double vec3d_t[3];
//typedef std::array<long double, 2> vec2q_t;
//typedef std::array<long double, 3> vec3q_t;

template <typename T=long double>
struct vec2
{
	typedef T TT;
	T raw[2];
	long double dotP(const vec2<T> &v) const {
		return raw[0] * v.raw[0] + raw[1] * v.raw[1];
	}
	long double mag() const {
		return sqrtl(dotP(*this));
	}
	long double phase() const {
		return atanl(raw[1]/raw[0]);
	}
	vec2<T> sum(const vec2<T> &v) const {
		return {raw[0]+v.raw[0], raw[1]+v.raw[1]};
	}
	vec2<long double> normalize() const {
		const long double mag = this->mag();
		return {raw[0]/mag, raw[1]/mag};
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
	typedef T TT;
	T raw[3];

	vec2<T> drop_z() const {
		return {raw[0], raw[1]};
	}
	long double dotP(const vec3<T> &v) const {
		return raw[0] * v.raw[0] + raw[1] * v.raw[1] + raw[2] * v.raw[2];
	}
	long double mag() const {
		return sqrtl(dotP(*this));
	}
	// not really such a thing in spherical coordinates, but there is
	// if we ignore z and pretend it's cylindrical, which is generally
	// a Good Enough Approximation(tm) of a planet-star-type orbit
	long double phase() const {
		return atan2l(raw[1], raw[0]);
	}
	vec3<T> normalize() const {
		const long double mag = this->mag();
		return {raw[0]/mag, raw[1]/mag, raw[2]/mag};
	}
};
typedef vec3<double> vec3d_t;
typedef vec3<long double> vec3q_t;

template <typename T=long double>
struct mat3x3
{
	typedef T TT;
	T elems[3][3];

	// https://gssc.esa.int/navipedia/index.php/Transformation_between_Terrestrial_Frames
	// ????
	vec3q_t mul(const vec3q_t &v) const {
		return {
			elems[0][0] * v.raw[0] + elems[0][1] * v.raw[1] + elems[0][2] * v.raw[2],
			elems[1][0] * v.raw[0] + elems[1][1] * v.raw[1] + elems[1][2] * v.raw[2],
			elems[2][0] * v.raw[0] + elems[2][1] * v.raw[1] + elems[2][2] * v.raw[2],
		};
	}
	static constexpr mat3x3 R_0(long double α, long double θ_1, long double θ_2, long double θ_3) {
		//const mat3x3 r_1 = R_1(θ_1);
		//const mat3x3 r_2 = R_2(θ_2);
		//const mat3x3 r_3 = R_3(θ_3);
		return {
			{   α, -θ_3,  θ_2},
			{ θ_3,    α, -θ_1},
			{-θ_2,  θ_1,    α},
		};
	}
	static constexpr mat3x3 R_1(long double θ) {
		const long double sin_θ = sinl(θ);
		const long double cos_θ = cosl(θ);
		return {
			{1.0l,   0.0l,  0.0l},
			{0.0l,  cos_θ, sin_θ},
			{0.0l, -sin_θ, cos_θ},
		};
	}
	static constexpr mat3x3 R_2(long double θ) {
		const long double sin_θ = sinl(θ);
		const long double cos_θ = cosl(θ);
		return {
			{cos_θ, 0.0l, -sin_θ},
			{ 0.0l, 1.0l,   0.0l},
			{sin_θ, 0.0l,  cos_θ},
		};
	}
	static constexpr mat3x3 R_3(long double θ) {
		const long double sin_θ = sinl(θ);
		const long double cos_θ = cosl(θ);
		return {
			{ cos_θ, sin_θ, 0.0l},
			{-sin_θ, cos_θ, 0.0l},
			{  0.0l,  0.0l, 1.0l},
		};
	}
};

typedef mat3x3<long double> mat3x3q_t;

template <typename T=long double>
struct funmat3x3
{
	typedef T TT;
	typedef std::function<TT(TT)> matfun;
	static constexpr matfun ident = [](const TT &t) -> TT { return t; };
	static constexpr matfun zero = [](const TT &t) -> TT { return TT(0.0); };
	static constexpr matfun one = [](const TT &t) -> TT { return TT(1.0); };
	static constexpr matfun sine = [](const TT &θ) -> TT { return sinl(θ); };
	static constexpr matfun cosine = [](const TT &θ) -> TT { return cosl(θ); };
	static constexpr auto sincos = [](const TT &θ) -> vec2q_t { return {sinl(θ), cosl(θ)}; }; // premature optimization

	matfun elems[3][3] = {
		{ident, ident, ident},
		{ident, ident, ident},
		{ident, ident, ident},
	};
};

template <typename Vec_T>
struct coordspace {};

template <typename Mat_T>
struct matxfrm {};

template <typename Vec_T, typename Mat_T>
struct spacexfrm {};

struct cartesian2dspace : public coordspace<vec2q_t>
{
	struct point : public vec2q_t
	{
		TT x() const { return this->raw[0]; }
		TT y() const { return this->raw[1]; }
	};
};

struct cylindrical2dspace : public coordspace<vec2q_t>
{
	struct point : public vec2q_t
	{
		TT r() const { return this->raw[0]; }
		TT θ() const { return this->raw[1]; }
	};
};

struct spacexfrm2d : public spacexfrm<vec2q_t, mat3x3q_t>
{
	cylindrical2dspace::point operator()(const cartesian2dspace::point &p) {
		return cylindrical2dspace::point {p.mag(), p.phase()};
	}
	cartesian2dspace::point operator()(const cylindrical2dspace::point &p) {
		return cartesian2dspace::point {p.r() * cosl(p.θ()), p.r() * sinl(p.θ())};
	}
};

struct cartesian3dspace : public coordspace<vec3q_t>
{
	struct point : public vec3q_t
	{
		TT x() const { return this->raw[0]; }
		TT y() const { return this->raw[1]; }
		TT z() const { return this->raw[2]; }
	};
};

struct spherical3dspace : public coordspace<vec3q_t>
{
	struct point : public vec3q_t
	{
		TT r() const { return this->raw[0]; }
		TT θ() const { return this->raw[1]; }
		TT ø() const { return this->raw[2]; }
	};
};


#endif /* _PAULYC_LALGEBRA_HPP_ */
