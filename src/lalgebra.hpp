/**
 * part of newmoon, moon phase calculator
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

#ifndef _PAULYC_LALGEBRA_HPP_
#define _PAULYC_LALGEBRA_HPP_

#include <functional>
#include <string>
#include <cstdint>
#include <cmath>

// magic number copied from somewhere(??). Not eligible for copyright protection as a simple statement of fact.
// this is really a binary256 but who tf knows how many bits you'll ever get in a long double anyway,
// it gets truncated to 80 bits on some i386/amd64 platforms, which makes it something of a short long double since
// the chip does all double FP in 80-bit anyway and rounds it back off to 64.
// who knows, maybe some will whole-ass it with all 256, or make a long long double later.
static constexpr long double MMM_PI    = 0x3.243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89p0;
// might come in handy for solving the QM diffEQs or something idk but it's fun
static constexpr long double MMM_PI_PI = 0x3.243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89p1;
static constexpr long double MMM_2_PI  = 2.0l * MMM_PI;

template <std::size_t N, typename T>
struct Nvec
{
	typedef Nvec<N, T> ThisT;
	typedef T PrimT;
	typedef PrimT PrimArrayT[N];
	PrimArrayT data;

	constexpr Nvec() {
		for (std::size_t i = 0; i < N; ++i) {
			data[i] = T(0.0l);
		}
	}
	constexpr Nvec(T initval) {
		for (std::size_t i = 0; i < N; ++i) {
			data[i] = initval;
		}
	}

	PrimT& operator[](std::size_t i) {
		if (i < N) {
			return data[i];
		} else {
			throw std::range_error("OOR");
		}
	}

	constexpr Nvec operator+(const Nvec &that) const {
		Nvec sum;
		for (std::size_t i = 0; i < N; ++i) {
			sum[i] = (*this)[i] + that[i];
		}
		return sum;
	}

	Nvec& operator+=(const Nvec &that) {
		for (std::size_t i = 0; i < N; ++i) {
			(*this)[i] += that[i];
		}
		return *this;
	}

	constexpr Nvec operator-() const {
		Nvec difference;
		for (std::size_t i = 0; i < N; ++i) {
			difference[i] = -(*this)[i];
		}
		return difference;
	}

	constexpr T operator*(const Nvec &that) const {
		T product;
		for (std::size_t i = 0; i < N; ++i) {
			product += (*this)[i] * that[i];
		}
		return product;
	}
};

template <typename T>
struct Nvec<0,T> {
	constexpr T dotP(const Nvec<0, T> &that) const {
		return T(0.0l);
	}
};
/*
template <typename T>
struct Nvec<3,T> {
	constexpr T crossP(const Nvec<3, T> &that) const {
		Nvec product;

	}
};
*/
typedef Nvec<0, long double> v0q_t;
typedef Nvec<1, long double> v1q_t;
typedef Nvec<2, long double> v2q_t;
typedef Nvec<3, long double> v3q_t;

template <typename T>
static const Nvec<0, T> TheZeroVector;

template <std::size_t N, std::size_t M, typename T=long double>
struct NxMmatrix
{
	typedef T TT;
	typedef Nvec<N, T> NT;
	typedef Nvec<M, T> MT;
	typedef Nvec<N, MT> MVecT;
	MVecT data;
	MT& operator[](std::size_t n) {
		if (n < N) {
			 return (*this)[n];
		} else {
			throw std::range_error("OOR");
		}
	}
	NT Nvec(std::size_t m) const {
		if (m < M) {
			NT nvec;
			for (auto n = 0; n < N; ++n) {
				nvec[n] = (*this)[n][m];
			}
			return nvec;
		} else {
			throw std::range_error("OOR");
		}
	}
	NT operator*(const MT &v) const {
		NT product;
		for (int n = 0; n < N; ++n) {
			const MT &mvec = (*this)[n];
			for (int m = 0; m < M; ++m) {
				product[n] += mvec[m] * v[m];
			}
		}
		return product;
	}

	// rows in m must be equal to cols in *this but let's keep it simple
	NxMmatrix operator*(const NxMmatrix &that) const {
		NxMmatrix product;
		for (auto n = 0; n < N; ++n) {
			const MT &mvec = (*this)[n];
			for (auto m = 0; m < M; ++m) {
				product[n] += mvec[m] * that[m];
			}
		}
		return product;
	}
};
template <typename T=long double>
struct vec2
{
	typedef T TT;
	T raw[2];

	T dotP(const vec2<T> &v) const {
		return raw[0] * v.raw[0] + raw[1] * v.raw[1];
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
	T dotP(const vec3<T> &v) const {
		return raw[0] * v.raw[0] + raw[1] * v.raw[1] + raw[2] * v.raw[2];
	}
	T mag() const {
		return sqrtl(dotP(*this));
	}
};
typedef vec3<double> vec3d_t;
typedef vec3<long double> vec3q_t;

template <typename T=long double>
struct mat3x3
{
	typedef T TT;
	T elems[3][3];

	constexpr mat3x3() {
		for (std::size_t i = 0; i < 3; ++i) {
			for (std::size_t j = 0; i < 3; ++j) {
				elems[i][j] = T(0.0l);
			}
		}
	}

	constexpr mat3x3(const T (&m)[3][3]) {
		for (std::size_t i = 0; i < 3; ++i) {
			for (std::size_t j = 0; i < 3; ++j) {
				elems[i][j] = m[i][j];
			}
		}
	}
	constexpr mat3x3(const mat3x3<T> &that) {
		for (std::size_t i = 0; i < 3; ++i) {
			for (std::size_t j = 0; i < 3; ++j) {
				this->elems[i][j] = that.elems[i][j];
			}
		}
	}
	constexpr mat3x3& operator=(const mat3x3 &that) {
		for (std::size_t i = 0; i < 3; ++i) {
			for (std::size_t j = 0; i < 3; ++j) {
				this->elems[i][j] = that.elems[i][j];
			}
		}
		return *this;
	}

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

template <typename Vec_T>
struct coordspace {};

template <typename Mat_T>
struct xfrmmatrix {};

template <typename SpaceA_T, typename SpaceB_T, typename TransformationMatrixT=mat3x3q_t>
struct spacexfrm {};

struct cartesian2dvec : public vec2q_t
{
	constexpr TT x() const { return this->raw[0]; }
	constexpr TT y() const { return this->raw[1]; }
		TT mag() const {
		return sqrtl(dotP(*this));
	}
	TT phase() const {
		return atanl(y()/x());
	}
	cartesian2dvec sum(const cartesian2dvec &v) const {
		return {x()+v.x(), y()+v.y()};
	}
	cartesian2dvec normalize() const {
		const TT mag = this->mag();
		return {x()/mag, y()/mag};
	}
	TT angle(const cartesian2dvec &v) const {
		return acosl(dotP(v)/(mag() * v.mag()));
	}
};

struct cylindrical2dvec : public vec2q_t
{
	constexpr TT r() const { return this->raw[0]; }
	constexpr TT θ() const { return this->raw[1]; }
	constexpr TT mag() const {
		return r();
	}
	constexpr TT phase() const {
		return θ();
	}
	cylindrical2dvec sum(const cylindrical2dvec &v) const;
	constexpr vec2q_t normalize() const {
		return {1.0l, phase()};
	}
	TT angle(const cylindrical2dvec &v) const {
		const TT diff = fmodl(phase() - v.phase(), MMM_2_PI);
		if (diff < 0.0l) {
			return diff + MMM_2_PI;
		} else {
			return diff;
		}
	}
};

struct spacexfrm2d : public spacexfrm<cartesian2dvec, cylindrical2dvec, mat3x3q_t>
{
	static cylindrical2dvec cart2cyl(const cartesian2dvec &p) {
		return cylindrical2dvec {p.mag(), p.phase()};
	}
	static cartesian2dvec cyl2cart(const cylindrical2dvec &p) {
		return cartesian2dvec {p.r() * cosl(p.θ()), p.r() * sinl(p.θ())};
	}
};

cylindrical2dvec cylindrical2dvec::sum(const cylindrical2dvec &v) const
{
	return spacexfrm2d::cart2cyl(spacexfrm2d::cyl2cart(*this).sum(spacexfrm2d::cyl2cart(v)));
}

struct cartesian3dvec : public vec3q_t
{
	constexpr TT x() const { return this->raw[0]; }
	constexpr TT y() const { return this->raw[1]; }
	constexpr TT z() const { return this->raw[2]; }
	TT mag() const {
		return sqrtl(dotP(*this));
	}
	TT phase() const {
		return atanl(y()/x());
	}
	cartesian3dvec sum(const cartesian3dvec &v) const {
		return {x()+v.x(), y()+v.y(), z()+v.z()};
	}
	cartesian3dvec normalize() const {
		const TT mag = this->mag();
		return {x()/mag, y()/mag, z()/mag};
	}
	TT angle(const cartesian3dvec &v) const {
		return acosl(dotP(v)/(mag() * v.mag()));
	}
};

struct spherical3dvec : public vec3q_t
{
	constexpr TT r() const { return this->raw[0]; }
	constexpr TT θ() const { return this->raw[1]; }
	constexpr TT ø() const { return this->raw[2]; }
	constexpr TT mag() const {
		return r();
	}
	//spherical3dvec sum(const spherical3dvec &v) const;
	constexpr vec3q_t normalize() const {
		return {1.0l, θ(), ø()};
	}
	// not really such a thing in spherical coordinates, but there is
	// if we ignore z and pretend it's cylindrical, which is generally
	// a Good Enough Approximation(tm) of a planet-star-type orbit
	constexpr TT phase() const {
		return θ();
	}
	// ignoring ø for now keep it simple see above
	TT angle(const spherical3dvec &v) const {
		const TT diff = fmodl(phase() - v.phase(), MMM_2_PI);
		if (diff < 0.0l) {
			return diff + MMM_2_PI;
		} else {
			return diff;
		}
	}
};

struct spacexfrm3d : public spacexfrm<cartesian3dvec, spherical3dvec, mat3x3q_t>
{
	static spherical3dvec cart2sph(const cartesian3dvec &p) {
		return spherical3dvec {p.mag(), atan2l(sqrtl(p.x()*p.x()+p.y()*p.y()), p.z()), atan2l(p.y(), p.x())};
	}
	static cartesian3dvec sph2cart(const spherical3dvec &p) {
		const long double sin_θ = sinl(p.θ());
		const long double cos_θ = cosl(p.θ());
		const long double sin_ø = sinl(p.ø());
		const long double cos_ø = cosl(p.ø());
		// https://www.web-formulas.com/Math_Formulas/Linear_Algebra_Transform_from_Cartesian_to_Spherical_Coordinate.aspx
		const long double m[3][3] = {
			{sin_θ*cos_ø, cos_θ*cos_ø, -sin_ø,},
			{sin_θ*sin_ø, cos_θ*sin_ø,  cos_ø,},
			{      cos_θ,      -sin_θ,   0.0l,},
		};
		mat3x3q_t xfrm(m);
		const vec3q_t q = xfrm.mul(p);
		return cartesian3dvec {q.raw[0], q.raw[1], q.raw[2]};
	}
};

template <typename VecT, typename T=long double>
struct funmat {};

template <typename VecT, typename T=long double>
struct funmat3x3 : public funmat<VecT, T>
{
	typedef T TT;
	typedef VecT VecTT;
	typedef mat3x3q_t XfrmMatrixT;
	typedef std::function<T(const VecTT&, int, int)> coeffun;

	static constexpr coeffun ident = [](const VecTT &v, int r, int c) -> T { return r == c ? v.raw[r] : T(0.0); };
	static constexpr coeffun zero = [](const VecT &v, int r, int c) -> T { return T(0.0); };
	static constexpr coeffun one = [](const VecT &v, int r, int c) -> T { return T(1.0); };
	static constexpr coeffun sine = [](const VecT &v, int r, int c) -> T { return r == c ? sinl(v.raw[r]) : T(0.0); };
	static constexpr coeffun cosine = [](const VecT &v, int r, int c) -> T { return r == c ? cosl(v.raw[r]) : T(0.0); };

	static constexpr coeffun sinθ = [](const spherical3dvec&v, int, int) -> TT { return sinl(v.θ()); };
	static constexpr coeffun sinø = [](const spherical3dvec&v, int, int) -> TT { return sinl(v.ø()); };
	static constexpr coeffun cosθ = [](const spherical3dvec&v, int, int) -> TT { return cosl(v.θ()); };
	static constexpr coeffun cosø = [](const spherical3dvec&v, int, int) -> TT { return cosl(v.ø()); };
	static constexpr coeffun sinθcosø = [](const spherical3dvec&v, int, int) -> TT { return sinl(v.θ()) * cosl(v.ø()); };
	static constexpr coeffun cosθcosø = [](const spherical3dvec&v, int, int) -> TT { return cosl(v.θ()) * cosl(v.ø()); };
	static constexpr coeffun sinθsinø = [](const spherical3dvec&v, int, int) -> TT { return sinl(v.θ()) * sinl(v.ø()); };
	static constexpr coeffun cosθsinø = [](const spherical3dvec&v, int, int) -> TT { return cosl(v.θ()) * sinl(v.ø()); };

	static constexpr XfrmMatrixT identitymatrix = {
		{  one, zero, zero },
		{ zero,  one, zero },
		{ zero, zero,  one },
	};

	static constexpr XfrmMatrixT R_1_mtrx ={
		{  one,  zero, zero },
		{ zero,  cosθ, sinθ },
		{ zero, -sinθ, cosθ },
	};

	static constexpr XfrmMatrixT R_2_mtrx = {
		{ cosθ, zero, -sinθ },
		{ zero,  one,  zero },
		{ sinθ, zero,  cosθ },
	};

	static constexpr XfrmMatrixT R_3_mtrx = {
		{  cosθ, sinθ, zero },
		{ -sinθ, cosθ, zero },
		{  zero, zero,  one },
	};

	static constexpr XfrmMatrixT sph2cart_mtrx = {
		{ sinθcosø, cosθcosø, -sinø },
		{ sinθsinø, cosθsinø,  cosø },
		{     cosθ,    -sinθ,  zero },
	};

	XfrmMatrixT xfrmmtrx = identitymatrix;

	void setxfrmmtrx(const XfrmMatrixT &mtrx) {
		xfrmmtrx = mtrx;
	}

	VecT mul(const VecT &v) const {
		return {
			xfrmmtrx.elems[0][0](v, 0, 0) * v.raw[0] + xfrmmtrx.elems[0][1](v, 0, 1) * v.raw[1] + xfrmmtrx.elems[0][2](v, 0, 2) * v.raw[2],
			xfrmmtrx.elems[1][0](v, 1, 0) * v.raw[0] + xfrmmtrx.elems[1][1](v, 1, 1) * v.raw[1] + xfrmmtrx.elems[1][2](v, 1, 2) * v.raw[2],
			xfrmmtrx.elems[2][0](v, 2, 0) * v.raw[0] + xfrmmtrx.elems[2][1](v, 2, 1) * v.raw[1] + xfrmmtrx.elems[2][2](v, 2, 2) * v.raw[2],
		};
	}

	XfrmMatrixT mul(const XfrmMatrixT &m) const {
		XfrmMatrixT product;
		//product[0][0] = [](const VecT &v, int r, int c) -> T { return xfrmmtrx.elems[r][c](v,r,c) * m.elems[0][0](v,r,c); };

		return {
			//[](const VecT &v, int r, int c) -> T { return xfrmmtrx.elems[0][0](v,r,c) * m.elems[0][0](v,r,c); };
		};
	}
};

#endif /* _PAULYC_LALGEBRA_HPP_ */
