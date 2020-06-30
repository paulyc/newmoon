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

#ifndef PAULYC_EPHEMSHELPER_HPP
#define PAULYC_EPHEMSHELPER_HPP

#include "jpl_int.h"
#include "jpleph.h"
#include "lalgebra.hpp"
#include "jd_clock.hpp"

class JPLEphems
{
    static constexpr size_t MAX_CONSTANTS = 1024;
public:
    enum Point {
        Mercury =  1,
        Venus   =  2,
        Earth   =  3,
        Mars    =  4,
        Jupiter =  5,
        Saturn  =  6,
        Uranus  =  7,
        Neptune =  8,
        Pluto   =  9,
        Moon    = 10,
        Sun     = 11,
        SolarSystemBarycenter = 12,
        EarthMoonBarycenter   = 13,
        Nutations             = 14,
        Librations            = 15,
        LunarMantleOmega      = 16,
        TT_TDB                = 17,
    };

    struct State
    {
        double pv[6];
        cartesian3dvec position() const {
            return {{static_cast<__float128>(pv[0]), static_cast<__float128>(pv[1]), static_cast<__float128>(pv[2])}};
        }
        cartesian3dvec velocity() const {
            return {{static_cast<__float128>(pv[3]), static_cast<__float128>(pv[4]), static_cast<__float128>(pv[5])}};
        }
    };
    struct NutationState
    {
        double pv[6];
        //radians - delta psi
        __float128 nutationInLongitude() const {
            return static_cast<__float128>(pv[0]);
        }
        //radians - delta epsilon
        __float128 nutationInObliquity() const {
            return static_cast<__float128>(pv[1]);
        }
        //radians/day
        __float128 nutationInLongitudeRate() const {
            return static_cast<__float128>(pv[2]);
        }
        //radians/day
        __float128 nutationInObliquityRate() const {
            return static_cast<__float128>(pv[3]);
        }
    };
    JPLEphems() : _ephdata(nullptr) {}
    ~JPLEphems()
    {
        if (_ephdata != nullptr) {
            jpl_close_ephemeris(_ephdata);
            _ephdata = nullptr;
        }
    }

    void init(const std::string &filename)
    {
        _ephdata = static_cast<jpl_eph_data*>(jpl_init_ephemeris(filename.c_str(), _names, _values));
        if (_ephdata == nullptr) {
            throw std::runtime_error("jpl_init_ephemeris returned code %d"_fmt.format(jpl_init_error_code()));
        }
    }
    bool initialized() const { return _ephdata != nullptr; }
    State get_state(double jdt, Point center, Point ref)
    {
        State result;
		if (!initialized()) {
			throw std::runtime_error("try calling JPLEphems::init() first");
		}
        int res = jpl_pleph(_ephdata, jdt, ref, center, result.pv, 0);
        if (res != 0) {
            throw std::runtime_error("jpl_pleph returned code %d"_fmt.format(res));
        }
        return result;
    }
/*
    OBLIQUITY OF THE ECLIPTIC, NUTATION AND LATITUDES
    OF THE ARCTIC AND ANTARCTIC CIRCLES

    NUTATION COMPUTED USING THE IAU 2000B SERIES

    CALENDAR DATE = 2020 Febuary 21 Friday (Gregorian)
    TIME      HMS =  06:30:05 = +0.27089120 day
    ± Delta T HMS = +00:00:00 = +0.00000000 day

    JD NUMBER FOR DATE, TIME AND GIVEN DELTA T
    JD = 2458900.77089120 TT

    TIME VARIABLE CORRESPONDING TO JD
    T = +0.2013900312 = Julian centuries reckoned from J2000.0
    t = +0.0201390031 = Julian millennia reckoned from J2000.0

    ---------------------------------------------------------------------------
    OBLIQUITY OF THE ECLIPTIC

    Eps Mean = 23.4366725232°  = 23° 26' 12.021" (Laskar)
    Nutation = -0.0001406109°  =         -0.506" (IAU 2000B)
    Eps True = 23.4365319123°  = 23° 26' 11.515"

    ---------------------------------------------------------------------------
    LATITUDES OF ARCTIC(+) AND ANTARCTIC(−) CIRCLES

    ± 66° 33' 48.485"  =  ± 66.5634680877°

    static __float128 obliquity_of_ecliptic(double jd2000)
    {
        return 23.4393 - 3.563E-7 * jd2000;
    }
*/
    NutationState get_nutations(double jdt)
	{
        NutationState result;
        if (!initialized()) {
            throw std::runtime_error("try calling JPLEphems::init() first"_fmt.format());
		}

        int res = jpl_pleph(_ephdata, jdt, Nutations, 0, result.pv, 0);
        if (res != 0) {
            throw std::runtime_error("jpl_pleph returned code %d"_fmt.format(res));
        }
        return result;
	}

    State get_librations(double jdt)
    {
        State result;
        if (!initialized()) {
            throw std::runtime_error("try calling JPLEphems::init() first"_fmt.format());
        }
        int res = jpl_pleph(_ephdata, jdt, Librations, 0, result.pv, 0);
        if (res != 0) {
            throw std::runtime_error("jpl_pleph returned code %d"_fmt.format(res));
        }
        return result;
    }
private:
    char _names[MAX_CONSTANTS][6];
    double _values[MAX_CONSTANTS];
    jpl_eph_data *_ephdata;
};

#endif /* PAULYC_EPHEMSHELPER_HPP */
