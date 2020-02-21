/*
* ephemshelper.hpp - part of newmoon, moon phase calculator
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

#ifndef _PAULYC_EPHEMSHELPER_H_
#define _PAULYC_EPHEMSHELPER_H_

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
        union {
            double pv[6];
            vec3d_t position;
            vec3d_t velocity;
        };

        long double angle(const State &that) const {
			return position.drop_z().normalize().angle(that.position.drop_z().normalize());

        }
        vec2q_t magphase(const State &other) const {
			return position.drop_z().normalize().sum(other.position.drop_z().normalize()).magphase();
        }
		vec2q_t ra_dec() const {
			const long double distance = position.mag();
			long double ra = position.phase();
			if (ra < 0.0l) {
				ra += 2.0l * MMM_PI;
			}
			return {ra, asinl(position.raw[1]/distance)};
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
            throw std::runtime_error(string_format("jpl_init_ephemeris returned code %d", jpl_init_error_code()));
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
            throw std::runtime_error(string_format("jpl_pleph returned code %d", res));
        }
        return result;
    }
	State get_nutations(double jdt)
	{
		State result;
		if (!initialized()) {
			throw std::runtime_error("try calling JPLEphems::init() first");
		}
		int res = jpl_pleph(_ephdata, jdt, Nutations, Earth, result.pv, 0);
		if (res != 0) {
            throw std::runtime_error(string_format("jpl_pleph returned code %d", res));
        }
        return result;
	}
private:
    char _names[MAX_CONSTANTS][6];
    double _values[MAX_CONSTANTS];
    jpl_eph_data *_ephdata;
};

#endif /* _PAULYC_EPHEMSHELPER_H_ */
