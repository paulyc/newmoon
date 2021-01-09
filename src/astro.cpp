/**
 * astro.cpp
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

#include "astro.hpp"

typedef std::function<long double(JPLEphems&, const jd_clock::time_point&)> f_type;

struct result {
    double jd_now;
    cartesian3dvec moonpos;
    long double moonR;
    cartesian3dvec sunpos;
    long double sunR ;
    JPLEphems::NutationState ns;
    long double α_moon;
    long double δ_moon;
    long double α_sun;
    long double δ_sun;
    long double dα;
    long double Δψ;
    long double Δɛ;
    long double arg;
};
std::ostream &operator<<(std::ostream &os, const result &res) {
    os << "arg " << res.arg << " α_moon " << res.α_moon << " α_sun " << res.α_sun << " dα " << res.dα;
    return os;
}

void moonSunAngle(JPLEphems &ephems, const jd_clock::time_point &jd, result &res) {
    res.jd_now = static_cast<double>(jd_clock::duration(jd.time_since_epoch()).count());
    res.moonpos = ephems.get_state(res.jd_now, JPLEphems::EarthMoonBarycenter, JPLEphems::Moon).position();
    res.moonR = res.moonpos.mag();
    res.sunpos = ephems.get_state(res.jd_now, JPLEphems::Earth, JPLEphems::Sun).position();
    res.sunR = res.sunpos.mag();
    res.ns = ephems.get_nutations(res.jd_now);
    res.α_moon = atanq(res.moonpos.y()/res.moonpos.x());// asinq(moonpos.y() / (moonR*cosq(δ_moon)));
    res.δ_moon = atanq(res.moonpos.z()/(res.moonpos.y()*sinq(res.α_moon)));

    res.α_sun = atanq(res.sunpos.y()/res.sunpos.x());//asinq(sunpos.y() / (moonR*cosq(δ_sun)));
    res.δ_sun = atanq(res.sunpos.z()/(res.sunpos.y()*sinq(res.α_sun)));

    res.dα = res.α_moon - res.α_sun;

    res.Δψ = res.ns.nutationInLongitude();
    res.Δɛ = res.ns.nutationInObliquity();
    res.arg = res.moonpos.angle(res.sunpos);

#if 0
    const long double ɛ_0 = 0.40904635907; //23.43663 deg in radians
    const long double ɛ = ɛ_0;// + Δɛ;

    const long double Δα_moon = (cosq(ɛ)  + sinq(ɛ)*sinq(α_moon)*tanq(δ_moon)) * Δψ - cosq(α_moon)*tan(δ_moon)*Δɛ;
    const long double Δδ_moon = cosq(α_moon)*sinq(ɛ)*Δψ + sinq(α_moon)*Δɛ;
    const long double Δα_sun = (cosq(ɛ)  + sinq(ɛ)*sinq(α_sun)*tanq(δ_sun)) * Δψ - cosq(α_sun)*tan(δ_sun)*Δɛ;
    const long double Δδ_sun = cosq(α_sun)*sinq(ɛ)*Δψ + sinq(α_sun)*Δɛ;
    const long double arg = (α_moon - Δα_moon) - (α_sun - Δα_sun);
    JPLEphems::State ls = ephems.get_librations(jd_now);
#endif
}

inline int signum(long double x) {
    if (x >= 0.0q) {
        return 1;
    } else {
        return -1;
    }
}
std::chrono::system_clock::time_point minFinder(JPLEphems &ephems, jd_clock::time_point &jd) {
    std::chrono::system_clock::time_point mintp = std::chrono::system_clock::now();
    long double minangle = MMM_2_PI;
    result res;
    moonSunAngle(ephems, jd, res);
    int sign_dα = signum(res.dα);
    for (int i = 0; i < 29*24*60; ++i) {
        jd += jd_clock::duration(60.0l/jd_clock::SECONDS_PER_JDAY);

        moonSunAngle(ephems, jd, res);
        int sign = signum(res.dα);
        if (sign_dα != sign) {
            //std::cout << res << std::endl;
            if (sign_dα == 1) {
                sign_dα = -1;
            } else {
                mintp = jd_clock::to_system_clock(jd);
                //std::cout << "minangle (newmoon) " << minangle << " at mintp " << mintp << std::endl;
                return mintp;
            }
        }
        /*const long double argabs = fabsq(res.arg);
        if (argabs < minangle) {
            minangle = argabs;
            mintp = jd_clock::to_system_clock(jd);
        } else if (argabs > MMM_PI/4.0q && minangle < MMM_PI/8.0q) {
            break;
        }*/
    }
    std::cout << " none found " << mintp << std::endl;
    return mintp;
};

#if CALCULUS_WORKS
auto minFinder2 = [](JPLEphems &ephems, jd_clock::time_point &jd, f_type moonSunAngle) -> long double {
    double t_begin = jd.time_since_epoch().count();
    double t_end = t_begin + 29 * 24 * 60 * 60;
    return github::paulyc::min_x([&ephems, &moonSunAngle](long double t) -> long double {return moonSunAngle(ephems, jd_clock::from_unixtime(t));}, t_begin, t_end, 1.0).value();
};
#endif
