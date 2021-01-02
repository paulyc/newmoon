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

long double moonSunAngle(JPLEphems &ephems, const jd_clock::time_point &jd) {
    const double jd_now = static_cast<double>(jd_clock::duration(jd.time_since_epoch()).count());
    cartesian3dvec moonpos = ephems.get_state(jd_now, JPLEphems::EarthMoonBarycenter, JPLEphems::Moon).position();
    const long double moonR = moonpos.mag();
    cartesian3dvec sunpos = ephems.get_state(jd_now, JPLEphems::Earth, JPLEphems::Sun).position();
    const long double sunR = sunpos.mag();
    JPLEphems::NutationState ns = ephems.get_nutations(jd_now);
    const long double α_moon = atanq(moonpos.y()/moonpos.x());// asinq(moonpos.y() / (moonR*cosq(δ_moon)));
    const long double δ_moon = atanq(moonpos.z()/(moonpos.y()*sinq(α_moon)));

    const long double α_sun = atanq(sunpos.y()/sunpos.x());//asinq(sunpos.y() / (moonR*cosq(δ_sun)));
    const long double δ_sun = atanq(sunpos.z()/(sunpos.y()*sinq(α_sun)));

    const long double Δψ = ns.nutationInLongitude();
    const long double Δɛ = ns.nutationInObliquity();
    const long double arg = moonpos.angle(sunpos);
#if 1
    return arg;
#else
    const long double ɛ_0 = 0.40904635907; //23.43663 deg in radians
    const long double ɛ = ɛ_0;// + Δɛ;

    const long double Δα_moon = (cosq(ɛ)  + sinq(ɛ)*sinq(α_moon)*tanq(δ_moon)) * Δψ - cosq(α_moon)*tan(δ_moon)*Δɛ;
    const long double Δδ_moon = cosq(α_moon)*sinq(ɛ)*Δψ + sinq(α_moon)*Δɛ;
    const long double Δα_sun = (cosq(ɛ)  + sinq(ɛ)*sinq(α_sun)*tanq(δ_sun)) * Δψ - cosq(α_sun)*tan(δ_sun)*Δɛ;
    const long double Δδ_sun = cosq(α_sun)*sinq(ɛ)*Δψ + sinq(α_sun)*Δɛ;
    const long double arg = (α_moon - Δα_moon) - (α_sun - Δα_sun);
    JPLEphems::State ls = ephems.get_librations(jd_now);
#endif
};

std::chrono::system_clock::time_point minFinder(JPLEphems &ephems, jd_clock::time_point &jd) {
    std::chrono::system_clock::time_point mintp = std::chrono::system_clock::now();
    long double minangle = MMM_2_PI;
    for (int i = 0; i < 29*24*60; ++i) {
        jd += jd_clock::duration(60.0l/jd_clock::SECONDS_PER_JDAY);

        const long double arg = moonSunAngle(ephems, jd);
        const long double argabs = fabsq(arg);
        if (argabs < minangle) {
            minangle = argabs;
            mintp = jd_clock::to_system_clock(jd);
        } else if (argabs > MMM_PI/4.0q && minangle < MMM_PI/8.0q) {
            break;
        }
    }
    std::cout << "minangle (newmoon) " << minangle << " at mintp " << mintp << std::endl;
    return mintp;
};

#if CALCULUS_WORKS
auto minFinder2 = [](JPLEphems &ephems, jd_clock::time_point &jd, f_type moonSunAngle) -> long double {
    double t_begin = jd.time_since_epoch().count();
    double t_end = t_begin + 29 * 24 * 60 * 60;
    return github::paulyc::min_x([&ephems, &moonSunAngle](long double t) -> long double {return moonSunAngle(ephems, jd_clock::from_unixtime(t));}, t_begin, t_end, 1.0).value();
};
#endif