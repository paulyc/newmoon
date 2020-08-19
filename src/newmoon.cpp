/**
 * newmoon.cpp - part of newmoon, moon phase calculator
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

#include "lalgebra.hpp"
#include "calculus.hpp"
#include "jd_clock.hpp"
#include "ephemshelper.hpp"

#include <iostream>
#include <functional>

typedef std::function<__float128(JPLEphems&, const jd_clock::time_point&)> f_type;

f_type moonSunAngle = [](JPLEphems &ephems, const jd_clock::time_point &jd) -> __float128 {
    const double jd_now = static_cast<double>(jd_clock::duration(jd.time_since_epoch()).count());
    cartesian3dvec moonpos = ephems.get_state(jd_now, JPLEphems::EarthMoonBarycenter, JPLEphems::Moon).position();
    const __float128 moonR = moonpos.mag();
    cartesian3dvec sunpos = ephems.get_state(jd_now, JPLEphems::Earth, JPLEphems::Sun).position();
    const __float128 sunR = sunpos.mag();
    JPLEphems::NutationState ns = ephems.get_nutations(jd_now);
    const __float128 α_moon = atanq(moonpos.y()/moonpos.x());// asinq(moonpos.y() / (moonR*cosq(δ_moon)));
    const __float128 δ_moon = atanq(moonpos.z()/(moonpos.y()*sinq(α_moon)));

    const __float128 α_sun = atanq(sunpos.y()/sunpos.x());//asinq(sunpos.y() / (moonR*cosq(δ_sun)));
    const __float128 δ_sun = atanq(sunpos.z()/(sunpos.y()*sinq(α_sun)));

    const __float128 Δψ = ns.nutationInLongitude();
    const __float128 Δɛ = ns.nutationInObliquity();
    const __float128 arg = moonpos.angle(sunpos);
#if 1
    return arg;
#else
    const __float128 ɛ_0 = 0.40904635907; //23.43663 deg in radians
    const __float128 ɛ = ɛ_0;// + Δɛ;

    const __float128 Δα_moon = (cosq(ɛ)  + sinq(ɛ)*sinq(α_moon)*tanq(δ_moon)) * Δψ - cosq(α_moon)*tan(δ_moon)*Δɛ;
    const __float128 Δδ_moon = cosq(α_moon)*sinq(ɛ)*Δψ + sinq(α_moon)*Δɛ;
    const __float128 Δα_sun = (cosq(ɛ)  + sinq(ɛ)*sinq(α_sun)*tanq(δ_sun)) * Δψ - cosq(α_sun)*tan(δ_sun)*Δɛ;
    const __float128 Δδ_sun = cosq(α_sun)*sinq(ɛ)*Δψ + sinq(α_sun)*Δɛ;
    const __float128 arg = (α_moon - Δα_moon) - (α_sun - Δα_sun);
    JPLEphems::State ls = ephems.get_librations(jd_now);
#endif
};

auto minFinder = [](JPLEphems &ephems, std::chrono::system_clock::time_point &t, jd_clock::time_point &jd, f_type moonSunAngle) -> std::chrono::system_clock::time_point {
    std::chrono::system_clock::time_point mintp = std::chrono::system_clock::now();
    __float128 minangle = MMM_2_PI;
    for (int i = 0; i < 29*24*60; ++i) {
        jd += jd_clock::duration(60.0l/jd_clock::SECONDS_PER_JDAY);
        t += std::chrono::seconds(60);

        const __float128 arg = moonSunAngle(ephems, jd);
        const __float128 argabs = fabsq(arg);
        if (argabs < minangle) {
            minangle = argabs;
            mintp = t;
        } else if (argabs > MMM_PI/4.0q && minangle < MMM_PI/8.0q) {
            break;
        }
    }
    std::cout << "minangle (newmoon) " << static_cast<long double>(minangle) << " at mintp " << mintp << std::endl;
    return mintp;
};

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;
    
    JPLEphems ephems;

    char tzbuf[] = "TZ=UTC";
    putenv(tzbuf);
    tzset();

    std::chrono::system_clock::time_point t = std::chrono::system_clock::now();
    t -= std::chrono::seconds(15*86400);
    jd_clock::time_point jd = jd_clock::now();
    jd -= jd_clock::duration(15.0);

    //std::cout << "hello newmoon " << t << " JD = " << jd << std::endl;

    timespec ts = {0,100000000}; //5.000000000s
    try {
        ephems.init("ephem/lnxm13000p17000.431");
    } catch (const std::exception &ex0) {
        try {
            ephems.init("ephem/linux_p1550p2650.430t");
        } catch (const std::exception &ex1) {
            try {
                ephems.init("ephem/linux_p1550p2650.430");
            } catch (const std::exception &ex2) {
                std::cerr << "Couldn't initialize ephems: " <<
                             ex0.what() << ' ' <<
                             ex1.what() << ' ' <<
                             ex2.what() << std::endl;
                return 1;
            }
        }
    }

    // generate newmoons every ts seconds forever
    for (;;) {
        std::chrono::system_clock::time_point tp = minFinder(ephems, t, jd, moonSunAngle);
        std::cout << "\"new moon (ISO8601)\"" << '"' << tp << '"' << std::endl;
        nanosleep(&ts, nullptr);
    }

    return 0;
}
