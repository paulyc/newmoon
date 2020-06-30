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
#include "jd_clock.hpp"
#include "ephemshelper.hpp"

#include <iostream>

int run()
{
    JPLEphems ephems;

    char tzbuf[] = "TZ=UTC";
    putenv(tzbuf);
    tzset();

    std::chrono::system_clock::time_point t = std::chrono::system_clock::now();
    jd_clock::time_point jd = jd_clock::now();

    std::cout << "hello newmoon " << t << " JD = " << jd << std::endl;

    timespec ts = {5,0}; //5.000000000s
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

	auto newway = [&ephems, &t, &jd]() {
		std::chrono::system_clock::time_point mintp = std::chrono::system_clock::now();
        __float128 minangle = MMM_2_PI;
		for (int i = 0; i < 29*24*60; ++i) {
			jd += jd_clock::duration(60.0l/jd_clock::SECONDS_PER_JDAY);
			t += std::chrono::seconds(60);
			const double jd_now = static_cast<double>(jd_clock::duration(jd.time_since_epoch()).count());
            cartesian3dvec moonpos = ephems.get_state(jd_now, JPLEphems::EarthMoonBarycenter, JPLEphems::Moon).position();
            const __float128 moonR = moonpos.mag();
            cartesian3dvec sunpos = ephems.get_state(jd_now, JPLEphems::Earth, JPLEphems::Sun).position();
            const __float128 sunR = sunpos.mag();
            JPLEphems::NutationState ns = ephems.get_nutations(jd_now);
            //asinq(moonpos.z() / moonR);
            const __float128 α_moon = atanq(moonpos.y()/moonpos.x());// asinq(moonpos.y() / (moonR*cosq(δ_moon)));
            const __float128 δ_moon = atanq(moonpos.z()/(moonpos.y()*sinq(α_moon)));

            const __float128 α_sun = atanq(sunpos.y()/sunpos.x());//asinq(sunpos.y() / (moonR*cosq(δ_sun)));
            const __float128 δ_sun = atanq(sunpos.z()/(sunpos.y()*sinq(α_sun)));

            const __float128 Δψ = ns.nutationInLongitude();
            const __float128 Δɛ = ns.nutationInObliquity();
#define MATH_IS_HARD 1
#if !MATH_IS_HARD
            const __float128 ɛ_0 = 0.40904635907; //23.43663 deg in radians
            const __float128 ɛ = ɛ_0;// + Δɛ;

            const __float128 Δα_moon = (cosq(ɛ)  + sinq(ɛ)*sinq(α_moon)*tanq(δ_moon)) * Δψ - cosq(α_moon)*tan(δ_moon)*Δɛ;
            const __float128 Δδ_moon = cosq(α_moon)*sinq(ɛ)*Δψ + sinq(α_moon)*Δɛ;
            const __float128 Δα_sun = (cosq(ɛ)  + sinq(ɛ)*sinq(α_sun)*tanq(δ_sun)) * Δψ - cosq(α_sun)*tan(δ_sun)*Δɛ;
            const __float128 Δδ_sun = cosq(α_sun)*sinq(ɛ)*Δψ + sinq(α_sun)*Δɛ;
            const __float128 arg = (α_moon - Δα_moon) - (α_sun - Δα_sun);
            //if (fabsq(arg) < 0.0000000001) {
            //    break;
            //}
#else

           // yeah i have no clue but one of these angles might have something to do with phase
#define USE_THE_LIBRATION 0
#if USE_THE_LIBRATION
            JPLEphems::State ls = ephems.get_librations(jd_now);

#endif
            //works better than using RAs and nutations anyway! go figure
            const __float128 arg = moonpos.angle(sunpos);
#endif
            const __float128 argabs = fabsq(arg);
			if (argabs < minangle) {
				minangle = argabs;
				mintp = t;
            } else if (argabs > MMM_PI/4.0q && minangle < 0.0001) {
				break;
            }
		}
        std::cout << "minangle (newmoon) " << static_cast<long double>(minangle) << " at mintp " << mintp << std::endl;
	};

	auto oldway = [&ephems, &t, &jd]() {
		std::chrono::system_clock::time_point mintp = std::chrono::system_clock::now();
        __float128 minangle = MMM_2_PI;
		// just find the time when its the minimum
        __float128 lastmags[2] = {-3.0q};
        size_t lastmag_indx = 0;
		__float128 minmag = 3.0q;
        for (int i = 0; i < 29*24*60; ++i) {
            jd += jd_clock::duration(60.0l/jd_clock::SECONDS_PER_JDAY);
            t += std::chrono::seconds(60);
            const double jd_now = static_cast<double>(jd_clock::duration(jd.time_since_epoch()).count());
            const cartesian3dvec sm = ephems.get_state(jd_now, JPLEphems::Earth, JPLEphems::Moon).position();
            const cartesian3dvec ss = ephems.get_state(jd_now, JPLEphems::Sun, JPLEphems::EarthMoonBarycenter).position();
            const cartesian3dvec sum = sm.sum(ss);
            const __float128 mag = sum.mag();
            const __float128 lastmag = lastmags[lastmag_indx % 2];
            ++lastmag_indx;
            lastmags[lastmag_indx % 2] = mag;
            const __float128 dmag = mag - lastmag;
            //std::cout << "Magnitude " << magphase[0] << " Angle " << magphase[1] << " at date " << t << std::endl;
            if (mag < minmag && dmag < 0.0q) {
                minmag = mag;
                mintp = t;
            } else if (mag > 1.0q && minmag < 0.05q && dmag > 0.0q) {
                break;
            }
        }
        std::cout << "minmag (newmoon) " << static_cast<long double>(minmag) << " at mintp " << mintp << std::endl;
	};

    // generate newmoons every ts seconds forever
    for (;;) {
#if 1
		newway();
#else
		oldway();
#endif

        nanosleep(&ts, nullptr);
    }

    std::cout << "goodbye newmoon" << std::endl;
    return 0;
}

int test() {
    jd_clock::test_delta_t_lerp();
    auto jd = jd_clock::now();
    std::cerr << jd << std::endl;
    std::cerr << jd.time_since_epoch().count() << std::endl;
    const double jd_now = static_cast<double>(jd_clock::duration(jd.time_since_epoch()).count());
    std::cerr << jd_now << std::endl;
    return 0;
}

int main() {
    return test() || run();
}
