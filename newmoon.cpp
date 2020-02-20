/*
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
 */

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
        long double minangle = 2*3.14159;
		for (int i = 0; i < 29*24*60; ++i) {
			jd += jd_clock::duration(60.0l/jd_clock::SECONDS_PER_JDAY);
			t += std::chrono::seconds(60);
			const double jd_now = static_cast<double>(jd_clock::duration(jd.time_since_epoch()).count());
			JPLEphems::State sm = ephems.get_state(jd_now, JPLEphems::Earth, JPLEphems::Moon);
			JPLEphems::State ss = ephems.get_state(jd_now, JPLEphems::Earth, JPLEphems::Sun);
			const long double arg = sm.angle(ss);
			if (abs(arg) < minangle) {
				minangle = abs(arg);
				mintp = t;
			} else if (abs(arg) > 3.14159/4 && minangle < 3.14159/8) {
				break;
			}
		}
		std::cout << "minangle (newmoon) " << minangle << " at mintp " << mintp << std::endl;
	};

	auto oldway = [&ephems, &t, &jd]() {
		std::chrono::system_clock::time_point mintp = std::chrono::system_clock::now();
        long double minangle = 2*3.14159;
		// just find the time when its the minimum
        long double lastmags[2] = {-3.0l};
        size_t lastmag_indx = 0;
		long double minmag = 3.0l;
        for (int i = 0; i < 29*24*60; ++i) {
            jd += jd_clock::duration(60.0l/jd_clock::SECONDS_PER_JDAY);
            t += std::chrono::seconds(60);
            const double jd_now = static_cast<double>(jd_clock::duration(jd.time_since_epoch()).count());
            JPLEphems::State sm = ephems.get_state(jd_now, JPLEphems::Earth, JPLEphems::Moon);
            JPLEphems::State ss = ephems.get_state(jd_now, JPLEphems::Sun, JPLEphems::EarthMoonBarycenter);
            vec2q_t magphase = sm.ra_magphase(ss);
            const long double mag = magphase.x;
            const long double lastmag = lastmags[lastmag_indx % 2];
            ++lastmag_indx;
            lastmags[lastmag_indx % 2] = mag;
            const long double dmag = mag - lastmag;
            //std::cout << "Magnitude " << magphase[0] << " Angle " << magphase[1] << " at date " << t << std::endl;
            if (mag < minmag && dmag < 0.0l) {
                minmag = mag;
                mintp = t;
            } else if (mag > 1.0l && minmag < 0.05l && dmag > 0.0l) {
                break;
            }
        }
        std::cout << "minmag (newmoon) " << minmag << " at mintp " << mintp << std::endl;
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
