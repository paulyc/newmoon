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

#include "newmoon.hpp"
#include "ephemshelper.hpp"

#include <iostream>

namespace paulyc {

//NewMoon::NewMoon(NewMoon &&other) {
//    _ephems.swap(other._ephems);
//}

NewMoon& NewMoon::operator=(NewMoon &&rhs) {
    _ephems.swap(rhs._ephems);
    return *this;
}

void NewMoon::init()
{
    char tzbuf[] = "TZ=UTC";
    putenv(tzbuf);
    tzset();

    std::chrono::system_clock::time_point _t = std::chrono::system_clock::now();
    jd_clock::time_point _jd = jd_clock::now();

    std::cout << "hello newmoon current time = " << _t << " JD = " << _jd << std::endl;

    _ephems = std::make_unique<JPLEphems>();

    try {
        _ephems->init("ephem/lnxm13000p17000.431");
    } catch (const std::exception &ex0) {
        try {
            _ephems->init("ephem/linux_p1550p2650.430t");
        } catch (const std::exception &ex1) {
            try {
                _ephems->init("ephem/linux_p1550p2650.430");
            } catch (const std::exception &ex2) {
                std::cerr << "Couldn't initialize ephems: " <<
                             ex0.what() << ' ' <<
                             ex1.what() << ' ' <<
                             ex2.what() << std::endl;
            }
        }
    }
}

std::chrono::system_clock::time_point NewMoon::next_moon()
{
    //static timespec ts = {5,0}; //5.000000000s
    std::chrono::system_clock::time_point mintp = std::chrono::system_clock::now();
    __float128 minangle = MMM_2_PI;
    for (int i = 0; i < 29*24*60; ++i) {
        _jd += jd_clock::duration(60.0l/jd_clock::SECONDS_PER_JDAY);
        _t += std::chrono::seconds(60);
        const double jd_now = static_cast<double>(jd_clock::duration(_jd.time_since_epoch()).count());
        cartesian3dvec moonpos = _ephems->get_state(jd_now, JPLEphems::Earth, JPLEphems::Moon).position();
        cartesian3dvec sunpos = _ephems->get_state(jd_now, JPLEphems::EarthMoonBarycenter, JPLEphems::Sun).position();

        const __float128 arg = moonpos.angle(sunpos);
        const __float128 argabs = fabsl(arg);
        if (argabs < minangle) {
            minangle = argabs;
            mintp = _t;
        } else if (argabs > MMM_PI/4.0q && minangle < MMM_PI/16.0q) {
            break;
        }
    }
    std::cout << "minangle (newmoon) " << static_cast<long double>(minangle) << " at mintp " << mintp << std::endl;
    return mintp;
}

}

int main() {
    paulyc::NewMoon nm;
    nm.init();
    for (;;) {
        auto next = nm.next_moon();
        std::cout << next << std::endl;
    }
    return 0;
}
