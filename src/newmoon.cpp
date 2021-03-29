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

#include "tetrabiblos.hpp"

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;
    
    JPLEphems ephems;

    char tzbuf[] = "TZ=UTC";
    putenv(tzbuf);
    tzset();

    timespec ts = {0,100000000};
    try {
        ephems.init("ephem/lnxm13000p17000.431");
    } catch (const std::exception &ex) {
        std::cerr << "Couldn't initialize ephems: " << ex.what() << std::endl;
        return 1;
    }


    github::paulyc::tetrabiblos::Date today = github::paulyc::tetrabiblos::getDate(ephems, std::chrono::system_clock::now());
    std::cout << today << std::endl;

    jd_clock::time_point jd = jd_clock::now();
    jd -= jd_clock::duration(28.0);

    // generate newmoons every ts seconds forever
    for (;;) {
        std::chrono::system_clock::time_point tp = minFinder(ephems, jd);
        std::cout << "\"new moon (ISO8601)\"" << '"' << tp << '"' << std::endl;
        jd += jd_clock::duration(1.0);
        nanosleep(&ts, nullptr);
    }

    return 0;
}
