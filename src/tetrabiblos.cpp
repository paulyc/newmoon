/**
 * tetrabiblos.cpp
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
#include <chrono>
#include <ctime>

namespace github {
namespace paulyc {
namespace tetrabiblos {

using std::chrono::system_clock;

static int monthLength(std::tm *t) {
    switch (t->tm_mon) {
    case 8:
    case 3:
    case 7:
    case 10:
        return 30;
    case 1:
        return ((t->tm_year % 4) == 0 && (t->tm_year % 400) != 0) ? 29 : 28;
    default:
        return 31;
    }
}

static void fixTm(std::tm *t) {
    int mlen = monthLength(t);
    if (t->tm_mday > mlen) {
        t->tm_mday -= mlen;
        ++t->tm_mon;
        if (t->tm_mon > 11) {
            t->tm_mon = 0;
            ++t->tm_year;
        }
    }
}

Date getDate(JPLEphems &ephems, const system_clock::time_point &tp) {
    Date d = {false, 0, 0, NONE, 0};
    std::time_t tp_time = system_clock::to_time_t(tp);
    tm *tp_tm = gmtime(&tp_time);
    jd_clock::time_point jd = jd_clock::from_system_clock(tp);
    // find next new moon
    const system_clock::time_point nextNewMoon = minFinder(ephems, jd);
    const std::time_t nextNewMoonTs = system_clock::to_time_t(nextNewMoon);
    jd -= jd_clock::duration(30.0);
    const system_clock::time_point lastNewMoon = minFinder(ephems, jd);
    const std::time_t lastNewMoonTs = system_clock::to_time_t(lastNewMoon);
    tm * lastTm = gmtime(&lastNewMoonTs);
    tm * nextTm = gmtime(&nextNewMoonTs);
    if (tp_tm->tm_mon >= nextTm->tm_mon || (tp_tm->tm_mon == nextTm->tm_mon && tp_tm->tm_mday >= nextTm->tm_mday)) {
        // TODO
    }
    return d;
}

std::ostream& operator<<(std::ostream &os, const Date &d) {
    if (!d.valid) {
        os << "dayOfMonth: " << d.dayOfMonth << " month: " << d.month << " year: " << d.year << " cycle: " << d.precessionalCycle;
    } else {
        os << "invalid Date";
    }
    return os;
}

}
}
}
