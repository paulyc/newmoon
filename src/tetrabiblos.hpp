/**
 * tetrabiblos.hpp
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

#ifndef PAULYC_TETRABIBLOS_HPP
#define PAULYC_TETRABIBLOS_HPP

#include "astro.hpp"

namespace github {
namespace paulyc {
namespace tetrabiblos {

enum Month {
    NONE,
    CAPRICORNUS = 1,
    AQUARIUS,
    PISCES,
    ARIES,
    TAURUS,
    CANCER,
    GEMINI,
    VIRGO,
    LEO,
    LIBRA,
    SCORPIO,
    SAGITTARIUS,
};

struct Date {
    bool valid;
    int precessionalCycle;
    int year;
    Month month;
    int dayOfMonth;
};

Date getDate(JPLEphems &ephems, const std::chrono::system_clock::time_point &tp);
std::ostream& operator<<(std::ostream &os, const Date &d);

}
}
}

#endif /* PAULYC_TETRABIBLOS_HPP */
