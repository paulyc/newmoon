/**
 * newmoon.hpp - part of newmoon, moon phase calculator
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

#ifndef PAULYC_NEWMOON_HPP
#define PAULYC_NEWMOON_HPP

#include <memory>
#include <chrono>

#include "jd_clock.hpp"

class JPLEphems;

namespace paulyc {

class NewMoon
{
public:
    NewMoon() = default;
    ~NewMoon() = default;
    NewMoon(const NewMoon&) = delete;
    NewMoon& operator=(const NewMoon&) = delete;
    NewMoon(NewMoon&&);
    NewMoon& operator=(NewMoon&&);

    void init();
    std::chrono::system_clock::time_point next_moon();

private:
    std::unique_ptr<JPLEphems> _ephems;
    std::chrono::system_clock::time_point _t;
    jd_clock::time_point _jd;
};

}

#endif // PAULYC_NEWMOON_HPP
