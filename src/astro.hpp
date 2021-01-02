/**
 * astro.hpp
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

#ifndef PAULYC_ASTRO_HPP
#define PAULYC_ASTRO_HPP

#include "lalgebra.hpp"
#include "calculus.hpp"
#include "jd_clock.hpp"
#include "ephemshelper.hpp"

#include <iostream>
#include <functional>

typedef std::function<long double(JPLEphems&, const jd_clock::time_point&)> f_type;

long double moonSunAngle(JPLEphems &ephems, const jd_clock::time_point &jd);
std::chrono::system_clock::time_point minFinder(JPLEphems &ephems, jd_clock::time_point &jd);

#endif /* PAULYC_ASTRO_HPP */
