/**
 * lalgebra.cpp tests
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

#include <gtest/gtest.h>
#include "../src/jd_clock.hpp"

namespace {

TEST(jd_clock_test_suite, test_delta_t_lerp) {
    EXPECT_TRUE(jd_clock::delta_t_lerp(-700) > 17190);
    EXPECT_TRUE(jd_clock::delta_t_lerp(1250) > 200 && jd_clock::delta_t_lerp(1250) < 1570);
    EXPECT_TRUE(jd_clock::delta_t_lerp(2030) > 68.97);
}

}
