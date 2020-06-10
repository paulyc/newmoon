#ifndef NEWMOON_TEST_JD_CLOCK_HPP
#define NEWMOON_TEST_JD_CLOCK_HPP

#include <gtest/gtest.h>
#include "../src/jd_clock.hpp"
#include "../src/jd_clock.cpp"
#include "quadmath.h"

namespace {

TEST(JDClockTestSuite, TestTrue) {
    jd_clock::test_delta_t_lerp();
    auto jd = jd_clock::now();
    std::cerr << jd << std::endl;
    std::cerr << jd.time_since_epoch().count() << std::endl;
    const double jd_now = static_cast<double>(jd_clock::duration(jd.time_since_epoch()).count());
    std::cerr << jd_now << std::endl;
    EXPECT_EQ(true, true);
}

}

#endif
