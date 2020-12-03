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
#include "../src/lalgebra.hpp"

namespace {

static constexpr long double fuzzy_delta = 0.00001q;

template <typename T>
bool fuzzy_eq(T x, T y) {
    const long double delta = fuzzy_delta;
    if (fabsq(x - y) > delta) {
        std:: cerr << "fuzzy eq returning false on " << static_cast<long double>(x) <<
                      " diff " << static_cast<long double>(y) << " > " << static_cast<long double>(delta) << std::endl;
        return false;
    } else {
        return true;
    }
}

bool fuzzy_eq(const Nvec<3> &u, const Nvec<3> &v) {
    const long double delta = 0.00001q;
    for (std::size_t i = 0; i < 3; ++i) {
        if (fabsq(u.data[i] - v.data[i]) > delta) {
            std::cerr << "fuzzy eq returning false on " << static_cast<long double>(u.data[i]) <<
                         " diff " << static_cast<long double>(v.data[i]) << " > " << static_cast<long double>(delta) << std::endl;
            return false;
        }
    }
    return true;
}

TEST(LAlgebraTestSuite, TestTrue) {
    EXPECT_EQ(true, true);
}

TEST(LAlgebraTestSuite, TestConstants) {
    EXPECT_TRUE(fuzzy_eq(MMM_2_PI, 2.0q * MMM_PI));
    EXPECT_TRUE(fuzzy_eq(MMM_4_PI, 4.0q * MMM_PI));
}

class NvecTestFixture : public testing::Test
{
public:
    NvecTestFixture() : testing::Test(), u({{1.0q, 2.0q, 3.0q}}), v({{4.0q, 5.0q, 6.0q}}), w({{MMM_PI, MMM_PI, MMM_PI}}) {}
protected:
    Nvec<3> u,v,w;
};

TEST_F(NvecTestFixture, TestAdd) {
    EXPECT_TRUE(fuzzy_eq(u.add(v), Nvec<3>({5.0q, 7.0q, 9.0q})));
}

TEST_F(NvecTestFixture, TestAssign) {
    Nvec<3> w = u;
    EXPECT_TRUE(fuzzy_eq(w, u));
}

TEST_F(NvecTestFixture, TestAddAcc) {
    Nvec<3> w = u;
    w.addacc(v);
    EXPECT_TRUE(fuzzy_eq(w, u.add(v)));
}

TEST_F(NvecTestFixture, TestDotP) {
    long double dotP = u.dotP(v);
    EXPECT_TRUE(fuzzy_eq(dotP, 4.0q + 10.0q + 18.0q));
}

TEST_F(NvecTestFixture, TestSelfDotP) {
    long double dotP = u.dotP(u);
    EXPECT_TRUE(fuzzy_eq(dotP, 1.0q + 4.0q + 9.0q));
}

TEST_F(NvecTestFixture, TestMag) {
    long double mag = w.mag();
    EXPECT_TRUE(fuzzy_eq(mag, sqrtq(MMM_PI*MMM_PI + MMM_PI*MMM_PI + MMM_PI*MMM_PI)));
}

TEST_F(NvecTestFixture, TestZeroVec) {
    EXPECT_TRUE(fuzzy_eq(TheZeroVector<long double>.dotP(TheZeroVector<long double>), 0.0q));
}

class NxMmatrixTestFixture : public testing::Test
{
};

TEST_F(NxMmatrixTestFixture, TestTrue) {
    EXPECT_TRUE(true);
}

#if 0
TEST_F(NxMmatrixTestFixture, TestVecMul) {
    NxMmatrix<3,3> m({
        {1.0q, 2.0q, 3.0q},
        {4.0q, 5.0q, 6.0q},
        {7.0q, 8.0q, 9.0q},
    });
    NxMmatrix<3,3> n({
        {10.0q, 20.0q, 30.0q},
        {40.0q, 50.0q, 60.0q},
        {70.0q, 80.0q, 90.0q},
    });
    Nvec<3> u({1.0q, 2.0q, 3.0q});
    Nvec<3> v = m.mul(u);
    long double dotP0 = m.row(0).dotP(u);
    long double dotP1 = m.row(1).dotP(u);
    long double dotP2 = m.row(2).dotP(u);
    EXPECT_TRUE(fuzzy_eq(dotP0, 1.0q+4.0q+9.0q));
    EXPECT_TRUE(fuzzy_eq(dotP1, 4.0q+10.0q+18.0q));
    EXPECT_TRUE(fuzzy_eq(dotP2, 7.0q+16.0q+27.0q));
    Nvec<3> z = {{dotP0, dotP1, dotP2}};
    EXPECT_TRUE(fuzzy_eq(v, z));
}
#endif

}
