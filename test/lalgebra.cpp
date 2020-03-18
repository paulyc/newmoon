#ifndef NEWMOON_TEST_LALGEBRA_HPP
#define NEWMOON_TEST_LALGEBRA_HPP

#include <gtest/gtest.h>
#include "../src/lalgebra.hpp"
#include <quadmath.h>

namespace {

static constexpr __float128 fuzzy_delta = 0.00001q;

template <typename T>
bool fuzzy_eq(T x, T y) {
    const __float128 delta = fuzzy_delta;
    if (fabsq(x - y) > delta) {
        std:: cerr << "fuzzy eq returning false on " << static_cast<long double>(x) <<
                      " diff " << static_cast<long double>(y) << " > " << static_cast<long double>(delta) << std::endl;
        return false;
    } else {
        return true;
    }
}

bool fuzzy_eq(const Nvec<3> &u, const Nvec<3> &v) {
    const __float128 delta = 0.00001q;
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
    EXPECT_TRUE(fuzzy_eq(MMM_2_PI, 2.0q*MMM_PI));
    EXPECT_TRUE(fuzzy_eq(MMM_PI_PI, MMM_PI * MMM_PI));
}

class NvecTestFixture : public testing::Test
{
protected:
    static constexpr Nvec<3> u = {1.0q, 2.0q, 3.0q};
    static constexpr Nvec<3> v = {4.0q, 5.0q, 6.0q};
    static constexpr Nvec<3> w = {MMM_PI, MMM_PI, MMM_PI};
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
    __float128 dotP = u.dotP(v);
    EXPECT_TRUE(fuzzy_eq(dotP, 4.0q + 10.0q + 18.0q));
}

TEST_F(NvecTestFixture, TestSelfDotP) {
    __float128 dotP = u.dotP(u);
    EXPECT_TRUE(fuzzy_eq(dotP, 1.0q + 4.0q + 9.0q));
}

TEST_F(NvecTestFixture, TestMag) {
    __float128 mag = w.mag();
    EXPECT_TRUE(fuzzy_eq(mag, sqrtq(MMM_PI*MMM_PI + MMM_PI*MMM_PI + MMM_PI*MMM_PI)));
}

TEST_F(NvecTestFixture, TestZeroVec) {
    EXPECT_TRUE(fuzzy_eq(TheZeroVector<__float128>.dotP(TheZeroVector<__float128>), 0.0q));
}

class NxMmatrixTestFixture : public testing::Test
{
protected:
    static constexpr Nvec<3> u = {1.0q, 2.0q, 3.0q};
    static constexpr NxMmatrix<3,3> m = {{
        {1.0q, 2.0q, 3.0q},
        {4.0q, 5.0q, 6.0q},
        {7.0q, 8.0q, 9.0q},
        }};
    static constexpr NxMmatrix<3,3> n = {{
        {10.0q, 20.0q, 30.0q},
        {40.0q, 50.0q, 60.0q},
        {70.0q, 80.0q, 90.0q},
        }};
};

TEST_F(NxMmatrixTestFixture, TestTrue) {
    EXPECT_TRUE(true);
}

TEST_F(NxMmatrixTestFixture, TestVecMul) {
    Nvec<3> v = m.mul(u);
    __float128 dotP0 = m.row(0).dotP(v);
    __float128 dotP1 = m.row(1).dotP(v);
    __float128 dotP2 = m.row(2).dotP(v);
    EXPECT_TRUE(fuzzy_eq(dotP0, 4.0q+9.0q));
    EXPECT_TRUE(fuzzy_eq(dotP1, 4.0q+10.0q+18.0q));
    EXPECT_TRUE(fuzzy_eq(dotP2, 7.0q+16.0q+27.0q));
    Nvec<3> u = {dotP0, dotP1, dotP2};
    EXPECT_TRUE(fuzzy_eq(v, u));
}

}

#endif /* NEWMOON_TEST_LALGEBRA_HPP */
