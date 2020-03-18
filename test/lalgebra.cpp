#ifndef NEWMOON_TEST_LALGEBRA_HPP
#define NEWMOON_TEST_LALGEBRA_HPP

#include <gtest/gtest.h>
#include "../src/lalgebra.hpp"

namespace {

static constexpr __float128 fuzzy_delta = 0.00001q;

template <typename T>
bool fuzzy_eq(T x, T y) {
    const long double delta = fuzzy_delta;
    if (fabsl(x - y) > delta) {
        std:: cerr << "fuzzy eq returning false on " << static_cast<long double>(x) << "diff" << static_cast<long double>(y) << " > " << delta << std::endl;
        return false;
    } else {
        return true;
    }
}

bool fuzzy_eq(const Nvec<3> &u, const Nvec<3> &v) {
    const long double delta = 0.00001l;
    for (std::size_t i = 0; i < 3; ++i) {
        if (fabsl(u.data[i] - v.data[i]) > delta) {
            std::cerr << "fuzzy eq returning false on " << static_cast<long double>(u.data[i]) <<
                         "diff" << static_cast<long double>(v.data[i]) << " > " << delta << std::endl;
            return false;
        }
    }
    return true;
}

#define EXPECT_FUZZY_EQ(a,b) EXPECT_TRUE(fuzzy_eq((a), (b)))

TEST(LAlgebraTestSuite, TestTrue) {
    EXPECT_EQ(true, true);
}

TEST(LAlgebraTestSuite, TestConstants) {
    EXPECT_EQ(MMM_2_PI, 2.0q*MMM_PI);
    EXPECT_EQ(MMM_PI_PI, MMM_PI * MMM_PI);
}

class NvecTestFixture : public testing::Test
{
protected:
    static constexpr Nvec<3> u = {1.0q, 2.0q, 3.0q};
    static constexpr Nvec<3> v = {4.0q, 5.0q, 6.0q};
    static constexpr Nvec<3> w = {MMM_PI, MMM_PI, MMM_PI};
};

TEST_F(NvecTestFixture, TestAdd) {
    EXPECT_FUZZY_EQ(u.add(v), Nvec<3>({5.0q, 7.0q, 9.0q}));
}

TEST_F(NvecTestFixture, TestAssign) {
    Nvec w = u;
    EXPECT_FUZZY_EQ(w, u);
}

TEST_F(NvecTestFixture, TestAddAcc) {
    Nvec w = u;
    w.addacc(v);
    EXPECT_FUZZY_EQ(w, u.add(v));
}

TEST_F(NvecTestFixture, TestDotP) {
    __float128 dotP = u.dotP(v);
    EXPECT_FUZZY_EQ(dotP, 4.0q + 10.0q + 18.0q);
}

TEST_F(NvecTestFixture, TestSelfDotP) {
    __float128 dotP = u.dotP(u);
    EXPECT_FUZZY_EQ(dotP, 1.0q + 4.0q + 9.0q);
}

TEST_F(NvecTestFixture, TestMag) {
    __float128 mag = w.mag();
    EXPECT_FUZZY_EQ(mag, sqrtl(MMM_PI*MMM_PI + MMM_PI*MMM_PI + MMM_PI*MMM_PI));
}

TEST_F(NvecTestFixture, TestMag_PI_PI) {
    __float128 mag = w.mag();
    EXPECT_FUZZY_EQ(mag, sqrtl(3.0q * MMM_PI_PI));
}

TEST_F(NvecTestFixture, TestZeroVec) {
    EXPECT_FUZZY_EQ(TheZeroVector<__float128>.dotP(TheZeroVector<__float128>), 0.0q);
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
    EXPECT_FUZZY_EQ(dotP0, 4.0q+9.0q);
    EXPECT_FUZZY_EQ(dotP1, 4.0q+10.0q+18.0q);
    EXPECT_FUZZY_EQ(dotP2, 7.0q+16.0q+27.0q);
    Nvec<3> u = {dotP0, dotP1, dotP2};
    EXPECT_FUZZY_EQ(v, u);
}

}

#endif /* NEWMOON_TEST_LALGEBRA_HPP */
