#ifndef NEWMOON_TEST_LALGEBRA_HPP
#define NEWMOON_TEST_LALGEBRA_HPP

#include <gtest/gtest.h>
#include "../src/lalgebra.hpp"

namespace {

static constexpr long double fuzzy_delta = 0.00001l;

template <typename T>
bool fuzzy_eq(T x, T y) {
    const long double delta = fuzzy_delta;
    if (fabsl(x - y) > delta) {
        std:: cerr << "fuzzy eq returning false on " << x << "diff" << y << " > " << delta << std::endl;
        return false;
    } else {
        return true;
    }
}

bool fuzzy_eq(const Nvec<3> &u, const Nvec<3> &v) {
    const long double delta = 0.00001l;
    for (std::size_t i = 0; i < 3; ++i) {
        if (fabsl(u.data[i] - v.data[i]) > delta) {
            std::cerr << "fuzzy eq returning false on " << u.data[i] << 'diff' << v.data[i] << " > " << delta << std::endl;
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
    EXPECT_EQ(MMM_2_PI, 2.0l*MMM_PI);
    EXPECT_EQ(MMM_4_PI, 4.0l*MMM_PI);
    EXPECT_EQ(MMM_HALF_PI, 0.5l*MMM_PI);
    EXPECT_EQ(MMM_QUARTER_PI, 0.25l*MMM_PI);
    EXPECT_EQ(MMM_PI_PI, MMM_4_PI * MMM_QUARTER_PI);
}

class NvecTestFixture : public testing::Test
{
protected:
    static constexpr Nvec<3> u = {1.0l, 2.0l, 3.0l};
    static constexpr Nvec<3> v = {4.0l, 5.0l, 6.0l};
    static constexpr Nvec<3> w = {MMM_PI, MMM_PI, MMM_PI};
};

TEST_F(NvecTestFixture, TestAdd) {
    EXPECT_FUZZY_EQ(u.add(v), Nvec<3>({5.0l, 7.0l, 9.0l}));
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
    long double dotP = u.dotP(v);
    EXPECT_FUZZY_EQ(dotP, 4.0l + 10.0l + 18.0l);
}

TEST_F(NvecTestFixture, TestSelfDotP) {
    long double dotP = u.dotP(u);
    EXPECT_FUZZY_EQ(dotP, 1.0l + 4.0l + 9.0l);
}

TEST_F(NvecTestFixture, TestMag) {
    long double mag = w.mag();
    EXPECT_FUZZY_EQ(mag, sqrtl(MMM_PI*MMM_PI + MMM_PI*MMM_PI + MMM_PI*MMM_PI));
}

TEST_F(NvecTestFixture, TestMag_PI_PI) {
    long double mag = w.mag();
    EXPECT_FUZZY_EQ(mag, sqrtl(3.0l * MMM_PI_PI));
}

TEST_F(NvecTestFixture, TestZeroVec) {
    EXPECT_FUZZY_EQ(TheZeroVector<long double>.dotP(TheZeroVector<long double>), 0.0l);
}

class NxMmatrixTestFixture : public testing::Test
{
protected:
    static constexpr Nvec<3> u = {1.0l, 2.0l, 3.0l};
    static constexpr NxMmatrix<3,3> m = {{
        {1.0l, 2.0l, 3.0l},
        {4.0l, 5.0l, 6.0l},
        {7.0l, 8.0l, 9.0l},
        }};
    static constexpr NxMmatrix<3,3> n = {{
        {10.0l, 20.0l, 30.0l},
        {40.0l, 50.0l, 60.0l},
        {70.0l, 80.0l, 90.0l},
        }};
};

TEST_F(NxMmatrixTestFixture, TestTrue) {
    EXPECT_TRUE(true);
}

TEST_F(NxMmatrixTestFixture, TestVecMul) {
    Nvec<3> v = m.mul(u);
    long double dotP0 = m.row(0).dotP(v);
    long double dotP1 = m.row(1).dotP(v);
    long double dotP2 = m.row(2).dotP(v);
    EXPECT_FUZZY_EQ(dotP0, 4.0l+9.0l);
    EXPECT_FUZZY_EQ(dotP1, 4.0l+10.0l+18.0l);
    EXPECT_FUZZY_EQ(dotP2, 7.0l+16.0l+27.0l);
    Nvec<3> u = {dotP0, dotP1, dotP2};
    EXPECT_FUZZY_EQ(v, u);
}

}

#endif /* NEWMOON_TEST_LALGEBRA_HPP */
