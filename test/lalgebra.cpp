#ifndef NEWMOON_TEST_LALGEBRA_HPP
#define NEWMOON_TEST_LALGEBRA_HPP

#include <gtest/gtest.h>
#include "../src/lalgebra.hpp"

namespace {

static constexpr long double fuzzy_delta = 0.00001l;

template <typename T>
bool fuzzy_eq(T x, T y) {
    const long double delta = fuzzy_delta;
    return fabsl(x - y) < delta;
}

bool fuzzy_eq(const Nvec<3> &u, const Nvec<3> &v) {
    const long double delta = 0.00001l;
    for (std::size_t i = 0; i < 3; ++i) {
        if (fabsl(u.data[i] - v.data[i]) > delta) {
            return false;
        }
    }
    return true;
}

#define EXPECT_FUZZY_EQ(a,b) EXPECT_TRUE(fuzzy_eq((a), (b)))

TEST(LAlgebraTestSuite, TestTrue) {
    EXPECT_EQ(true, true);
}

class NvecTestFixture : public testing::Test
{
protected:
    static constexpr Nvec<3> u = {1.0l, 2.0l, 3.0l};
    static constexpr Nvec<3> v = {4.0l, 5.0l, 6.0l};
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

TEST(NxMmatrix, TestTrue) {
    NxMmatrix<3, 3> mat;
    Nvec<3> vec = {1.0l, 2.0l, 3.0l};
    Nvec<3> res = mat.mul(vec);
    EXPECT_TRUE(fuzzy_eq(vec, res));
}

}

#endif /* NEWMOON_TEST_LALGEBRA_HPP */
