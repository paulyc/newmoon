#ifndef NEWMOON_TEST_LALGEBRA_HPP
#define NEWMOON_TEST_LALGEBRA_HPP

#include <gtest/gtest.h>
#include "../src/lalgebra.hpp"

namespace {

TEST(LAlgebraTestSuite, TestTrue) {
    EXPECT_EQ(true, true);
}

bool fuzzy_eq(Nvec<3> &u, Nvec<3> &v) {
    const long double delta = 0.00001l;
    for (std::size_t i = 0; i < 3; ++i) {
        if (fabsl(u.data[i] - v.data[i]) > delta) {
            return false;
        }
    }
    return true;
}

TEST(NxMmatrix, TestTrue) {
    NxMmatrix<3, 3> mat;
    Nvec<3> vec = {1.0l, 2.0l, 3.0l};
    Nvec<3> res = mat.mul(vec);
    EXPECT_TRUE(fuzzy_eq(vec, res));
}

}

#endif /* NEWMOON_TEST_LALGEBRA_HPP */
