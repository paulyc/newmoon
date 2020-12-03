#include <gtest/gtest.h>

//static_assert(sizeof(long double) == sizeof(__float128));

TEST(MainTestSuite, LongDoubleIsFloat128) {
    EXPECT_EQ(sizeof(long double), sizeof(__float128));
}

int main(int argc, char *argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
