#include <gtest/gtest.h>

#define private public
#define protected public


#include <chrono>
#include <iostream>

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}