#include <iostream>
#include "spdlog/spdlog.h"

int main(int argc, const char **argv)
{
    std::cout << "Helo world!" << std::endl;
    std::cout << "ABCD!" << std::endl;
    spdlog::info("Hudhhsud");
    spdlog::warn("Hudhhsud");
    spdlog::debug("Hudhhsud");
    spdlog::error("Hudhhsud");
    spdlog::critical("Hudhhsud");
    std::cout << "Second!" << std::endl;

    return 0;
}