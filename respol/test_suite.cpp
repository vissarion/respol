#include <iostream>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

void test_res(std::string const& strin, std::string const& expected_result){
    std::string str = "./res_enum_d -v 0 < resultant_examples/";
    str = str + strin;
    auto str_result = exec(str.c_str());

    if (str_result == expected_result)
        std::cout << "SUCCESS" << std::endl;
    else
        std::cout << "FAIL" << std::endl;
}

void test_dis(std::string const& strin, std::string const& expected_result){
    std::string str = "./res_enum_d -v 0 -d < discriminant_examples/";
    str = str + strin;
    auto str_result = exec(str.c_str());

    if (str_result == expected_result)
        std::cout << "SUCCESS" << std::endl;
    else
        std::cout << "FAIL" << std::endl;
}

void test_sec(std::string const& strin, std::string const& expected_result){
    std::string str = "./res_enum_d -v 0 -s < secondary_examples/";
    str = str + strin;
    auto str_result = exec(str.c_str());

    if (str_result == expected_result)
        std::cout << "SUCCESS" << std::endl;
    else
        std::cout << "FAIL" << std::endl;
}

int main(){

    test_res("all_points_4.tmp", "20\n");
    test_res("bicubic_cayley.tmp", "6\n");
    test_res("cayley4_small.tmp", "19\n");
    test_res("cayley4_special.tmp", "23\n");
    test_res("cayley4.tmp", "119\n");
    test_res("CH_points_4.tmp", "20\n");
    test_res("example1mega.txt", "6\n");
    test_res("example1_2mega.txt", "4\n");

    test_dis("example2mega.txt", "8\n");
    test_dis("example3disc.txt", "7\n");

    test_sec("cube3.txt", "74\n");
    test_sec("simple.txt", "4\n");

    return 0;
}
