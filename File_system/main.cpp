#include <iostream>
#include <fstream>
#include <cstdlib>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

int main()
{
    fs::create_directories("sandbox/1/2/a");
    fs::create_directory("sandbox/1/2/b");
    fs::permissions("sandbox/1/2/b", fs::perms::remove_perms | fs::perms::others_all);
    fs::create_directory("sandbox/1/2/c", "sandbox/1/2/b");
    std::system("ls -l sandbox/1/2");
    fs::remove_all("sandbox");
}
