
#include <io/path.h>

using namespace std;

namespace edda{

std::string getPath(const std::string &filepath)
{
        size_t i = filepath.find_last_of("/\\");
        if (i==string::npos)
                return "";
        else
                return filepath.substr(0, i+1);
}
std::string getFilename(const std::string &filepath)
{
        size_t i = filepath.find_last_of("/\\");
        if (i==string::npos)
                return filepath;
        else
                return filepath.substr(i+1);

}
string removeFileExtension(const string &filename)
{
        size_t i = filename.find_first_of('.');
        return filename.substr(0, i);
}

bool isFilenameOnly(const string &filename)
{
        return filename.find_first_of("/\\")==string::npos;
}

string getFileExtension(const string &filename)
{
        size_t i = filename.find_last_of('.');
        return filename.substr(i+1);
}

}
