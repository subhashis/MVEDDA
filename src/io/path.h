#ifndef PATH_H_
#define PATH_H_

#include <string>
#include "edda_export.h"

namespace edda{
std::string EDDA_EXPORT getPath(const std::string &file);
std::string EDDA_EXPORT getFilename(const std::string &filepath);
std::string EDDA_EXPORT removeFileExtension(const std::string &filename);
bool EDDA_EXPORT isFilenameOnly(const std::string &filename);
std::string EDDA_EXPORT getFileExtension(const std::string &filename);
}
#endif


