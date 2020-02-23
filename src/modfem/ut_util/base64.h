#include <string>

std::string& base64_encode(unsigned char const* , unsigned int len, std::string &ret);
std::string& base64_decode(std::string const& s, std::string& ret);
