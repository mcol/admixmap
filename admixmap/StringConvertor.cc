#include "StringConvertor.h"

using namespace std;

string StringConvertor::dequote(const string& str)
{
    if (str.size() >= 2 && str[0] == '"' && str[str.size() - 1] == '"') {
        return str.substr(1, str.size() - 2);
    }

    return str;
}

int StringConvertor::toInt(const string& str)
{
    return atoi(dequote(str).c_str());
}

double StringConvertor::toFloat(const string& str)
{
    return atof(dequote(str).c_str());
}

pair<int, int> StringConvertor::toIntPair(const string& str)
{
    const std::string tmp = dequote(str);

    // Try to find comma first
    string::const_iterator sci = find(tmp.begin(), tmp.end(), ',');
                    
    // ... the try find slash.
    if (sci == tmp.end()) sci = find(tmp.begin(), tmp.end(), '/');

    // If nothinf is found return default value.
    if (sci == tmp.end()) return make_pair(0, 0);
    
    const string allele0_str(tmp.begin(), sci);
    const string allele1_str(sci + 1, tmp.end());
    return make_pair(atoi(allele0_str.c_str()), atoi(allele1_str.c_str()));
}

bool StringConvertor::isMissingValue(const string& str)
{
    static string MISSING[] = {"#", ".", "NA"};
    return MISSING + 3 != find(MISSING, MISSING + 3, str);
}
