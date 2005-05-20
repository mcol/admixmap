#include "StringConvertor.h"

using namespace std;

string StringConvertor::dequote(const string& str)
{
  if (str.size() >= 2 && (str[0] == '"' || str[0]== '\'') && (str[str.size() - 1] == '"' || (str[str.size() - 1] == '\''))) {
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

//pair<int, int> StringConvertor::toIntPair(const string& str)
// void StringConvertor::toIntPair(unsigned int *a,const string& str)
// {
//   cout<<"a= "<<a[0]<<" "<<a[1]<<endl;

//     const std::string tmp = dequote(str);
//     cout<<"tmp = "<<tmp<<endl;
//     a[0]=a[1]=0;
//     // Try to find comma first
//     string::const_iterator sci = find(tmp.begin(), tmp.end(), ',');
                    
//     // ... then try find slash.
//     if (sci == tmp.end()) sci = find(tmp.begin(), tmp.end(), '/');

//     // If nothing is found return default value.
//     if (sci == tmp.end()) {//return make_pair(0, 0);
//       a[0] = a[1] = 0;
//     }
//     else{
//       //cout<<tmp.begin()<<" "<<sci<<" "<<tmp.end()<<endl;
//       const string allele0_str(tmp.begin(), sci);
//       const string allele1_str(sci + 1, tmp.end());
//       //return make_pair(atoi(allele0_str.c_str()), atoi(allele1_str.c_str()));
 
//       a[0] = atoi(allele0_str.c_str());
//       a[1] = atoi(allele1_str.c_str());
//     }
//   cout<<"a= "<<a[0]<<" "<<a[1]<<endl;
//     //return a;
// }

void StringConvertor::toIntPair(unsigned int *a,const string& str)
{

    const std::string tmp = dequote(str);

    a[0]=a[1]=0;
  
    int i = tmp.find_first_of(",/"); 
    if( tmp.find_first_of(",/") < tmp.length()){
      a[0] = atoi(tmp.substr(0,i).c_str());
      a[1] = atoi(tmp.substr(i+1,tmp.length()-i).c_str());
    }
    else{
      a[0] = a[1] = 0;
    }
                    
}

void StringConvertor::toIntPair(unsigned short *a,const string& str)
{

    const std::string tmp = dequote(str);

    a[0]=a[1]=0;
  
    int i = tmp.find_first_of(",/"); 
    if( tmp.find_first_of(",/") < tmp.length()){
      a[0] = atoi(tmp.substr(0,i).c_str());
      a[1] = atoi(tmp.substr(i+1,tmp.length()-i).c_str());
    }
    else{
      a[0] = a[1] = 0;
    }
                    
}
bool StringConvertor::isMissingValue(const string& str)
{
    static string MISSING[] = {"#", ".", "NA"};
    return MISSING + 3 != find(MISSING, MISSING + 3, str);
}
