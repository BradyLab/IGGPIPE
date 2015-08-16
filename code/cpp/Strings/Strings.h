/*
 *  Strings.h
 *
 *  Created by Ted Toal on 10/8/10.
 *
 *  String utility functions.
 */

#ifndef Strings_h
#define Strings_h

#ifndef stdhdr_h
#include "stdhdr.h"
#endif

#ifndef _STRING_H_
#include <string.h>
#endif

#ifndef _GLIBCXX_SSTREAM
#include <sstream>
#endif

#ifndef _GLIBCXX_STRING
#include <string>
#endif

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

using namespace std;

// A class to add enhancements to class string.
class S : public string
    {
    // Global empty string.
    public: static const string E;

    // Search this for searchStr and replace it with replaceStr, returning a
    // reference to the modified string.  If global is true, replace all occurrences.
    public: string& findReplace (const string& searchStr, const string& replaceStr,
        size_t* numFound = NULL, size_t startPos = 0, bool global = false);

    // Strip extension from filename and return it in 'extension' returning stripped file name.
    public: static string removeExtension(string filename, string& extension);

    // Constructors for each standard C++ type.
    S() : string() {}
    S(const string& str) : string(str) {}
    S(const char* str) : string(str) {}
    explicit S(const bool& val) : string() { ostringstream oss; oss << boolalpha << val; assign(oss.str()); }
    explicit S(const short& val) : string() { ostringstream oss; oss << val; assign(oss.str()); }
    explicit S(const unsigned short& val) : string() { ostringstream oss; oss << val; assign(oss.str()); }
    explicit S(const int& val) : string() { ostringstream oss; oss << val; assign(oss.str()); }
    explicit S(const unsigned int& val) : string() { ostringstream oss; oss << val; assign(oss.str()); }
    explicit S(const long& val) : string() { ostringstream oss; oss << val; assign(oss.str()); }
    explicit S(const unsigned long& val) : string() { ostringstream oss; oss << val; assign(oss.str()); }
    explicit S(const float& val) : string() { ostringstream oss; oss << val; assign(oss.str()); }
    explicit S(const double& val) : string() { ostringstream oss; oss << val; assign(oss.str()); }
    explicit S(const long double& val) : string() { ostringstream oss; oss << val; assign(oss.str()); }

    };

/*
    Some very useful standard exception class extensions.
*/
class invalidArgument : public invalid_argument
    {
    public: explicit invalidArgument(const S& arg, const S& val=S::E, const S& func=S::E, const S& clss=S::E)
        : invalid_argument(string(arg+": '"+val+"' in "+clss+"::"+func+"()")) {}
    };

class outOfRange : public out_of_range
    {
    public: explicit outOfRange(const S& arg, const S& val=S::E, const S& func=S::E, const S& clss=S::E)
        : out_of_range(string(arg+": '"+val+"' in "+clss+"::"+func+"()")) {}
    };

class runtimeError : public runtime_error
    {
    public: explicit runtimeError(const S& arg, const S& val=S::E, const S& func=S::E, const S& clss=S::E)
        : runtime_error(string(arg+": '"+val+"' in "+clss+"::"+func+"()")) {}
    };

#endif
