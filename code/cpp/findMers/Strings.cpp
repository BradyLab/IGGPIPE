/*
 *  Strings.cpp
 *
 *  Created by Ted Toal on 10/8/10.
 */

#include "Strings.h"

const string S::E;

// Search and replace.
string& S::findReplace (const string& searchStr, const string& replaceStr,
    size_t* numFound, size_t startPos, bool global)
    {
    size_t num = 0;
    size_t len = searchStr.length();
    while (true)
        {
        size_t pos = find(searchStr, startPos);
        if (pos == string::npos)
            break;
        num++;
        replace(pos, len, replaceStr);
        startPos = pos+len;
        if (!global)
            break;
        }
    if (numFound != NULL)
        *numFound = num;
    return(*this);
    }

string S::removeExtension(string filename, string& extension)
    {
	extension = "";
	size_t i = filename.rfind('.');
	if (i == string::npos)
		return(filename);
	extension = filename.substr(i);
	if (i == 0)
		return("");
	return(filename.substr(0, i));
    }

