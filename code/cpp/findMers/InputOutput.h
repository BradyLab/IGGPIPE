/*
 *  InputOutput.h
 *
 *  Created by Ted Toal on 5/18/11.
 *
 *  Input/output utility functions.
 */

#ifndef InputOutput_h
#define InputOutput_h

#ifndef stdhdr_h
#include "stdhdr.h"
#endif

#ifndef Strings_h
#include "Strings.h"
#endif

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif

#ifndef _GLIBCXX_IOMANIP
#include <iomanip>
#endif

using namespace std;

// Text file output.
//
// ofstream is a public base class, so all ofstream functions are available.
//
// Example of use: OutputFile o("x.txt"); o << "Hello world." << endl; o.Close();
class OutputFile : public ofstream
    {
    public: string fileName;
    public: OutputFile() : ofstream() {}
    public: OutputFile(const char* filename) : ofstream(), fileName(filename)
        { open(filename, ios::trunc); }
    public: OutputFile(string filename) : ofstream(), fileName(filename)
        { open(filename.c_str(), ios::trunc); }
    public: bool Create(const char* filename, bool text = true)
        { fileName = filename; open(filename, (text ? ios::trunc : ios::trunc|ios::binary)); return(is_open()); }
    public: bool Create(string filename, bool text = true)
        { fileName = filename; open(filename.c_str(), (text ? ios::trunc : ios::trunc|ios::binary)); return(is_open()); }
    public: bool OpenAppend(const char* filename, bool text = true)
        { fileName = filename; open(filename, (text ? ios::app : ios::app|ios::binary)); return(is_open()); }
    public: bool OpenAppend(string filename, bool text = true)
        { fileName = filename; open(filename.c_str(), (text ? ios::app : ios::app|ios::binary)); return(is_open()); }
    public: void Close(void) { close(); }
    public: ~OutputFile() { close(); }
    
    // Functions for outputting name=value in a form consistent with the R language source() function.
    // Applying source() to the file created this way defines variables in R named "name" with value "value",
    // providing an easy way to transfer specific variable values from C++ to R.
    public: void Rsource_int(const char* name, int value)
        { *this << name << " = " << value << endl; }
    public: void Rsource_unsigned(const char* name, unsigned value)
        { *this << name << " = " << value << endl; }
    public: void Rsource_double(const char* name, double value)
        { *this << name << " = " << value << endl; }
    public: void Rsource_string(const char* name, const char* value)
        { *this << name << " = "; putString(value); *this << endl; }
    // For vectors, if elementNames is NULL, a simple vector is output, else a structure is output.
    public: void Rsource_int_vector(const char* name, size_t size, const char** elementNames, int* values);
    public: void Rsource_unsigned_vector(const char* name, size_t size, const char** elementNames, unsigned* values);
    public: void Rsource_double_vector(const char* name, size_t size, const char** elementNames, double* values);
    public: void Rsource_string_vector(const char* name, size_t size, const char** elementNames, const char** values);
    // A different way to output vectors, by specifying each element's name, type, and value individually after
    // the "name" argument, i.e. name, name1, type1, value1, name2, type2, value2, ..., NULL.  Last element name must
    // be NULL to terminate the list.  The type arguments must be one of "int", "unsigned", "double", or "string".
    // This function always outputs a structure.
    public: void Rsource_vector(const char* name, ...);
    // Helper functions used by above.
    protected: void putPreVector(const char* name, const char** elementNames);
    protected: void putPostVector(size_t size, const char** elementNames);
    protected: void putString(const char* value); // Handles strings with ", \, and other special chars in them.
    };

// Text file input.  Similar to above class for file output.
//
// Example of use: line s; InputFile i("x.txt"); getline(i, s); i.Close();
class InputFile : public ifstream
    {
    public: string fileName;
    public: InputFile() : ifstream() {}
    public: InputFile(const char* filename) : ifstream(), fileName(filename)
        { open(filename); }
    public: InputFile(string filename) : ifstream(), fileName(filename)
        { open(filename.c_str()); }
    public: bool Open(const char* filename, bool text = true)
        { fileName = filename; open(filename, (text ? ios::in : ios::in|ios::binary)); return(is_open()); }
    public: bool Open(string filename, bool text = true)
        { fileName = filename; open(filename.c_str(), (text ? ios::in : ios::in|ios::binary)); return(is_open()); }
    public: void Close(void) { close(); }
    public: ~InputFile() { Close(); }
    };
    
#endif
