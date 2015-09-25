/*
 *  InputOutput.cpp
 *
 *  Created by Ted Toal on 5/18/11.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "InputOutput.h"

void OutputFile::Rsource_int_vector(const char* name, size_t size, const char** elementNames, int* values)
    {
    putPreVector(name, elementNames);
    for (size_t i = 0; i < size; i++)
        {
        if (i > 0)
            *this << ",";
        *this << values[i];
        }
    putPostVector(size, elementNames);
    }

void OutputFile::Rsource_unsigned_vector(const char* name, size_t size, const char** elementNames, unsigned* values)
    {
    putPreVector(name, elementNames);
    for (size_t i = 0; i < size; i++)
        {
        if (i > 0)
            *this << ",";
        *this << values[i];
        }
    putPostVector(size, elementNames);
    }

void OutputFile::Rsource_double_vector(const char* name, size_t size, const char** elementNames, double* values)
    {
    putPreVector(name, elementNames);
    for (size_t i = 0; i < size; i++)
        {
        if (i > 0)
            *this << ",";
        *this << values[i];
        }
    putPostVector(size, elementNames);
    }

void OutputFile::Rsource_string_vector(const char* name, size_t size, const char** elementNames, const char** values)
    {
    putPreVector(name, elementNames);
    for (size_t i = 0; i < size; i++)
        {
        if (i > 0)
            *this << ",";
        putString(values[i]);
        }
    putPostVector(size, elementNames);
    }

void OutputFile::Rsource_vector(const char* name, ...)
    {
    *this << name << " = structure(c(";
    // Loop twice through the argument list, first time outputting element values, second time element names.
    for (int i = 0; i < 2; i++)
        {
        va_list args;
        va_start(args, name);
        const char* elementName = va_arg(args, const char*);
        // Loop through argument list.  "first" lets us know when to put commas after a previous value.
        for (bool first = true; elementName != NULL; first = false)
            {
            if (!first)
                *this << ",";
            // Second pass outputs element names.
            if (i == 1) // second pass
                putString(elementName);
            // Both passes must get value arguments, but only second pass outputs them.
            const char* type = va_arg(args, const char*);
            if (strcmp(type, "int") == 0)
                {
                int v = va_arg(args, int);
                if (i == 0) *this << v;
                }
            else if (strcmp(type, "unsigned") == 0)
                {
                unsigned v = va_arg(args, unsigned);
                if (i == 0) *this << v;
                }
            else if (strcmp(type, "double") == 0)
                {
                double v = va_arg(args, double);
                if (isnan(v))
                    v = 0; // Underflow causes this.  What is the RIGHT way to handle this???
                if (i == 0) *this << v;
                }
            else if (strcmp(type, "string") == 0)
                {
                const char* v = va_arg(args, const char*);
                if (i == 0) putString(v);
                }
            elementName = va_arg(args, const char*); // Next element name.
            }
        va_end(args);
        if (i == 0) *this << "),.Names=c(";
        }
    *this << "))" << endl;
    }
    
void OutputFile::putPreVector(const char* name, const char** elementNames)
    {
    *this << name << " = ";
    if (elementNames != NULL)
        *this << "structure(";
    *this << "c(";
    }

void OutputFile::putPostVector(size_t size, const char** elementNames)
    {
    if (elementNames != NULL)
        {
        *this << "),.Names=c(";
        for (size_t i = 0; i < size; i++)
            {
            if (i > 0)
                *this << ",";
            putString(elementNames[i]);
            }
        *this << ")";
        }
    *this << ")" << endl;
    }

void OutputFile::putString(const char* value)
    {
    *this << "\"";
    while (*value != 0)
        {
        switch (*value)
            {
            case '\\': *this << "\\\\"; break; // this outputs \\
            case '\"': *this << "\\\""; break; // this outputs \"
            case '\n': *this << "\\n"; break;
            case '\r': *this << "\\r"; break;
            case '\t': *this << "\\t"; break;
            case '\b': *this << "\\b"; break;
            case '\a': *this << "\\a"; break;
            case '\f': *this << "\\f"; break;
            case '\v': *this << "\\v"; break;
            // We won't deal with other special characters.
            default: *this << *value; break;
            }
        value++;
        }
    *this << "\"";
    }

