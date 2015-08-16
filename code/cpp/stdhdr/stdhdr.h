/*
 *  stdhdr.h
 *
 *  Created by Ted Toal on 10/8/10.
 *
 *  Standard header file for Mac OSX/XCode environment.
 */

#ifndef stdhdr_h
#define stdhdr_h

#ifndef _GLIBCXX_STDEXCEPT
#include <stdexcept>
#endif

/*
    Assume std namespace.
*/

using namespace std;

/*
    Define some common mathematical constant values and macros.
*/

#ifndef PI
#define PI ((double) 3.141592654)
#endif
#ifndef DEG_PER_RAD
#define DEG_PER_RAD 57.29577951
#endif
#ifndef ROUND
#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#endif
#ifndef MAX
#define MAX(x,y) ((x)>=(y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x)<=(y) ? (x) : (y))
#endif

/*
    Define a single type name for a long double.
*/

typedef long double ldouble;

/*
    Define pseudo "sizeof" macros applying to array types.
*/

#define elementsize(array)  (sizeof(*array))
#define arraysize(array)    (sizeof(array)/sizeof(*array))
#define arrayend(array)     (array+sizeof(array)/sizeof(*array))

/*
    Define pseudo program control keyword.
*/

#define loop for(;;)

/*
    Define some standard "byte" types.
*/

typedef unsigned char ubyte;
typedef unsigned short ubyte2;
typedef unsigned long ubyte4;
typedef unsigned long long ubyte8;

/*
    Define common logical constant values.
*/

#ifndef YES
#define YES     1
#endif

#ifndef ON
#define ON      1
#endif

#ifndef SET
#define SET     1
#endif

#ifndef HIGH
#define HIGH    1
#endif

#ifndef ENABLE
#define ENABLE  1
#endif

#ifndef NO
#define NO      0
#endif

#ifndef OFF
#define OFF     0
#endif

#ifndef CLEAR
#define CLEAR   0
#endif

#ifndef LOW
#define LOW     0
#endif

#ifndef DISABLE
#define DISABLE 0
#endif

/*
    Define a macro to return the least significant byte of its argument.
*/

#define LSB(x) ((ubyte)((x)&0xFF))

/*
    Define macros for combining and splitting bytes and words.
*/

#ifndef LO_BYTE
#define LO_BYTE(w)          ((ubyte)(w))
#endif
#ifndef HI_BYTE
#define HI_BYTE(w)          ((ubyte)((uint)(w) >> 8))
#endif
#ifndef LO_WORD
#define LO_WORD(l)          ((ubyte2)(l))
#endif
#ifndef HI_WORD
#define HI_WORD(l)          ((ubyte2)((ubyte4)(l) >> 16))
#endif
#ifndef BYTE1_OF_LONG
#define BYTE1_OF_LONG(l)    ((ubyte)(l))
#endif
#ifndef BYTE2_OF_LONG
#define BYTE2_OF_LONG(l)    ((ubyte)((ubyte4)(l) >> 8))
#endif
#ifndef BYTE3_OF_LONG
#define BYTE3_OF_LONG(l)    ((ubyte)((ubyte4)(l) >> 16))
#endif
#ifndef BYTE4_OF_LONG
#define BYTE4_OF_LONG(l)    ((ubyte)((ubyte4)(l) >> 24))
#endif
#ifndef MAKE_LONG
#define MAKE_LONG(low, high)    ((long)(((ubyte2)(low)) | (((ubyte4)((ubyte2)(high))) << 16)))
#endif
#ifndef MAKE_UBYTE2
#define MAKE_UBYTE2(low, high)  ((ubyte2)(((ubyte2)((ubyte)(low))) | (((ubyte2)((ubyte)(high))) << 8)))
#endif
#ifndef MAKE_UBYTE4
#define MAKE_UBYTE4(low, high)  ((ubyte4)MAKE_LONG(low, high))
#endif

/*
    Enumeration that is useful for functions that must be able to access the
    first, next, previous, or last element in some last.
*/
enum eWhichInSeries
    {
    OBJECT_FIRST,
    OBJECT_NEXT,
    OBJECT_PREVIOUS,
    OBJECT_LAST
    };

/*
    Macros to increment/decrement an enumeration type.
*/
#define eINCR(e, eType) ((e) = (eType) (((unsigned) (e)) + 1))
#define eDECR(e, eType) ((e) = (eType) (((unsigned) (e)) - 1))

/*
    Macros to add an integer to/subtract one from an enumeration type and
    return the result as an enumeration type.
*/
#define ePLUS(e, i, eType) ((eType) (((unsigned) (e)) + (i)))
#define eMINUS(e, i, eType) ((eType) (((unsigned) (e)) - (i)))

/*
    Macro to return the difference between two enumeration types as an integer.
*/
#define eDIFF(e1, e2) ((unsigned) (e1) - (unsigned) (e2))

#endif
