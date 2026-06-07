/** \file
 *
 *
 * \brief Typedefs for data types so that sizes are known.
 *
 *
 */

#include<cstdint>
#include<stack>

#ifndef Windows_DataType
#define Windows_DataType

/// 64-bit data type
typedef uint64_t    QWORD;

/// 32-bit data type
typedef uint32_t    DWORD;

/// Boolean
typedef int BOOL;

/// 8-bit data type
typedef unsigned char BYTE;

/// 16-bit data type
typedef uint16_t    WORD;

/// 32-bit float
typedef float FLOAT;

/// long, presumably meaning 64-bit integer.
typedef long LONG;

/// 64-bit double
typedef double DOUBLE;


// The following are needed for Radarsat-2 because the SLC files
// store their values as 2s complement signed data.

/// signed 8-bit for reading 2s complement signed data.
typedef char       SIGNEDBYTE;

/// signed 16-bit for reading 2s complement signed data.
typedef short      SIGNEDWORD;

/// 2D point.
typedef struct
{
	int Row;
	int Col;
} POINTex;


/// 2D point, double.
typedef struct
{
	double Row;
	double Col;
} POINTDex;

/// Rectangle, int, (x1, y1) = top left, (x2, y2) = bottom right
typedef struct
{
    int x1;
    int y1;
    int x2;
    int y2;
} Rect_i;

/// Rectangle, double, (x1, y1) = top left, (x2, y2) = bottom right
typedef struct
{
    double x1;
    double y1;
    double x2;
    double y2;
} Rect_d;

typedef Rect_i Rect;

/// Defines offsets from center pixel to get to each of 8 pixels in neighbourhood.
/**
 * The order that the neighbours appear in this array in is specified below:
 * -------
 * |1|2|3|
 * -------
 * |0| |4|
 * -------
 * |7|6|5|
 * -------
 */
const static POINTex neighbor8[8]=
{
	{ 0,-1},
	{-1,-1},
	{-1, 0},
	{-1, 1},
	{ 0, 1},
	{ 1, 1},
	{ 1, 0},
	{ 1,-1},
};


/// Defines offsets from center pixel to get to each of 4 pixels in neighbourhood.
/**
 * The order that the neighbours appear in this array in is specified below:
 * -------
 * | |1| |
 * -------
 * |0| |2|
 * -------
 * | |3| |
 * -------
 */
const static POINTex neighbor4[4]=
{
	{ 0,-1},
	{-1, 0},
	{ 0, 1},
	{ 1, 0},
};


/// Defines offsets from center pixel to get to each of 5 pixels in neighbourhood.
/**
 * The order that the neighbours appear in this array in is specified below:
 * -------
 * | |1| |
 * -------
 * |0|4|2|
 * -------
 * | |3| |
 * -------
 * Note that this includes the center pixel itself.
 */
const static POINTex neighbor5[5]=
{
	{ 0,-1},
	{-1, 0},
	{ 0, 1},
	{ 1, 0},
	{ 0, 0},
};


/// Defines offsets from center pixel to get to each of 16 pixels in neighbourhood.
/**
 * The order that the neighbours appear in this array in is specified below:
 * ----------------
 * | 2| 3| 4| 5| 6|
 * ----------------
 * | 1|  |  |  | 7|
 * ----------------
 * | 0|  |  |  | 8|
 * ----------------
 * |15|  |  |  | 9|
 * ----------------
 * |14|13|12|11|10|
 * ----------------
 */
const static POINTex neighbor16[16]=
{
	{ 0,-2},
	{-1,-2},
	{-2,-2},
	{-2,-1},
	{-2, 0},
	{-2, 1},
	{-2, 2},
	{-1, 2},
	{ 0, 2},
	{ 1, 2},
	{ 2, 2},
	{ 2, 1},
	{ 2, 0},
	{ 2,-1},
	{ 2,-2},
	{ 1,-2},
};


#endif
