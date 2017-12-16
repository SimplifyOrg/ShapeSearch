/* Header file stlparser.h  */
#ifndef STLPARSER_H
#define STLPARSER_H

# include <cstdlib>
# include <iostream>
#include <string>
#include <cstring>
#include <vector>

using namespace std;

//export functions
#if defined(_MSC_VER)
    //  Microsoft
    #define EXPORT __declspec(dllexport)
    #define IMPORT __declspec(dllimport)
#elif defined(__GNUC__)
    //  GCC
    #define EXPORT __attribute__((visibility("default")))
    #define IMPORT
#else
    //  do nothing and hope for the best?
    #define EXPORT
    #define IMPORT
    #pragma warning Unknown dynamic link import/export semantics.
#endif
//export function complete

class v3
{
public:

	v3();
	v3(char* bin);
	v3(double x, double y, double z)
	{
		m_x = x; m_y = y; m_z = z;
	}

	~v3();
private:
	double m_x, m_y, m_z;
};

class triangle
{
public:

	triangle(v3 p1, v3 p2, v3 p3, v3 n)
	{
		m_p1 = p1; m_p2 = p2; m_p3 = p3; m_n = n;
	}
	~triangle();
private:
	v3 m_p1, m_p2, m_p3, m_n;
};

EXPORT void read_bin_stl(string fname, vector<triangle> &v);

#endif /* !STLPARSER_H */
