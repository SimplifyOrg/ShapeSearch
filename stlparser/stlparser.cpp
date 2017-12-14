/* Source file stlparser.cpp */


#include "stlparser.h"
#include <fstream>

v3::v3()
{
	m_x = 0.0;
	m_y = 0.0;
	m_z = 0.0;
}

v3::v3(char* facet)
{
    char f1[4] = {facet[0],
        facet[1],facet[2],facet[3]};

    char f2[4] = {facet[4],
        facet[5],facet[6],facet[7]};

    char f3[4] = {facet[8],
        facet[9],facet[10],facet[11]};

    float xx = *((float*) f1 );
    float yy = *((float*) f2 );
    float zz = *((float*) f3 );

    m_x = double(xx);
    m_y = double(yy);
    m_z = double(zz);
}

void read_bin_stl(string fname, vector<triangle> &v)
{
		char * cstr = new char [fname.length()+1];
		std::strcpy (cstr, fname.c_str());
    //!!
    //don't forget ios::binary
    //!!
    ifstream stlFile (cstr, ios::in | ios::binary);

    char header_info[80] = "";
    char nTri[4];
    //unsigned long nTriLong;
		char node[12];

    //read 80 byte header
    if (stlFile)
		{
        stlFile.read (header_info, 80);
				stlFile.read (nTri, 4);
				float num_facets = *((float*) nTri );

				for(int i = 0; i < num_facets; ++i)
				{
					stlFile.read (node, 12);
					v3 *v1 = new v3(node);
					stlFile.read (node, 12);
					v3 *v2 = new v3(node);
					stlFile.read (node, 12);
					v3 *v31 = new v3(node);
					stlFile.read (node, 12);
					v3 *v4 = new v3(node);
					triangle *t = new triangle(*v2,*v31,*v4,*v1);
					v.push_back(*t);
					char att[2];
					stlFile.read (att, 2);
				}
		}
}
