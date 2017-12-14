# include <cstdlib>
# include <iostream>
#include <armadillo>
# include <iomanip>
# include <cmath>
#include <stdlib.h>
#include <vector>
#include <dirent.h>


using namespace std;
using namespace arma;
#include "stla_io.h"

//Kmeans start

/*var sse = {};
for (var k = 1; k <= maxK; ++k) {
    sse[k] = 0;
    clusters = kmeans(dataset, k);
    clusters.forEach(function(cluster) {
        mean = clusterMean(cluster);
        cluster.forEach(function(datapoint) {
            sse[k] += Math.pow(datapoint - mean, 2);
        });
    });
}
*/

#define EPSILON 0.000001

double DOT(vector<double> a, vector<double> b)
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

vector<double> CROSS(vector<double> a,vector<double> b)
{
	vector<double> c = {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};

	return c;
}

void SUB(vector<double> & ret, vector <double> a, vector <double> b)
{
	//std::vector<double>::iterator it1;
	for(int i = 0; i < a.size(); ++i)
	{
		ret.push_back(a[i] - b[i]);
	}
}

int triangle_intersection( const vector<double>   V1,  // Triangle vertices
                           const vector<double>   V2,
                           const vector<double>   V3,
                           const vector<double>    O,  //Ray origin
                           const vector<double>    D,  //Ray direction
                                 double* out )
{
  vector<double> e1, e2;  //Edge1, Edge2
  vector<double> P, Q, T;
  double det, inv_det, u, v;
  double t;

  //Find vectors for two edges sharing V1
  SUB(e1, V2, V1);
  SUB(e2, V3, V1);

  //Begin calculating determinant - also used to calculate u parameter
  P = CROSS(D, e2);
  //if determinant is near zero, ray lies in plane of triangle or ray is parallel to plane of triangle
  det = DOT(e1, P);
  //NOT CULLING
  if(det > -EPSILON && det < EPSILON) return 0;
  inv_det = 1.f / det;

  //calculate distance from V1 to ray origin
  SUB(T, O, V1);
  //Calculate u parameter and test bound
  u = DOT(T, P) * inv_det;
  //The intersection lies outside of the triangle
  if(u < 0.f || u > 1.f) return 0;

  //Prepare to test v parameter
  Q = CROSS(T, e1);

  //Calculate V parameter and test bound
  v = DOT(D, Q) * inv_det;
  //The intersection lies outside of the triangle
  if(v < 0.f || u + v  > 1.f) return 0;

  t = DOT(e2, Q) * inv_det;

  if(t > EPSILON) { //ray intersection
    //*out = t;
    return 1;
  }

  // No hit, no win
  return 0;
}

mat vec2mat(vector<vector<double> >&vec){
    int rows = vec.size();
    int cols = vec[0].size();
    mat A(rows, cols);
    for(int i = 0; i<rows; i++){
        for(int j=0; j<cols; j++){
            A(i, j) = vec[i][j];
        }
    }
    return A;
}

double getDistance(rowvec a, rowvec b){
    rowvec temp = a - b;
    return norm(temp, 2);
}

double which_is_nearest(vector<rowvec>& centroids, rowvec pt){
    double minDistance = 0;
    int minLabel = 0;
    for(int i=0; i<centroids.size(); i++){
        double tempDistance = getDistance(centroids[i], pt);
        if(i == 0|| tempDistance < minDistance){
            minDistance = tempDistance;
            minLabel = i;
        }
    }
    return minLabel;
}

double getDistortionFunction(mat data, vector<vector<int> >& cluster, vector<rowvec>& centroids){

    int nSamples = data.n_rows;
    int nDim = data.n_cols;
    double SumDistortion = 0.0;
    for(int i = 0; i < cluster.size(); i++){
        for(int j = 0; j < cluster[i].size(); j++){
            double temp = getDistance(data.row(cluster[i][j]), centroids[i]);
            SumDistortion += temp;
        }
    }
    return SumDistortion;
}

static double uniform(void)
{ int z;
  static const int m1 = 2147483563;
  static const int m2 = 2147483399;
  const double scale = 1.0/m1;

  static int s1 = 0;
  static int s2 = 0;

  if (s1==0 || s2==0) /* initialize */
  { unsigned int initseed = (unsigned int) time(0);
    srand(initseed);
    s1 = rand();
    s2 = rand();
  }

  do
  { int k;
    k = s1/53668;
    s1 = 40014*(s1-k*53668)-k*12211;
    if (s1 < 0) s1+=m1;
    k = s2/52774;
    s2 = 40692*(s2-k*52774)-k*3791;
    if(s2 < 0) s2+=m2;
    z = s1-s2;
    if(z < 1) z+=(m1-1);
  } while (z==m1); /* To avoid returning 1.0 */

  return z*scale;
}


static int binomial(int n, double p)
{
	const double q = 1 - p;
  if (n*p < 30.0) /* Algorithm BINV */
  { const double s = p/q;
    const double a = (n+1)*s;
    double r = exp(n*log(q)); /* pow() causes a crash on AIX */
    int x = 0;
    double u = uniform();
    while(1)
    { if (u < r) return x;
      u-=r;
      x++;
      r *= (a/x)-s;
    }
  }
  else /* Algorithm BTPE */
  { /* Step 0 */
    const double fm = n*p + p;
    const int m = (int) fm;
    const double p1 = floor(2.195*sqrt(n*p*q) -4.6*q) + 0.5;
    const double xm = m + 0.5;
    const double xl = xm - p1;
    const double xr = xm + p1;
    const double c = 0.134 + 20.5/(15.3+m);
    const double a = (fm-xl)/(fm-xl*p);
    const double b = (xr-fm)/(xr*q);
    const double lambdal = a*(1.0+0.5*a);
    const double lambdar = b*(1.0+0.5*b);
    const double p2 = p1*(1+2*c);
    const double p3 = p2 + c/lambdal;
    const double p4 = p3 + c/lambdar;
    while (1)
    { /* Step 1 */
      int y;
      int k;
      double u = uniform();
      double v = uniform();
      u *= p4;
      if (u <= p1) return (int)(xm-p1*v+u);
      /* Step 2 */
      if (u > p2)
      { /* Step 3 */
        if (u > p3)
        { /* Step 4 */
          y = (int)(xr-log(v)/lambdar);
          if (y > n) continue;
          /* Go to step 5 */
          v = v*(u-p3)*lambdar;
        }
        else
        { y = (int)(xl+log(v)/lambdal);
          if (y < 0) continue;
          /* Go to step 5 */
          v = v*(u-p2)*lambdal;
        }
      }
      else
      { const double x = xl + (u-p1)/c;
        v = v*c + 1.0 - fabs(m-x+0.5)/p1;
        if (v > 1) continue;
        /* Go to step 5 */
        y = (int)x;
      }
      /* Step 5 */
      /* Step 5.0 */
      k = abs(y-m);
      if (k > 20 && k < 0.5*n*p*q-1.0)
      { /* Step 5.2 */
        double rho = (k/(n*p*q))*((k*(k/3.0 + 0.625) + 0.1666666666666)/(n*p*q)+0.5);
        double t = -k*k/(2*n*p*q);
        double A = log(v);
        if (A < t-rho) return y;
        else if (A > t+rho) continue;
        else
        { /* Step 5.3 */
          double x1 = y+1;
          double f1 = m+1;
          double z = n+1-m;
          double w = n-y+1;
          double x2 = x1*x1;
          double f2 = f1*f1;
          double z2 = z*z;
          double w2 = w*w;
          if (A > xm * log(f1/x1) + (n-m+0.5)*log(z/w)
                + (y-m)*log(w*p/(x1*q))
                + (13860.-(462.-(132.-(99.-140./f2)/f2)/f2)/f2)/f1/166320.
                + (13860.-(462.-(132.-(99.-140./z2)/z2)/z2)/z2)/z/166320.
                + (13860.-(462.-(132.-(99.-140./x2)/x2)/x2)/x2)/x1/166320.
                + (13860.-(462.-(132.-(99.-140./w2)/w2)/w2)/w2)/w/166320.)
             continue;
          return y;
        }
      }
      else
      { /* Step 5.1 */
        int i;
        const double s = p/q;
        const double aa = s*(n+1);
        double f = 1.0;
        for (i = m; i < y; f *= (aa/(++i)-s));
        for (i = y; i < m; f /= (aa/(++i)-s));
        if (v > f) continue;
        return y;
      }
    }
  }
  /* Never get here */
  return -1;
}

static void randomassign (int nclusters, int nelements, vector<int> &clusterid)
{
	 int i, j;
  int k = 0;
  double p;
  int n = nelements-nclusters;
  /* Draw the number of elements in each cluster from a multinomial
   * distribution, reserving ncluster elements to set independently
   * in order to guarantee that none of the clusters are empty.
   */
  for (i = 0; i < nclusters-1; i++)
  { p = 1.0/(nclusters-i);
    j = binomial(n, p);
    n -= j;
    j += k+1; /* Assign at least one element to cluster i */
    for ( ; k < j; k++) clusterid.push_back(i);
  }
  /* Assign the remaining elements to the last cluster */
  for ( ; k < nelements; k++) clusterid.push_back(i);

  /* Create a random permutation of the cluster assignments */
  for (i = 0; i < nelements; i++)
  { j = (int) (i + (nelements-i)*uniform());
    k = clusterid[j];
    clusterid[j] = clusterid[i];
    clusterid[i] = k;
  }

  return;
}

double kmeans(mat data, int K, vector<vector<int> > &cluster, vector<rowvec> &centroids){

    //mat data = vec2mat(input_vectors);
    int nSamples = data.n_rows;
    int nDim = data.n_cols;
    double lastDistortion = 0.0;
    vector<int> tclusterid;
    randomassign (K, nSamples, tclusterid);
    //randomly select k samples to initialize k centroid
    for(int k = 0; k < K; k ++){
        //int rand_int = rand() % nSamples;
        centroids.push_back(data.row(tclusterid[k]));
    }

    //iteratively find better centroids

    for(int k = 0; k < K; k ++){
        vector<int > vectemp;
        cluster.push_back(vectemp);
    }
    int counter = 0;
    while(1){
        for(int k = 0; k < K; k ++){
            cluster[k].clear();
        }
        bool converge = true;
        //for every sample, find which cluster it belongs to,
        //by comparing the distance between it and each clusters' centroid.
        for(int m = 0; m < nSamples; m++){
            int which = which_is_nearest(centroids, data.row(m));
            cluster[which].push_back(m);
        }
        //for every cluster, re-calculate its centroid.
        for(int k = 0; k < K; k++){
            int cluster_size = cluster[k].size();
            rowvec vectemp = zeros<rowvec>(nDim);
            for(int i = 0; i < cluster_size; i++){
                vectemp = vectemp + data.row(cluster[k][i]) / (double)cluster_size;
            }
            //if(getDistance(centroids[k], vectemp) >= 1e-6) converge = false;
            centroids[k] = vectemp;
        }
        double currentDistortion = getDistortionFunction(data, cluster, centroids);

        if(pow((currentDistortion - lastDistortion), 2) >= 1e-6)
        {
        	converge = false;
        }
        lastDistortion = currentDistortion;
        if(converge) break;
        ++counter;
    }
    double distortion = getDistortionFunction(data, cluster, centroids);
    return distortion;
}

//Kmeans ends




double GetBRandomNumTWZeroOne()
{
	return ((double) rand() / (RAND_MAX));
}

int GetRandomNumber(int rand_max)
{
	return (rand()%rand_max);
}

bool CreateHistogramImproved(int face_num, int node_num, double * node_xyz, int * face_node, double * face_normal, vector <double> &hist)
{
	int randFacet = 0, j = 0;
	double randPoint[2][3];
	vector <double> dist, distIN, distOUT, distMIXED;
	vector <vector<double>> nodeTriangle;
	vector <double> rayDirection;
	rayDirection.resize(3);
	int size = 0.1 * face_num;
	if(size < 1000)
		size = face_num;
	dist.resize(size);
	/*distIN.resize(size);
	distOUT.resize(size);
	distMIXED.resize(size);*/
	hist.resize(3);
	std::vector<double>::iterator it;
	for(it = hist.begin(); it != hist.end(); ++it)
		*it = 0;

	for (int i = 0; i < size; ++i)
	{
		randFacet = GetRandomNumber(face_num);
		int facetTwo = randFacet;
		int facetOne = 0;
		for(j = 0; j < 2 ; ++j)
		{
			int tempRandFacet = randFacet;
			while(tempRandFacet == randFacet && j == 1)
					tempRandFacet = GetRandomNumber(face_num);
			facetOne = tempRandFacet;
			int node1 = face_node[0+tempRandFacet*3];
			int node2 = face_node[1+tempRandFacet*3];
			int node3 = face_node[2+tempRandFacet*3];
			double X1 = node_xyz[0+node1*3];
			double Y1 = node_xyz[1+node1*3];
			double Z1 = node_xyz[2+node1*3];

			double X2 = node_xyz[0+node2*3];
			double Y2 = node_xyz[1+node2*3];
			double Z2 = node_xyz[2+node2*3];

			double X3 = node_xyz[0+node3*3];
			double Y3 = node_xyz[1+node3*3];
			double Z3 = node_xyz[2+node3*3];
			double r1 = GetBRandomNumTWZeroOne();
			double r2 = GetBRandomNumTWZeroOne();
			randPoint[j][0] = (1-sqrt(r1))*X1+sqrt(r2)*(1-r2)*X2+sqrt(r1)*r2*X3;
			randPoint[j][1] = (1-sqrt(r1))*Y1+sqrt(r2)*(1-r2)*Y2+sqrt(r1)*r2*Y3;
			randPoint[j][2] = (1-sqrt(r1))*Z1+sqrt(r2)*(1-r2)*Z2+sqrt(r1)*r2*Z3;


		}
		dist[i] = sqrt(pow((randPoint[0][0] - randPoint[1][0]),2) + pow((randPoint[0][1] - randPoint[1][1]),2) + pow((randPoint[0][2] - randPoint[1][2]),2));
		//direction of line
		rayDirection[0] = randPoint[1][0] - randPoint[0][0];
		rayDirection[1] = randPoint[1][1] - randPoint[0][1];
		rayDirection[2] = randPoint[1][2] - randPoint[0][2];

		double magRayDir = sqrt(pow(rayDirection[0],2) + pow(rayDirection[1],2) + pow(rayDirection[2],2) );

		rayDirection[0] = rayDirection[0]/magRayDir;
		rayDirection[1] = rayDirection[1]/magRayDir;
		rayDirection[2] = rayDirection[2]/magRayDir;

		double n1_x = face_normal[0+facetOne*3];
		double n1_y = face_normal[1+facetOne*3];
		double n1_z = face_normal[2+facetOne*3];

		double n2_x = face_normal[0+facetTwo*3];
		double n2_y = face_normal[1+facetTwo*3];
		double n2_z = face_normal[2+facetTwo*3];
		int numOfFaces = 0;
		for(int k = 0; k < face_num; ++k)
		{
			int node1 = face_node[0+k*3];
			int node2 = face_node[1+k*3];
			int node3 = face_node[2+k*3];
			vector <double> V1, V2, V3;
			V1.resize(3);
			V2.resize(3);
			V3.resize(3);

			V1[0] = node_xyz[0+node1*3];
			V1[1] = node_xyz[1+node1*3];
			V1[2] = node_xyz[2+node1*3];

			V2[0] = node_xyz[0+node2*3];
			V2[1] = node_xyz[1+node2*3];
			V2[2] = node_xyz[2+node2*3];

			V3[0] = node_xyz[0+node3*3];
			V3[1] = node_xyz[1+node3*3];
			V3[2] = node_xyz[2+node3*3];
			vector<double> start = {randPoint[0][0], randPoint[0][1], randPoint[0][2]};
			double * out = NULL;
			int isIntersect = triangle_intersection(V1, V2, V3, start, rayDirection, out);
			if(isIntersect == 1)
			{
				++numOfFaces;
				if(numOfFaces >= 4)
					break;
			}
		}

		double dp = n2_x*rayDirection[0] + n2_y*rayDirection[1] + n2_z*rayDirection[2];
		// decision making, distance belongs to which category
		if(dp <= 0 && numOfFaces < 4)
			//IN distance
			distIN.push_back(dist[i]);
		else if(numOfFaces < 4)
			//OUT distance
			distOUT.push_back(dist[i]);
		else
			//MIXED distance
			distMIXED.push_back(dist[i]);
	}

	double all = dist.size();
	double inPer = distIN.size();
	double outPer = distOUT.size();
	double mixedPer = distMIXED.size();

	hist[0] = (inPer/all)*100;
	hist[1] = (outPer/all)*100;
	hist[2] = (mixedPer/all)*100;

	return true;

}

bool CreateHistogram(int face_num, int node_num, double * node_xyz, int * face_node, vector <double> &hist )
{
	int randFacet = 0, j = 0;
	double randPoint[2][3];
	vector <double> dist;
	dist.resize(1024);
	hist.resize(1024);
	std::vector<double>::iterator it;
	for(it = hist.begin(); it != hist.end(); ++it)
		*it = 0;

	for (int i = 0; i < 1024; ++i)
	{
		randFacet = GetRandomNumber(face_num);
		for(j = 0; j < 2 ; ++j)
		{
			int tempRandFacet = randFacet;
			while(tempRandFacet == randFacet && j == 1)
				tempRandFacet = GetRandomNumber(face_num);
			int node1 = face_node[0+tempRandFacet*3];
			int node2 = face_node[1+tempRandFacet*3];
			int node3 = face_node[2+tempRandFacet*3];
			double X1 = node_xyz[0+node1*3];
			double Y1 = node_xyz[1+node1*3];
			double Z1 = node_xyz[2+node1*3];

			double X2 = node_xyz[0+node2*3];
			double Y2 = node_xyz[1+node2*3];
			double Z2 = node_xyz[2+node2*3];

			double X3 = node_xyz[0+node3*3];
			double Y3 = node_xyz[1+node3*3];
			double Z3 = node_xyz[2+node3*3];
			double r1 = GetBRandomNumTWZeroOne();
			double r2 = GetBRandomNumTWZeroOne();
			randPoint[j][0] = (1-sqrt(r1))*X1+sqrt(r2)*(1-r2)*X2+sqrt(r1)*r2*X3;
			randPoint[j][1] = (1-sqrt(r1))*Y1+sqrt(r2)*(1-r2)*Y2+sqrt(r1)*r2*Y3;
			randPoint[j][2] = (1-sqrt(r1))*Z1+sqrt(r2)*(1-r2)*Z2+sqrt(r1)*r2*Z3;


		}
		dist[i] = sqrt(pow((randPoint[0][0] - randPoint[1][0]),2) + pow((randPoint[0][1] - randPoint[1][1]),2) + pow((randPoint[0][2] - randPoint[1][2]),2));
		//dist = dist/(double)1024;
		//++hist[(int)dist];
		//cout<<dist<<endl;
	}
	vector <double> tempDist = dist;

	std::sort(tempDist.begin(), tempDist.end());

	//Normalise the values between 0 - 1023
	for(int k = 0 ; k < 1024; ++k)
	{
		dist[k] = 1023 * ((dist[k] - tempDist[0])/(tempDist[1023] - tempDist[0]));
	}

	//Create histogram
	for(int t = 0; t < 1024; ++t)
	{
		++hist[(int)dist[t]];
	}
	//std::vector<double>::iterator it1;
	//for(it1 = hist.begin(); it1 != hist.end(); ++it1)
		//	cout<<*it1<<endl;;
	return true;
}

void test03 ( string input_file_name, vector <double> & hist )
{
  bool error;
  int *face_node;
  double *face_normal;
  int face_num;

  int node_num;
  double *node_xyz;
  int solid_num;
  int text_num;

  //cout << "\n";
  //cout << "TEST03\n";
  //cout << "  STLA_READ reads an object in an ASCII STL file.\n";
  //cout << "Input file Name: " << input_file_name << endl;

  stla_size ( input_file_name, &solid_num, &node_num, &face_num, &text_num );

  face_node = new int[3*face_num];
  face_normal = new double[3*face_num];
  node_xyz = new double[3*node_num];

  error = stla_read ( input_file_name, node_num, face_num, node_xyz,
    face_node, face_normal );

  if ( error )
  {
    cout << "\n";
    cout << "  STLA_READ returned error flag.\n";
    return;
  }

  //error = CreateHistogram(face_num, node_num, node_xyz, face_node, hist);
  error = CreateHistogramImproved(face_num, node_num, node_xyz, face_node, face_normal, hist);

  //stla_face_normal_print ( face_num, face_normal );
  /*int randFacet = GetRandomNumber(face_num);

  stla_size_print ( input_file_name, solid_num, node_num, face_num, text_num );

  stla_face_node_print ( face_num, face_node );
  stla_face_normal_print ( face_num, face_normal );
  stla_node_xyz_print ( node_num, node_xyz );*/

  delete [] face_node;
  delete [] face_normal;
  delete [] node_xyz;

  return;
}

vector<string> getAllFilesInDir(const char* dirpath)
{
	DIR *dir;
	struct dirent *ent;
	vector<string> files;
	if ((dir = opendir (dirpath)) != NULL) {
	  /* print all the files and directories within directory */
	  while ((ent = readdir (dir)) != NULL) {
	    //printf ("%s\n", ent->d_name);
		  files.push_back(ent->d_name);
	  }
	  closedir (dir);
	} else {
	  /* could not open directory */
	  perror ("");
	  //return EXIT_FAILURE;
	}
	return files;
}

vector<string> searchSimilar(const char* dirpath, string file, int num_clusters)
{
	//string dirpath = "/media/abhishek/54F649860A9E68BC/Code/nodejs/node-ffi-master/example/smartSearch/stl/";
	vector<string> fileNames = getAllFilesInDir(dirpath);
	mat totalSpace(fileNames.size()-2, 3);
	for(int i = 2; i < fileNames.size();++i)
	{
		vector <double> histogram;
		test03(dirpath+fileNames[i], histogram);
		for(int j=0; j<3; ++j){
			totalSpace(i-2, j) = histogram[j];
			        }
	}

	vector<vector<int> > cluster;
	vector<rowvec> centroids;
	kmeans(totalSpace, num_clusters, cluster, centroids);

	vector <double> histogram1;
	test03(file, histogram1);

	mat data(1, 3);
	for(int j=0; j<3; j++){
					            data(0, j) = histogram1[j];
					        }
	int which = which_is_nearest(centroids, data.row(0));
	vector<string> similarSeq;
	vector<int> cen = cluster[which];
	for(int k = 0; k < cluster[which].size(); ++k)
	{
		similarSeq.push_back(fileNames[cluster[which][k]+2]);
	}

	return similarSeq;
}


int main(int argc, char *argv[])
{

	//cout << "File: " << argv[1]<<endl;
	int number = atoi(argv[3]);
	vector<string> similar = searchSimilar(argv[1], argv[2], number);
	string dirpath = "http://localhost:3000/";
	for(int i = 0;i < similar.size(); ++i)
	{
		cout <<"'"<< similar[i];
		if(i != similar.size()- 1)
		{
			cout<<"',";
		}
		else
		{
			cout<<"'";
		}
	}
	/*srand((unsigned)time(0));
	vector <double> histogram1;
	test03("humanoid.stl", histogram1);
	vector <double> histogram2;
	test03("humanoid.stl", histogram2);
	vector <double> histogram3;
	test03("sphere.stl", histogram3);
	vector <double> histogram4;
	test03("sphere.stl", histogram4);
	//vector <double> histogram5;
	//test03("sphere.stl", histogram5);
	//vector <double> histogram6;
	//test03("ULA_Support_Bracket_VB05-ascii.stl", histogram6);

	mat A(4, 3);
	for(int j=0; j<3; j++){
	            A(0, j) = histogram1[j];
	        }

	for(int j=0; j<3; j++){
		            A(1, j) = histogram2[j];
		        }
	for(int j=0; j<3; j++){
		            A(2, j) = histogram3[j];
		        }

	for(int j=0; j<3; j++){
			            A(3, j) = histogram4[j];
			        }

	//cout<<"Before: "<<endl;
	//A.print();
	vector<vector<int> > cluster;
	int num_clusters = 2;
	vector<rowvec> centroids;
	kmeans(A, num_clusters, cluster, centroids);
	cout<<"Cluster Sequence: "<<endl;
	//A.print();
	for(int i = 0; i < num_clusters ; ++i)
	{
		std::vector<int>::iterator it1;
		for(it1 = cluster[i].begin(); it1 != cluster[i].end(); ++it1)
				cout<<*it1<<endl;
		cout<<"next:"<<endl;
	}

	/*vector <double> histogram6;
	test03("ULA_Support_Bracket_VB05-ascii.stl", histogram6);
	mat data(1, 1024);
	for(int j=0; j<1024; j++){
					            data(0, j) = histogram6[j];
					        }
	int which = which_is_nearest(centroids, data.row(0));
	cluster[which].push_back(5);

	cout<<"New cluster Sequence: "<<endl;
		//A.print();
	for(int i = 0; i < num_clusters ; ++i)
	{
		std::vector<int>::iterator it1;
		for(it1 = cluster[i].begin(); it1 != cluster[i].end(); ++it1)
				cout<<*it1<<endl;
		cout<<"next:"<<endl;
	}*/
	return 0;
}
