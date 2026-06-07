/** \file
 *
 * \brief PolygonProjectionUtil class
 *
 * \author Peter Yu
 *
 */

#ifndef POLYGONPROJECTIONUTIL_H
#define POLYGONPROJECTIONUTIL_H

#define DEGS_PER_RAD 0.0174532925199

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

#include "DataType.h"

#define DEGS_PER_RAD 0.0174532925199

using namespace std;

/// Represents a ground control point.
struct GroundControlPoint
{
	//changed by Ross from ints
	float pixel;
	float line;
	double lat; // latitude in DEGREES
	double lon; // longitude in DEGREES
};

/// Lambert Conic Conformal Projection Parameters
class LCCProjectionInfo
{
public:
	
	LCCProjectionInfo();

	/// Initializes internal values, taking parameters in radians (params will be stored in radians)
	void InitializeValuesRadians(double pa,
			double pf,
			double pe,
			double pphi_1,
			double pphi_2,
			double pf_phi,
			double pf_lambda,
			double pf_easting,
			double pf_northing,
			float pulx,
			float puly,
			float psizex,
			float psizey);

	/// Initializes internal values, taking parameters in degrees (params will be stored in radians)
	void InitializeValuesDeg(double pa,
			double pf,
			double pe,
			double pphi_1_deg,
			double pphi_2_deg,
			double pf_phi_deg,
			double pf_lambda_deg,
			double pf_easting,
			double pf_northing,
			float pulx,
			float puly,
			float psizex,
			float psizey);

public:
	/// semi-major axis
	double a;
	/// flattening
	double f;
	/// eccentricity
	double e;
	/// first parallel in radians
	double phi_1;
	/// second parallel in radians
	double phi_2;
	/// false latitude in radians
	double f_phi;
	/// false longtitude in radians
	double f_lambda;
	/// false easting
	double f_easting;
	/// false northing
	double f_northing;
	/// upper left corner x coord
	float ulx;
	/// upper left corner y coord
	float uly;
	/// size of pixel in x.
	float size_x;	
	/// size of pixel in y.
	float size_y;	
};

/// Converts latitude / longitude to image coordinates according to the ground control points
class PolygonProjectionUtil
{
public:

	PolygonProjectionUtil();

	/// Sets the projection info data structure to use.
	void SetLCCProjInfo(const LCCProjectionInfo& pinfo);

	/// sets the GCP points to use.
	void SetGCP(const vector<GroundControlPoint>& itp);

	/// Converts lambda and phi (longitude and latitude in radians, respectively) to LCC coordinates based on projection info in pinfo.
	static void lcc(double &p_easting, double &p_northing, const double lambda, 
			const double phi, const LCCProjectionInfo& pinfo);

	/// Converts the polygon defined by the given list of lat / lon coordinates (in degrees) to image coordinates.
	vector<POINTex> PolygonLatLongToImageCoords(vector<double> lat, vector<double> lon);

	// Converts the polygon defined by the given list of row/col coordinates to lat/lon (in degrees).
	POINTex PolygonImageToLatLongCoords(const vector<float> row, const vector<float> col, double* lat, double* lon);

	/// Converts the polygon defined by the given list of lat / lon coordinates (in degrees) to LCC coordinates.
	vector<POINTex> PolygonLatLongToLCC(vector<double> lat, vector<double> lon);

	/// Converts the polygon defined by the given LCC coordinates to /// Converts the polygon defined by the given list of lat / lon coordinates (in degrees).
	void PolygonLCCToLatLong(const double easting, const double northing, double &lat, double &lon);

	//returns true if imagetiepoints is initialized
	bool isInitialized();

	~PolygonProjectionUtil();


protected:

	/// Initialize transform coefficients based on projection info.
	void InitializeTransformations();

	/// Fits the GCPs to get relationship between image and geographic coordinates.
	/** Returns the coefs in a. 
	 * x are the easting values of each GCP.
	 * y are the northing values of each GCP.
	 * colgcp are either the line or column image coordinate of each GCP
	 */
	void get_pol2_coef(vector<double>& a, const vector<double>& x, const vector<double>& y, const vector<double>& colgcp);

	/// Not sure what this function does, included in case its ever needed again.
	void get_pol1_coef(vector<double> &a, const vector<double>& x1, const vector<double>& x2, const vector<double>& colgcp);

	/// determinant of matrix.
	static double determinant(double **mat, int n);

	/// Inverse of matrix.
	static void matinv(double **mat, int s);

protected:

	/// Coefficients that converts LCC coordinates to image coordinates.
	vector<double> x_coef;	//Found by fitting the GCPs.
	/// Coefficients that converts LCC coordinates to image coordinates.
	vector<double> y_coef;	//Found by fitting the GCPs.

	//Ross Feb 2016
	/// Coefficients that converts image coordinates to LLC coordinates.
	vector<double> r_coef, c_coef; //Found by fitting the GCPs.

	/// LCC Parameters
	LCCProjectionInfo projinfo;

	/// GCP list.
	vector<GroundControlPoint> imagetiepoints;

	/// LCC parameter that must be found from the GCPs
	/**
	 * set to correct value when InitializeTransformations is called.
	 */
	double mineast;

	/// LCC parameter that must be found from the GCPs
	/**
	 * set to correct value when InitializeTransformations is called.
	 */
	double minnorth;
};

#endif
