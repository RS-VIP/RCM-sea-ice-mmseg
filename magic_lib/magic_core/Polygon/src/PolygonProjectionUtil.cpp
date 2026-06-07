#include "PolygonProjectionUtil.h"
#include <fstream>

//vector<GroundControlPoint> PolygonProjectionUtil::imagetiepoints;

/// Lambert Conic Conformal Projection Parameters
LCCProjectionInfo::LCCProjectionInfo() : a(6378206.400),
                                         f(1.0/294.97870),
                                         e(0.08227185),
                                         phi_1(49.0*DEGS_PER_RAD),
                                         phi_2(77.0*DEGS_PER_RAD),
                                         f_phi(40.0*DEGS_PER_RAD),
                                         f_lambda(-100.0*DEGS_PER_RAD),
                                         f_easting(0.0),
                                         f_northing(0.0),
                                         ulx(1),
                                         uly(1),
                                         size_x(1),
                                         size_y(1)
{
}

/// Initializes internal values, taking parameters in radians (params will be stored in radians)
void LCCProjectionInfo::InitializeValuesRadians(double pa,
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
                                                float psizey)
{
    a = pa; f = pf; e = pe;
    phi_1 = pphi_1;
    phi_2 = pphi_2;
    f_phi = pf_phi;
    f_lambda = pf_lambda;
    f_easting = pf_easting;
    f_northing = pf_northing;
    ulx = pulx;
    uly = puly;
    size_x = psizex;
    size_y = psizey;
}

/// Initializes internal values, taking parameters in degrees (params will be stored in radians)
void LCCProjectionInfo::InitializeValuesDeg(double pa,
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
                                            float psizey)
{
    a = pa; f = pf; e = pe;
    phi_1 = pphi_1_deg * DEGS_PER_RAD;
    phi_2 = pphi_2_deg * DEGS_PER_RAD;
    f_phi = pf_phi_deg * DEGS_PER_RAD;
    f_lambda = pf_lambda_deg * DEGS_PER_RAD;
    f_easting = pf_easting;
    f_northing = pf_northing;
    ulx = pulx;
    uly = puly;
    size_x = psizex;
    size_y = psizey;
}

/// Converts latitude / longitude to image coordinates according to the ground control points
PolygonProjectionUtil::PolygonProjectionUtil()
{
    //imagetiepoints.clear();
    InitializeTransformations();
}

void PolygonProjectionUtil::SetLCCProjInfo(const LCCProjectionInfo& pinfo)
{
    projinfo = pinfo;
}

void PolygonProjectionUtil::SetGCP(const vector<GroundControlPoint>& itp)
{
    imagetiepoints.clear();
    imagetiepoints = itp;
    mineast = 1000000000.0;
    minnorth = 1000000000.0;
    x_coef.resize(6);
    y_coef.resize(6);
    r_coef.resize(6);
    c_coef.resize(6);
    InitializeTransformations();
}


/// Converts lambda and phi (longitude and latitude in radians, respectively) to LCC coordinates based on projection info in pinfo.
void PolygonProjectionUtil::lcc(double &p_easting, double &p_northing, const double lambda,
                                const double phi, const LCCProjectionInfo& pinfo)
{

    double m1, m2, t, t1, t2, tf, e2, pi, n, F, r, rf, theta;
    e2 = pow(pinfo.e,2.0); pi = 3.1415926535897932;

    m1 = cos(pinfo.phi_1)/pow(1-(e2*pow(sin(pinfo.phi_1),2)),0.5);
    m2 = cos(pinfo.phi_2)/pow(1-(e2*pow(sin(pinfo.phi_2),2)),0.5);
    t = (tan((pi/4)-(phi/2))/pow((1-(pinfo.e*sin(phi)))/(1+(pinfo.e*sin(phi))),(pinfo.e/2)));
    t1 = tan((pi/4)-(pinfo.phi_1/2))/pow((1-(pinfo.e*sin(pinfo.phi_1)))/(1+(pinfo.e*sin(pinfo.phi_1))),(pinfo.e/2));
    t2 = tan((pi/4)-(pinfo.phi_2/2))/pow((1-(pinfo.e*sin(pinfo.phi_2)))/(1+(pinfo.e*sin(pinfo.phi_2))),(pinfo.e/2));
    tf = tan((pi/4)-(pinfo.f_phi/2))/pow((1-(pinfo.e*sin(pinfo.f_phi)))/(1+(pinfo.e*sin(pinfo.f_phi))),(pinfo.e/2));
    n = (log(m1)-log(m2))/(log(t1) - log(t2));
    F = m1/(n*pow(t1,n));
    r = pinfo.a*F*pow(t,n);
    rf = pinfo.a*F*pow(tf,n);
    theta = n*(lambda-pinfo.f_lambda);
    p_northing = pinfo.f_northing+rf-(r*cos(theta));
    p_easting = pinfo.f_easting+(r*sin(theta));
}

/// Converts the polygon defined by the given list of lat / lon coordinates (in degrees) to image coordinates.
vector<POINTex> PolygonProjectionUtil::PolygonLatLongToImageCoords(vector<double> lat, vector<double> lon)
{
    vector<POINTex> ret;

    if (x_coef.size() < 6 || y_coef.size() < 6)
        return ret;

    ret.resize(lat.size());

    POINTex point;

    double phi;
    double lambda;
    double easting;
    double northing;

    // Convert the points to image coordinates.
    for (size_t j=0; j < lat.size(); j++)
    {
        phi =  lat[j] * DEGS_PER_RAD;
        lambda =  lon[j] * DEGS_PER_RAD;

        lcc(easting, northing, lambda, phi, projinfo);
        easting = (float) easting - (float) mineast;
        northing = (float) northing - (float) minnorth;

        point.Col = (int) floor((x_coef[0] + (x_coef[1] * easting) + (x_coef[2] * northing) + (x_coef[3] * easting * northing) + (x_coef[4] * easting * easting) + (x_coef[5] * northing * northing))+0.5);
        point.Row = (int) floor((y_coef[0] + (y_coef[1] * easting) + (y_coef[2] * northing) + (y_coef[3] * easting * northing) + (y_coef[4] * easting * easting) + (y_coef[5] * northing * northing))+0.5);

        // add to points list.
        ret[j] = point;
    }
    return ret;
}

//Ross Feb 2016
// Converts the polygon defined by the given list of row/col coordinates to lat/lon (in degrees)

POINTex PolygonProjectionUtil::PolygonImageToLatLongCoords(vector<float> row, vector<float> col, double* lat, double* lon)
{
    POINTex ret;
    ret.Row = 0;
    ret.Col = 0;
    if (r_coef.size() < 6 || c_coef.size() < 6)
        return ret;

    double easting, northing;

    for (size_t j = 0; j < (size_t)row.size(); j++)
    {
        //Convert from image to LCC
        easting  = c_coef[0] + (c_coef[1] * col[j]) + (c_coef[2] * row[j]) + (c_coef[3] * col[j] * row[j]) + (c_coef[4] * col[j] * col[j]) + (c_coef[5] * row[j] * row[j]);
        northing = r_coef[0] + (r_coef[1] * col[j]) + (r_coef[2] * row[j]) + (r_coef[3] * col[j] * row[j]) + (r_coef[4] * col[j] * col[j]) + (r_coef[5] * row[j] * row[j]);

        easting += mineast;
        northing += minnorth;
        //Use to save easting/northing of a point
        if (j == 0)
        {
            ret.Col = (int)easting;
            ret.Row = (int)northing;
        }

        //Convert northing and easting to latitude and longitude
        PolygonLCCToLatLong(easting, northing, lat[j], lon[j]);
        //round lat and lon to 5 decimal places
        lat[j] = floor(lat[j]*100000 + 0.5)/100000;
        lon[j] = floor(lon[j]*100000 + 0.5)/100000;
    }
    return ret;
}

// Uses Lambert conformal projection to geographic from
//http://www.linz.govt.nz/data/geodetic-system/coordinate-conversion/projection-conversions/lambert-conformal-conic-geographic

void PolygonProjectionUtil::PolygonLCCToLatLong(const double easting, const double northing, double &lat, double &lon)
{
    const LCCProjectionInfo& pinfo = projinfo;
    double x_tilda = easting-pinfo.f_easting;
    double y_tilda = northing-pinfo.f_northing;
    double m1, m2, t_tilda, t1, t2, tf, e2, pi, n, F, rho, rf, theta, sgn_n;

    pi = 3.1415926535897932;
    e2 = pow(pinfo.e,2.0);
    m1 = cos(pinfo.phi_1)/pow(1-(e2*pow(sin(pinfo.phi_1),2)),0.5);
    m2 = cos(pinfo.phi_2)/pow(1-(e2*pow(sin(pinfo.phi_2),2)),0.5);
    t1 = tan((pi/4)-(pinfo.phi_1/2))/pow((1-(pinfo.e*sin(pinfo.phi_1)))/(1+(pinfo.e*sin(pinfo.phi_1))),(pinfo.e/2));
    t2 = tan((pi/4)-(pinfo.phi_2/2))/pow((1-(pinfo.e*sin(pinfo.phi_2)))/(1+(pinfo.e*sin(pinfo.phi_2))),(pinfo.e/2));
    n = (log(m1)-log(m2))/(log(t1) - log(t2));
    F = m1/(n*pow(t1,n));
    tf = tan((pi/4)-(pinfo.f_phi/2))/pow((1-(pinfo.e*sin(pinfo.f_phi)))/(1+(pinfo.e*sin(pinfo.f_phi))),(pinfo.e/2));
    rf = pinfo.a*F*pow(tf,n);

    if(n>=0)
        sgn_n = 1.0;
    else if(n<0)
        sgn_n = -1.0;
    else
        sgn_n = 0.0;

    rho = sgn_n*pow((pow(x_tilda,2)+pow(rf-y_tilda,2)),0.5);
    t_tilda = pow(rho/(pinfo.a*F),1/n);
    theta = atan(x_tilda/(rf-y_tilda));
    lat = 0.5*pi-2*atan(t_tilda);
    double lat_old;
    int count = 0;
    do
    {
        lat_old = lat;
        lat = 0.5*pi - 2*atan(t_tilda*(pow((1-pinfo.e*sin(lat))/(1+pinfo.e*sin(lat)),0.5*pinfo.e)));
        count++;
    } while (abs(lat - lat_old) > 1.0e-8 && count < 10);
    lat=lat*180/pi;
    lon = (pinfo.f_lambda + theta/n)*180/pi;
    return;
}

/// Converts the polygon defined by the given list of lat / lon coordinates (in degrees) to LCC coordinates.
vector<POINTex> PolygonProjectionUtil::PolygonLatLongToLCC(vector<double> lat, vector<double> lon)
{
    vector<POINTex> ret(lat.size());
    POINTex point;

    double phi;
    double lambda;
    double easting;
    double northing;

    // Convert the points to LCC
    for (size_t j=0; j < lat.size(); j++)
    {
        phi =  lat[j] * DEGS_PER_RAD;
        lambda =  lon[j] * DEGS_PER_RAD;

        lcc(easting, northing, lambda, phi, projinfo);
        point.Row =  (int) floor((double)(((float)northing-projinfo.uly)/(projinfo.size_y))+0.5);
        point.Col = (int) floor((double)(((float)easting-projinfo.ulx)/(projinfo.size_x))+0.5);

        // add to points list.
        ret[j] = point;
    }
    return ret;
}

bool PolygonProjectionUtil::isInitialized()
{
    if (imagetiepoints.empty())
        return false;
    return true;
}

PolygonProjectionUtil::~PolygonProjectionUtil()
{

}

/// Initialize transform coefficients based on projection info.
void PolygonProjectionUtil::InitializeTransformations()
{
    // allocate arrays for the GCPs.
    vector<double> east_gcp(imagetiepoints.size());
    vector<double> north_gcp(imagetiepoints.size());
    vector<double> xgcp(imagetiepoints.size());
    vector<double> ygcp(imagetiepoints.size());

    // LCC coordinates for each GCP.
    double easting = 0.0;
    double northing = 0.0;

    double phi;
    double lambda;

    // Convert each GCP to LCC coordinates.
    for(size_t i = 0; i < imagetiepoints.size(); i++)
    {
        phi =  imagetiepoints[i].lat * DEGS_PER_RAD;
        lambda =  imagetiepoints[i].lon * DEGS_PER_RAD;

        lcc(easting, northing, lambda, phi, projinfo);
        if(easting < mineast) mineast = easting;
        if(northing < minnorth) minnorth = northing;

        east_gcp[i] = easting;
        north_gcp[i] = northing;
        xgcp[i] = (double) imagetiepoints[i].pixel;
        ygcp[i] = (double) imagetiepoints[i].line;

    }

    // Subtract the min from the LCC coordinates.
    for(size_t i=0; i < imagetiepoints.size(); i++)
    {
        east_gcp[i] = east_gcp[i] - mineast;
        north_gcp[i] = north_gcp[i] - minnorth;

    }


    // calculate the fitting coefficients that converts LCC to image coordinates.
    get_pol2_coef(x_coef, east_gcp, north_gcp, xgcp);  // coefficients for the columns
    get_pol2_coef(y_coef, east_gcp, north_gcp, ygcp); // coefficients for the lines

    // calculate the fitting coefficients that converts image to LCC cooridnates
    get_pol2_coef(c_coef, xgcp, ygcp, east_gcp); //coefficients for longitude
    get_pol2_coef(r_coef, xgcp, ygcp, north_gcp); //coefficients for latitude
}


/// Fits the GCPs to get relationship between image and geographic coordinates.
/** Returns the coefs in a.
 * x are the easting values of each GCP.
 * y are the northing values of each GCP.
 * colgcp are either the line or column image coordinate of each GCP
 */
void PolygonProjectionUtil::get_pol2_coef(vector<double>& a, const vector<double>& x, const vector<double>& y, const vector<double>& colgcp)
{
    double	n, *z, sx, sy, sxy, sx2, sy2, sx2y, sxy2, sx3, sy3, sx2y2, sx3y, sxy3, sx4, sy4;
    double	sz, sxz, syz, sxyz, sx2z, sy2z, **xtx, xtz[6];
    long		i, j;

    int n_gcp = (int) x.size();

    a.resize(6);

    z = (double *) calloc(n_gcp, sizeof(double));
    sz = (double) 0.0; sx = (double) 0.0; sy = (double) 0.0; sxy = (double) 0.0; sx2 = (double) 0.0; sy2 = (double) 0.0;
    sx2y = (double) 0.0; sxy2 = (double) 0.0; sx3 = (double) 0.0; sy3 = (double) 0.0; sx2y2 = (double) 0.0;
    sx3y = (double) 0.0; sxy3 = (double) 0.0; sx4 = (double) 0.0; sy4 = (double) 0.0;
    sz = (double) 0.0; sxz = (double) 0.0; syz = (double) 0.0; sxyz = (double) 0.0; sx2z = (double) 0.0; sy2z = (double) 0.0;
    n = (double) n_gcp;
    xtx = (double **) calloc(6, sizeof(double *));

    for(i=0; i<6; i++) xtx[i] = (double *) calloc(6, sizeof(double));

    for(i=0; i<n_gcp; i++)
    {
        z[i] = (double) colgcp[i];
        sx = sx + x[i]; sy = sy + y[i];
        sxy = sxy + (x[i] * y[i]); sx2y2 = sx2y2 + (x[i] * x[i] * y[i] * y[i]);
        sx2 = sx2 + (x[i] * x[i]); sy2 = sy2 + (y[i] * y[i]);
        sx2y = sx2y + (x[i] * x[i] * y[i]); sxy2 = sxy2 + (x[i] * y[i] * y[i]);
        sx3 = sx3 + (x[i] * x[i] * x[i]); sy3 = sy3 + (y[i] * y[i] * y[i]);
        sx3y = sx3y + (x[i] * x[i] * x[i] * y[i]); sxy3 = sxy3 + (x[i] * y[i] * y[i] * y[i]);
        sx4 = sx4 + (x[i] * x[i] * x[i] * x[i]); sy4 = sy4 + (y[i] * y[i] * y[i] * y[i]);
        sz = sz + z[i]; sxz = sxz + (x[i] * z[i]); syz = syz + (y[i] * z[i]);
        sxyz = sxyz + (x[i] * y[i] * z[i]); sx2z = sx2z + (x[i] * x[i] * z[i]); sy2z = sy2z + (y[i] * y[i] * z[i]);

    }

    xtx[0][0] = n;		xtx[0][1] = sx;		xtx[0][2] = sy;		xtx[0][3] = sxy;	xtx[0][4] = sx2;	xtx[0][5] = sy2;
    xtx[1][0] = sx;		xtx[1][1] = sx2;	xtx[1][2] = sxy;	xtx[1][3] = sx2y;	xtx[1][4] = sx3;	xtx[1][5] = sxy2;
    xtx[2][0] = sy;		xtx[2][1] = sxy;	xtx[2][2] = sy2;	xtx[2][3] = sxy2;	xtx[2][4] = sx2y;	xtx[2][5] = sy3;
    xtx[3][0] = sxy;	xtx[3][1] = sx2y;	xtx[3][2] = sxy2;	xtx[3][3] = sx2y2;	xtx[3][4] = sx3y;	xtx[3][5] = sxy3;
    xtx[4][0] = sx2;	xtx[4][1] = sx3;	xtx[4][2] = sx2y;	xtx[4][3] = sx3y;	xtx[4][4] = sx4;	xtx[4][5] = sx2y2;
    xtx[5][0] = sy2;	xtx[5][1] = sxy2;	xtx[5][2] = sy3;	xtx[5][3] = sxy3;	xtx[5][4] = sx2y2;	xtx[5][5] = sy4;
    xtz[0] = sz;		xtz[1] = sxz;		xtz[2] = syz;		xtz[3] = sxyz;		xtz[4] = sx2z;		xtz[5] = sy2z;

    matinv(xtx, 6);

    for(i=0; i<6; i++)
    {
        a[i] = (double) 0.0;
        for(j=0; j<6; j++)
        {
            a[i] = a[i] + (xtx[i][j] * xtz[j]);
        }
    }

    free(z);
    for(i=0; i<6; i++)
        free(xtx[i]);
    free(xtx);

}

/// Not sure what this function does, included in case its ever needed again.
void PolygonProjectionUtil::get_pol1_coef(vector<double> &a, const vector<double>& x1, const vector<double>& x2, const vector<double>& colgcp)
{
    double	n, *y, sx1, sx2, sy, ssqx1, ssqx2, ssqy, sx1y, sx2y, sx1x2, mx1, mx2, my;
    double	**xtx, ixtx[2][2], xty[2], detxtx;
    int		i;

    int n_gcp = (int) x1.size();

    a.resize(6);
    for(i=0; i<6; i++) a[i] = (double) 0.0;

    y = (double *) calloc(n_gcp, sizeof(double));
    sy = (double) 0.0; sx1 = (double) 0.0; sx2 = (double) 0.0;
    ssqy = (double) 0.0; ssqx1 = (double) 0.0; ssqx2 = (double) 0.0;
    sx1y = (double) 0.0; sx2y = (double) 0.0; sx1x2 = (double) 0.0;
    n = (double) n_gcp;
    xtx = (double **) calloc(6, sizeof(double *));
    if (!(xtx)) {printf("Not enough memory to allocate mask buffer\n"); exit(1);}
    for(i=0; i<6; i++) xtx[i] = (double *) calloc(6, sizeof(double));

    for(i=0; i<n_gcp; i++)
    {
        y[i] = (double) colgcp[i];
        sx1 = sx1 + x1[i];
        sx2 = sx2 + x2[i];
        sy = sy + y[i];
        sx1y = sx1y + (x1[i] * y[i]);
        sx2y = sx2y + (x2[i] * y[i]);
        sx1x2 = sx1x2 + (x1[i] * x2[i]);
        ssqx1 = ssqx1 + (x1[i] *x1[i]);
        ssqx2 = ssqx2 + (x2[i] *x2[i]);
        ssqy = ssqy + (y[i] *y[i]);
    }
    mx1 = sx1 / n; mx2 = sx2 / n; my = sy / n;
    xtx[0][0] = ssqx1 - ((sx1 * sx1)/n); xtx[1][1] = ssqx2 - ((sx2 * sx2)/n);
    xtx[0][1] = sx1x2 - ((sx1 * sx2)/n); xtx[1][0] = xtx[0][1];
    xty[0] = sx1y-((sx1 * sy)/n); xty[1] = sx2y-((sx2 * sy)/n);
    detxtx = determinant(xtx, 2);
    ixtx[0][0] = xtx[1][1]/detxtx; ixtx[1][1] = xtx[0][0]/detxtx;
    ixtx[0][1] = xtx[0][1]*(-1)/detxtx; ixtx[1][0] = xtx[1][0]*(-1)/detxtx;
    a[1] = (ixtx[0][0] * xty[0]) + (ixtx[0][1] * xty[1]);
    a[2] = (ixtx[1][0] * xty[0]) + (ixtx[1][1] * xty[1]);
    a[0] = my - ((a[1] * mx1) + (a[2] * mx2));

}

/// determinant of matrix.
double PolygonProjectionUtil::determinant(double **mat, int n)
{
    int		i, j, k, l, p;
    double	pos, neg, temp1, temp2, **submat, det;

    if(n == 2)
    {
        return((mat[0][0] * mat[1][1]) - (mat[0][1] * mat[1][0]));
    }
    if(n == 3)
    {
        pos = 0.0; neg = 0.0;
        for(i=0; i<n; i++)
        {
            temp1 = 1.0; temp2 = 1.0;
            for(j=0; j<n; j++)
            {
                k = j+i; if(k>n-1) k = k - n;
                temp1 = temp1 * mat[j][k];
                l = i-(j); if(l<0) l = l + n;
                temp2 = temp2 * mat[j][l];
            }
            pos = pos + temp1;
            neg = neg + temp2;
        }
        return(pos-neg);
    }
    if(n > 3)
    {
        det = 0.0;
        submat = (double **) calloc(n-1, sizeof(double *));
        for(i=0; i<n-1; i++) submat[i] = (double *) calloc(n-1, sizeof(double));

        for(j=0; j<n; j++)
        {
            for(k=1; k<n; k++)
            {
                p = 0;
                for(l=0; l<n; l++)
                {
                    if(l != j) {submat[k-1][p] = mat[k][l]; p++;}
                }
            }
            det = det + (mat[0][j] * determinant(submat,n-1) * pow((double) (-1), (double) (j)));
        }

        for(i=0; i<n-1; i++)
            free(submat[i]);
        free(submat);

        return(det);
    }
    else
    {
        return(1.0);
    }
}

/// Inverse of matrix.
void PolygonProjectionUtil::matinv(double **mat, int s)
{
    double	**tmat, **imat, det;
    int		i, j, k, l, m, n;


    imat = (double **) calloc(s, sizeof(double *));
    for(i=0; i<s; i++) imat[i] = (double *) calloc(s, sizeof(double));

    tmat = (double **) calloc(s, sizeof(double *));
    for(i=0; i<s; i++) tmat[i] = (double *) calloc(s, sizeof(double));

    det = determinant(mat,s);

    for(i=0; i<s; i++)
    {
        for(j=0; j<s; j++)
        {
            m = 0;
            for(k=0; k<s; k++)
            {
                n = 0;
                for(l=0; l<s; l++)
                {
                    if(k != i && l != j) {tmat[m][n] = mat[k][l]; n++;}
                }
                if(k != i) m++;
            }
            imat[i][j] = determinant(tmat,s-1);
        }
    }
    for(i=0; i<s; i++)
    {
        for(j=0; j<s; j++)
        {
            mat[j][i] = pow(-1.0,(double) (i+j)) * imat[i][j]/det;
        }
    }

    for(i=0; i<s; i++)
        free(imat[i]);
    free(imat);

    for(i=0; i<s; i++)
        free(tmat[i]);
    free(tmat);

    return;
}
