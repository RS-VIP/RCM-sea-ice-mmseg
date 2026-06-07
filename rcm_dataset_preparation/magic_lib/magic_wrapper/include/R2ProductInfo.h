//
// Created by max on 2021-11-03.
//

#ifndef MAGIC_PORT_R2PRODUCTINFO_H
#define MAGIC_PORT_R2PRODUCTINFO_H

#include "r2xml.h"
//#include "r2tiff.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>

namespace py=pybind11;

py::dict GetProductInfo(std::string sProductPath, int blockaverage=1)
{
    R2ProductInfo r2_pinfo;
    std::vector<R2FileInfo> r2_finfos;
    std::vector<GroundControlPoint> r2_gcps;

    double latmin, lonmin, latmax, lonmax;

    R2XML::LoadR2XML(sProductPath, r2_pinfo, r2_finfos, r2_gcps, latmin, lonmin, latmax, lonmax, blockaverage);

    int N = r2_gcps.size();
    py::array_t<float, py::array::c_style> GCP({4,N});
    auto ra = GCP.mutable_unchecked<2>();
    for (int i=0; i<N; i++){
        GroundControlPoint gcp = r2_gcps[i];
        ra(0,i) = gcp.lat;
        ra(0,i) = gcp.lon;
        ra(0,i) = gcp.line;
        ra(0,i) = gcp.pixel;
    }


    py::dict productInfo;
    productInfo["GroundControlPoints"] = GCP;
    productInfo["PassDirection"] = r2_pinfo.passdirection;

    return productInfo;
}


py::dict GetCalibrationLut(std::string sLutPath){

    double offset;
    std::vector<double> gains;

    R2XML::LoadR2LUT(sLutPath, offset, gains);

    py::dict LUT;
    LUT["gains"] = gains;
    LUT["offset"] = offset;

    return LUT;

}



#endif //MAGIC_PORT_R2PRODUCTINFO_H
