/** \file
 *
 * \author Peter Yu
 * \ modified by Max Manning to use the tinyXML2 library in place of the windows API for linux compatibility
 *
 * \brief Extracts useful information from the Radarsat-2 XML files, such as GCPs and image parameters.
 *
 */

#pragma once

#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include "PolygonProjectionUtil.h"
#include "tinyxml2.h"
#include <exception>

using namespace std;
using namespace tinyxml2;

/// Stores information about the R2 product.
struct R2ProductInfo
{	string starttime; // Start time of the product.
    string productid; // product ID.
    string producttype; // Product type.
    string passdirection; // ascending or descending.
    int numberoflines; // lines in image(s)
    int numberofsamples; // # of samples in image.
};

/// Stores filename and pole information for each image file in the R2 product.
struct R2FileInfo
{
    string filename;
    string pole;
};

static void readr2product(const string sproductxmlpath)
{
    // Load the information in the XML document.
    XMLDocument xmldoc;
    xmldoc.LoadFile(sproductxmlpath.c_str());

    XMLElement* pRoot = xmldoc.FirstChildElement("product");

    XMLElement* pImgAttr = pRoot->FirstChildElement("imageAttributes");

    if (pImgAttr != nullptr){
        cout << "successful read!!" << endl;
    } else{
        cout << "failed to read the xml file :(((" << endl;
    }

    const char* pFullResImagery = pImgAttr->FirstChildElement("fullResolutionImageData")->GetText();
    const char* pole = pImgAttr->FirstChildElement("fullResolutionImageData")->Attribute("pole");

    cout << "Full resolution imagery file: " << pFullResImagery << endl;
    cout << "Pole: " << pole << endl;

}

/// Reads information from R2 product.xml and lut*.xml files.
class R2XML
{
public:
    /// Loads the R2 XML Product Data: product info, info about each file, and the GCP list.
    static void LoadR2XML(const string sproductxmlpath, R2ProductInfo &r2p,
                          vector<R2FileInfo> &r2fileinfos, vector<GroundControlPoint> &r2gcp,
                          double &lonmin, double &latmin,
                          double &lonmax, double &latmax, int blockaverage)
    {
        try
        {
            r2fileinfos.clear();
            r2gcp.clear();

            // Load the information in the XML document.
            XMLDocument xmldoc;
            xmldoc.LoadFile(sproductxmlpath.c_str());

            XMLElement* pRoot = xmldoc.FirstChildElement("product");

            string passdirection;

            XMLElement* pImgAttr = pRoot->FirstChildElement("imageAttributes");

            for(XMLElement* e = pImgAttr->FirstChildElement("fullResolutionImageData"); e; e = e->NextSiblingElement("fullResolutionImageData"))
            {
                R2FileInfo  temp;
                temp.filename = e->GetText();
                temp.pole = e->Attribute("pole");
                r2fileinfos.push_back(temp);
            }

            XMLElement* pRasterAttr = pImgAttr->FirstChildElement("rasterAttributes");

            r2p.numberofsamples = pRasterAttr->FirstChildElement("numberOfSamplesPerLine")->IntText();
            r2p.numberoflines   = pRasterAttr->FirstChildElement("numberOfLines")->IntText();

            XMLElement*  pSourceAttr = pRoot->FirstChildElement("sourceAttributes");
            r2p.starttime  = pSourceAttr->FirstChildElement("rawDataStartTime")->GetText();

            r2p.productid = pRoot->FirstChildElement("productId")->GetText();
            r2p.producttype = pRoot->FirstChildElement("imageGenerationParameters")->FirstChildElement("generalProcessingInformation")->FirstChildElement("productType")->GetText();

            XMLElement*  pOrbAttr = pSourceAttr->FirstChildElement("orbitAndAttitude")->FirstChildElement("orbitInformation");
            r2p.passdirection = pOrbAttr->FirstChildElement("passDirection")->GetText();

            double line, pixel, lat, lon;

            lonmin = 1000;
            latmin = 1000;
            lonmax = -1000;
            latmax = -1000;

            XMLElement* pGeoGrid = pImgAttr->FirstChildElement("geographicInformation")->FirstChildElement("geolocationGrid");

            for(XMLElement* e = pGeoGrid->FirstChildElement("imageTiePoint"); e; e = e->NextSiblingElement("imageTiePoint"))
            {
                lat = e->FirstChildElement("geodeticCoordinate")->FirstChildElement("latitude")->DoubleText();
                lon = e->FirstChildElement("geodeticCoordinate")->FirstChildElement("longitude")->DoubleText();
                line  = e->FirstChildElement("imageCoordinate")->FirstChildElement("line")->DoubleText();
                pixel = e->FirstChildElement("imageCoordinate")->FirstChildElement("pixel")->DoubleText();

                GroundControlPoint tempgcp;
                // pixel and line coordinate with block average adjustment.
                tempgcp.pixel = (float) (pixel / (double)blockaverage);
                tempgcp.line = (float) (line / (double)blockaverage);
                tempgcp.lat = lat;
                tempgcp.lon = lon;
                r2gcp.push_back(tempgcp);

                // find min / max of lat and lon.
                if (lat < latmin)
                    latmin = lat;
                if (lon < lonmin)
                    lonmin = lon;
                if (lat > latmax)
                    latmax = lat;
                if (lon > lonmax)
                    lonmax = lon;
            }
        }
        catch(...)
        {
            throw string("XML Error in " + sproductxmlpath);
        }
    }

    /// Load the calibration gains from the look up table XML file.
    static void LoadR2LUT(const string s_lutxmlpath, double &offset, vector<double> &gains)
    {
        try
        {
            // Load the LUT file.
            XMLDocument xmldoc;
            xmldoc.LoadFile(s_lutxmlpath.c_str());
            XMLElement* pLut = xmldoc.FirstChildElement("lut");

            if (!pLut) throw std::runtime_error("this is not good, I was unable to find the root node!");

            // clear the gains.
            gains.clear();

            offset = pLut->FirstChildElement("offset")->DoubleText();

            // Get the gains
            string gainstring = pLut->FirstChildElement("gains")->GetText();

            // Create an istringstream out of it.
            istringstream iss(gainstring, istringstream::in);
            double gaintemp;

            // process it as an input stream and store it in the vector of gains.
            while(iss>>gaintemp)
            {
                gains.push_back(gaintemp);
            }
        }
        catch (...)
        {
            throw string("IO Error: Cannot read " + s_lutxmlpath);
        }
    }
//
//    // Ross April 2016
//    // Function to read sentinel-1 xml
//
//    static void LoadS1XML (String^ productPath, vector<GroundControlPoint>& imagetiepoints, int width)
//    {
//        try
//        {
//            XmlDocument^ xmldoc = gcnew XmlDocument();
//            xmldoc->Load(productPath);
//
//            // Resolution info.
//            XmlNodeList^ items = xmldoc->GetElementsByTagName("numberOfSamples");
//            int numberofsamples = Int32::Parse(items->Item(0)->InnerText);
//            if (numberofsamples <= 0)
//                numberofsamples = 1;
//            int blockaverage = (int)ceil((double)numberofsamples/(double)width);
//
//            // get the image file names
//            items = xmldoc->GetElementsByTagName("geolocationGridPoint");
//            XmlNodeList^ children;
//            GroundControlPoint tempgcp;
//            for (int i = 0; i < items->Count; i++)
//            {
//                children = items->Item(i)->ChildNodes;
//                for (int j = 0; j < children->Count; j++)
//                {
//                    if (children->Item(j)->Name == "line")
//                        tempgcp.line = (float)(Single::Parse(children->Item(j)->InnerText)/(double)blockaverage);
//                    else if (children->Item(j)->Name == "pixel")
//                        tempgcp.pixel = (float)(Single::Parse(children->Item(j)->InnerText)/(double)blockaverage);
//                    else if (children->Item(j)->Name == "latitude")
//                        tempgcp.lat = Double::Parse(children->Item(j)->InnerText);
//                    else if (children->Item(j)->Name == "longitude")
//                        tempgcp.lon = Double::Parse(children->Item(j)->InnerText);
//                }
//                imagetiepoints.push_back(tempgcp);
//            }
//        }
//        catch (System::IO::IOException^)
//        {
//            const char* name = (const char*)(Runtime::InteropServices::Marshal::StringToHGlobalAnsi(productPath)).ToPointer();
//            throw string("IO Error: Cannot read " + (string)name);
//        }
//        catch (XmlException^)
//        {
//            const char* name = (const char*)(Runtime::InteropServices::Marshal::StringToHGlobalAnsi(productPath)).ToPointer();
//            throw string("XML Error in " + (string)name);
//        }
//    }

};