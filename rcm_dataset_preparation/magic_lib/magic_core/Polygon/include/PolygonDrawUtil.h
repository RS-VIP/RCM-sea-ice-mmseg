/** \file
 *
 * \brief PolygonDrawUtil class
 *
 * \author Peter Yu
 *
 */

#pragma once

#include <vector>
#include "DataType.h"
#include "GrayImage.h"

/// defines clipping borders: top of image, left of image, bottom of image and right of image.
enum ClippingBorder { CLIP_TOP, CLIP_LEFT, CLIP_BOTTOM, CLIP_RIGHT };

/// Contains several routines related to drawing polygons.
class PolygonDrawUtil
{
public:

	/// Clips polygons to width * height image window starting at [0,0]
	static vector<POINTex> ClipPolygon(const vector<POINTex>& polygon, int width, int height);

	
	/// Fills polygon on image with value val. Polygon is clipped to the image width and height.
	template<typename T>
	static std::vector<POINTex> ClippedFillPolygon(const std::vector<POINTex> &polygon, TGrayImage<T>* img, T val, int &minl, int &minc, int &maxl, int &maxc, TGrayImage<BYTE>* mask)
	{
		vector<POINTex> clippedpolygon;

		if (!img)
			return clippedpolygon;

		clippedpolygon = ClipPolygon(polygon, img->Width(), img->Height());
		FillPolygon(clippedpolygon, img, val, minl, minc, maxl, maxc, mask);

		return clippedpolygon;
	}

	/// Simple polygon fill.
	template<typename T>
	static void FillPolygon(const vector<POINTex> &polypoints, TGrayImage<T>* img, T val, int &minl, int &minc, int &maxl, int &maxc, TGrayImage<BYTE>* mask)
	{
		if (!img)
			return;

		int nImgMaxCol = img->Width() - 1;
		int nImgMaxRow = img->Height() - 1;

		T** mat = img->GetData();

		///////////////////////////////////////////////////////////
		//
		// To understand this polygon filling algorithm, please
		// refer to this web site
		// http://alienryderflex.com/polygon_fill/
		// 
		// Some parts are modified to compatible to our requirement
		// -- Sharif
		//
		//////////////////////////////////////////////////////////
		
		
		long MAX_POLY_CORNERS = (long) polypoints.size();

		int  nodes, *nodeX, /*pixelX,*/ swap;
		long i, j;
		long pixelY;

		// declare new pointer so that I dn't have to change the existing code
		// nodeX will only store the x coordinate of any node because Y is fixed 
		int* polyX = NULL;
		int* polyY = NULL;
		polyX = new int[MAX_POLY_CORNERS];
		polyY = new int[MAX_POLY_CORNERS];

		for (long ij =0; ij < MAX_POLY_CORNERS; ij++)
		{		
			polyX[ij] = (int) polypoints[ij].Col;
			polyY[ij] = (int) polypoints[ij].Row;

			// do some bound checking to make sure this doesn't crash!
			if (polyX[ij] < 0)
				polyX[ij] = 0;
			if (polyY[ij] < 0)
				polyY[ij] = 0;
			if (polyX[ij] > nImgMaxCol)
				polyX[ij] = nImgMaxCol;
			if (polyY[ij] > nImgMaxRow)
				polyY[ij] = nImgMaxRow;			

		}

		nodeX = new int[MAX_POLY_CORNERS];


		// assigning incoming pointer to new pointer variable
		//polyX = col;
		//polyY = lin;


		// inilialize to some extreme unexpected value
		minl = 1000000;
		minc = 1000000;
		maxl = -1000000;
		maxc = -1000000;

		// get bounding rectangle coordinate points
		for ( i = 0; i < MAX_POLY_CORNERS; i++)
		{
			if (polyY[i] < minl)
				minl = polyY[i];
			if (polyY[i] > maxl)
				maxl = polyY[i];
			if (polyX[i] < minc)
				minc = polyX[i];
			if (polyX[i] > maxc)
				maxc = polyX[i];
		}

		// get and initialize top and botom of bounding polygon
		int IMAGE_TOP = minl;
		int IMAGE_BOT = maxl;
		int IMAGE_RIGHT = maxc;
		int	IMAGE_LEFT  = minc;
		double  y; // it is miising in the original code

		int yi;// = polyY[i];
		int yj;// = polyY[j];
		int xi;
		int xj;
		double dx;
		double dy;
		double interpolatedX = 0;
		// another issue here, should I include the IMAGE_BOT row
		//for (pixelY=IMAGE_TOP; pixelY<=IMAGE_BOT; pixelY++) {
		
		//  Loop through the rows of the image.
		//for (pixelY = IMAGE_TOP; pixelY < IMAGE_BOT; pixelY++) -- original code
		for (pixelY = IMAGE_TOP; pixelY < IMAGE_BOT; pixelY++)
		{
			  y = pixelY;

			  //  Build a list of nodes.
			  nodes=0; j = MAX_POLY_CORNERS - 1;//polyCorners-1;
			  for (i = 0; i < MAX_POLY_CORNERS; i++)//polyCorners; i++)
			  {
				yi = polyY[i];
				yj = polyY[j];
				xi = polyX[i];
				xj = polyX[j];
				if (yi < pixelY && yj >= pixelY)
				{
					dy = (polyY[j]-polyY[i]);
					dx = (polyX[j]-polyX[i]);
					if (dx < 0)
						interpolatedX = (double)polyX[i] - (y - polyY[i])* (abs(dx)/dy);
					else
						interpolatedX = (double)polyX[i] + (y - polyY[i])* (abs(dx)/dy);

					nodeX[nodes++] = (int)interpolatedX;

				}
				else if (yj < pixelY && yi >= pixelY)
				{
					dy = (polyY[i]-polyY[j]);
					dx = (polyX[i]-polyX[j]);
					
					if (dx < 0)
						interpolatedX = (double)polyX[j] - (y - polyY[j])*(abs(dx)/dy);
					else
						interpolatedX = (double)polyX[i]  + (y - polyY[i])*(abs(dx)/dy);

					nodeX[nodes++] = (int)interpolatedX;


				}
				j=i; 
			  }

			  //  Sort the nodes, via a simple �Bubble?sort.
			  i=0;
			  while (i < nodes-1)
			  {
				  
				if (nodeX[i] > nodeX[i+1])
				{
				  swap=nodeX[i]; nodeX[i]=nodeX[i+1]; nodeX[i+1]=swap; if (i) i--;
				}
				else
				{
				  i++;
				}
			  }

			  //  Fill the pixels between node pairs.
			  for (i = 0; i < nodes; i +=2)
			  {
					//if   (nodeX[i  ] >= IMAGE_RIGHT) break; --  original code
					if  (nodeX[i+1] >= IMAGE_LEFT )// (nodeX[i+1] > IMAGE_LEFT )
					{
						  //if (nodeX[i  ]< IMAGE_LEFT )
						//	  nodeX[i  ]=IMAGE_LEFT ;
						  //if (nodeX[i+1]> IMAGE_RIGHT)
						//	  nodeX[i+1]=IMAGE_RIGHT;
						  for (j= nodeX[i]; j < nodeX[i+1]; j++)
						  {
							  if(mask == NULL)
							  {
							      //fillPixel(j,pixelY);
							      mat[pixelY][j] = val;
							  }
							  else if((*mask)((int)pixelY,(int)j) == 0)
							  {
								  mat[pixelY][j] = 0;
							  }
							  else
							  {
								  mat[pixelY][j] = val;
							  }
						  }
					}
			  }
		} 

		for ( i = 0; i < MAX_POLY_CORNERS; i++)
		{
			mat[polyY[i]][polyX[i]]=val;
		}

		delete [] nodeX;
		delete [] polyX;// = new int[MAX_POLY_CORNERS];
		delete [] polyY;// = new int[MAX_POLY_CORNERS];
	}

protected:

	/// Find an intersection between a border (horizontal or vertical line) at borderposition and line defined by (prev, cur).
	static POINTex Intersect(const POINTex& cur, const POINTex& prev, ClippingBorder border, int borderposition);

	/// determines if a point is inside the clipping border at borderposition.
	static bool IsInside(const POINTex& point, ClippingBorder border, int borderposition);

	/// Performs Sutherland-Hodgman Clipping on polypoints against the specified border at borderposition.
	/**
	 * border is CLIP_TOP, CLIP_LEFT, CLIP_BOTTOM or CLIP_RIGHT, borderposition is where the border is 
	 * (e.g. CLIP_TOP with borderposition = 5 means that the polygon will be clipped at row 5 so that 
	 * it only appears below row 5)
	 * polypoints2 will store the clipped polygon.
	 */
	static void SutherlandHodgmanClipping(const vector<POINTex>& polypoints, vector<POINTex>& polypoints2, ClippingBorder border, int borderposition);

};
