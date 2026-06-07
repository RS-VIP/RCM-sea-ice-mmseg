#include "PolygonDrawUtil.h"

/// Clips polygons to width * height image window starting at [0,0]
vector<POINTex> PolygonDrawUtil::ClipPolygon(const vector<POINTex>& polygon, int width, int height)
{
	// need two buffers to do clipping.
	vector<POINTex> polypoints;
	vector<POINTex> polypoints2;

	// clip each side
	SutherlandHodgmanClipping(polygon, polypoints, CLIP_TOP, 0);
	SutherlandHodgmanClipping(polypoints, polypoints2, CLIP_LEFT, 0);
	SutherlandHodgmanClipping(polypoints2, polypoints, CLIP_BOTTOM, height - 1);
	SutherlandHodgmanClipping(polypoints, polypoints2, CLIP_RIGHT, width - 1);

	return polypoints2;
}

/// Find an intersection between a border (horizontal or vertical line) at borderposition and line defined by (prev, cur).
POINTex PolygonDrawUtil::Intersect(const POINTex& cur, const POINTex& prev, ClippingBorder border, int borderposition)
{
	POINTex ret;
	
	if (border == CLIP_TOP || border == CLIP_BOTTOM) // top border or bottom border
	{
		double slope = ((double)(cur.Col - prev.Col)) / ((double) (cur.Row - prev.Row));
		ret.Row = borderposition;
		ret.Col = (int) ((ret.Row - cur.Row) * slope + cur.Col);		
		return ret;
	}
	else // left or right border.
	{
		double slope = ((double)(cur.Row - prev.Row)) / ((double)(cur.Col - prev.Col));
		ret.Col = borderposition;
		ret.Row = (int) ((ret.Col - cur.Col) * slope + cur.Row);		
		return ret;
	}
}

/// determines if a point is inside the clipping border at borderposition.
bool PolygonDrawUtil::IsInside(const POINTex& point, ClippingBorder border, int borderposition)
{
	if (border == CLIP_TOP) // top border
	{
		return point.Row >= borderposition;
	}
	if (border == CLIP_LEFT) // left border
	{
		return point.Col >= borderposition;
	}
	if (border == CLIP_BOTTOM) // bottom border
	{
		return point.Row < borderposition;
	}
	if (border == CLIP_RIGHT) // right border
	{
		return point.Col < borderposition;
	}
	else
		return false;
}

/// Performs Sutherland-Hodgman Clipping on polypoints against the specified border at borderposition.
/**
 * border is CLIP_TOP, CLIP_LEFT, CLIP_BOTTOM or CLIP_RIGHT, borderposition is where the border is 
 * (e.g. CLIP_TOP with borderposition = 5 means that the polygon will be clipped at row 5 so that 
 * it only appears below row 5)
 * polypoints2 will store the clipped polygon.
 */
void PolygonDrawUtil::SutherlandHodgmanClipping(const vector<POINTex>& polypoints, vector<POINTex>& polypoints2, ClippingBorder border, int borderposition)
{
	polypoints2.clear();

	if (polypoints.size() < 3)
		return;

	polypoints2.reserve(polypoints.size());
	
	size_t length = polypoints.size();

	size_t prev;

	for (size_t i=0; i < length; i++)
	{
		if (i > 0)
			prev = i - 1;
		else
			prev = length - 1;
		

		if (IsInside(polypoints[i], border, borderposition))
		{
			if (!IsInside(polypoints[prev], border, borderposition))
			{
				polypoints2.push_back(Intersect(polypoints[i], polypoints[prev], border, borderposition));
			}
			polypoints2.push_back(polypoints[i]);
		}
		else if (IsInside(polypoints[prev], border, borderposition))
		{
			polypoints2.push_back(Intersect(polypoints[i], polypoints[prev], border, borderposition));
		}
	}

}