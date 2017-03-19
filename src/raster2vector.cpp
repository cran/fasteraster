/*
 raster2vector function to convert raster image matrix into list of polygons
 Copyright (C) 2016-2017  Andy Bosyi
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

#define FILL_NA DBL_MAX
#define POLY_NA -1
#define FENCE_NONE 0
#define FENCE_X 1
#define FENCE_Y 2
#define FENCE_XY 3
#define COORDS(x, y) ((x) * leny + (y))
#define FENCE_UNSET_DX (fence[COORDS(x + 1, y)] &= FENCE_Y)
#define FENCE_UNSET_DY (fence[COORDS(x, y + 1)] &= FENCE_X)
#define FENCE_OPEN_RIGHT (!(fence[COORDS(x + 1, y)] & FENCE_X))
#define FENCE_OPEN_DOWN (!(fence[COORDS(x, y)] & FENCE_Y))
#define FENCE_OPEN_LEFT (!(fence[COORDS(x, y)] & FENCE_X))
#define FENCE_OPEN_UP (!(fence[COORDS(x, y + 1)] & FENCE_Y))

struct polygon_t
{
  int x;
  int y;
  double fill;
};

struct vector_t
{
  int x;
  int y;
};

//' @useDynLib fasteraster
//' @importFrom Rcpp sourceCpp
//NULL

//' Converts raster image matrix into list of polygons.
//' 
//' Takes a raster matrix and vectorize it by polygons. Typically the nature of artefact images is linear
//' and can be vectorized in much more efficient way than draw a series of 90 degrees lines.
//'
//' @param raster input matrix of numeric values.
//' @param from lower (greater than) margin for recognizable polygon values.
//' @param to upper (less than or equal) margin.
//' @param step values gradient.
//' @param precision linearization precision in matrix cells.
//' @param exclave polygons may have other inside.
//' @return list of integer matrixes that represent polygon coordinates in two columns.
//' @examples library(datasets)
//' raster2vector(volcano, 100, 200, 60, 1, FALSE)
//' raster2vector(volcano, 120, 200, 20)
//' @export
// [[Rcpp::export]]
List raster2vector(NumericMatrix raster, double from = 0, double to = 1, double step = 0.1, int precision = 1, bool exclave = true)
{
  int lenx0 = raster.nrow();
  int leny0 = raster.ncol();
  
  std::vector<polygon_t> polyheads;
  int lenx = lenx0 + 1;
  int leny = leny0 + 1;
  int* fence = new int[lenx * leny];
  int* polyxy = new int[lenx * leny];
  int polylen = 0;
  
  for (int x = lenx0; x >= 0; x--)
    for (int y = leny0; y >= 0; y--)
    {
      // in all cases set fence to the lower left corner
      fence[COORDS(x, y)] = FENCE_XY;
      // and unknown polygonID
      polyxy[COORDS(x, y)] = POLY_NA;

      // margin fences
      if (x == lenx0 || y == leny0)
        continue;

      // get fills for current cell
      double fill = raster[x + y * lenx0];
            
      // background
      if (fill <= from || fill > to)
        continue;

      //get neighbour polygons
      int polyidx = polyxy[COORDS(x + 1, y)];
      int polyidy = polyxy[COORDS(x, y + 1)];
      
      // get fills for this and neighbour cells
      fill = ceilf(fill / step) * step;
      double filldx = polyidx == POLY_NA ? FILL_NA : polyheads[polyidx].fill;
      double filldy = polyidy == POLY_NA ? FILL_NA : polyheads[polyidy].fill;

      if (fill != filldx && fill != filldy)
      {
        // new fill, no neighbours like that
        // create new polygon head
        polygon_t pol;
        pol.x = x;
        pol.y = y;
        pol.fill = fill;
        polyheads.push_back(pol);
        polyxy[COORDS(x, y)] = polyheads.size() - 1;
        polylen ++;
      }
      else if (fill == filldx && fill != filldy)
      {
        // right neighbour
        FENCE_UNSET_DX;
        polyxy[COORDS(x, y)] = polyidx;
      }
      else if (fill != filldx && fill == filldy)
      {
        // upper neighbour
        FENCE_UNSET_DY;
        polyxy[COORDS(x, y)] = polyidy;
      }
      else
      {
        // both neighbours - right and upper
        
        if (polyidx != polyidy)
        {
          // welding different polygons with the same fill
          polyheads[polyidy].fill = FILL_NA; // mark as removed
          polylen--;
          
          for (int i = leny0; i > y; i--)
            if (polyxy[COORDS(x, i)] == polyidy)
              polyxy[COORDS(x, i)] = polyidx;
          FENCE_UNSET_DX;
          FENCE_UNSET_DY;
          polyxy[COORDS(x, y)] = polyidx;
        }
        else
        {
          // same fills, same polygons
          
          if ((fence[COORDS(x + 1, y + 1)] & FENCE_Y) && !exclave)
          {
            // there is an internal polygon inside
            
            // keep the horizontal fence and weld to right
            FENCE_UNSET_DX;
            polyxy[COORDS(x, y)] = polyidx;
          }
          else
          {
            // weld to both
            FENCE_UNSET_DX;
            FENCE_UNSET_DY;
            polyxy[COORDS(x, y)] = polyidx;
          }
        }
        
      }
    }
  
  List polygons(polylen);
  std::vector<vector_t> pol;
  std::vector<vector_t> pol1;
  std::vector<vector_t> pol2;
  std::vector<vector_t> pol3;
  
  int polyidx = 0;
  for (int p = polyheads.size() - 1; p >= 0; p--)
  {
    polygon_t ph = polyheads[p];
    if (ph.fill == FILL_NA)
      continue; // omit joined pieces
   
    pol.clear();

    int x = ph.x;
    int y = ph.y;
    vector_t v, nv;
    v.x = 0;
    v.y = 0;
    
    while (true)
    {
      nv.x = 0;
      nv.y = 0;
      
      if (v.x == 0 && v.y == -1)
      {
        // down
        if (FENCE_OPEN_RIGHT)
          nv.x = 1;
        else
          if (FENCE_OPEN_DOWN)
            nv.y = -1;
          else
            if (FENCE_OPEN_LEFT)
              nv.x = -1;
            else
              nv.y = 1;
      }
      else if (v.x == -1 && v.y == 0)
      {
        // left
        if (FENCE_OPEN_DOWN)
          nv.y = -1;
        else
          if (FENCE_OPEN_LEFT)
            nv.x = -1;
          else
            if (FENCE_OPEN_UP)
              nv.y = 1;
            else
              nv.x = 1;
      }
      else if (v.x == 0 && v.y == 1)
      {
        // up
        if (FENCE_OPEN_LEFT)
          nv.x = -1;
        else
          if (FENCE_OPEN_UP)
            nv.y = 1;
          else
            if (FENCE_OPEN_RIGHT)
              nv.x = 1;
            else
              nv.y = -1;
      }
      else if (v.x == 1 && v.y == 0)
      {
        // right
        if (FENCE_OPEN_UP)
          nv.y = 1;
        else
          if (FENCE_OPEN_RIGHT)
            nv.x = 1;
          else
            if (FENCE_OPEN_DOWN)
              nv.y = -1;
            else
              nv.x = -1;
      }
      else
      {
        // first cell of the polygon (v.x == 0 && v.y == 0)
        if (FENCE_OPEN_RIGHT)
          nv.x = 1;
        else
          if (FENCE_OPEN_DOWN)
            nv.y = -1;
          else
            if (FENCE_OPEN_LEFT)
              nv.x = -1;
            else
              if (FENCE_OPEN_UP)
                nv.y = 1;
            else
            {
              // one-point polygon
              pol.push_back(nv);
              break;
            }
      }

      if (pol.size() > 0 && x == ph.x && y == ph.y && pol[0].x == nv.x && pol[0].y == nv.y)
        break; // not just started, same point, same vector

      pol.push_back(nv);
      
      x += nv.x;
      y += nv.y;
      v = nv;
    };

    // join 90 deg lines - optimize polygon
    pol1.clear();
    v.x = INT_MAX;
    v.y = INT_MAX;
    for (int i = 0; i < (int) pol.size(); i++)
    {
      nv = pol[i];
      if (v.x == nv.x && v.y == nv.y)
      {
        pol1.back().x += nv.x;
        pol1.back().y += nv.y;
      }
      else
      {
        pol1.push_back(nv);
        v = nv;
      }
    }
    
    //find inclined lines - optimize polygon
    pol2.clear();
    v.x = INT_MAX;
    v.y = INT_MAX;
    for (int i = 0; i < (int) pol1.size(); i++)
    {
      nv = pol1[i];
      if ((v.x >= -precision && v.x <= precision && v.x != 0 && nv.x == 0)
            || (v.y >= -precision && v.y <= precision && v.y != 0 && nv.y == 0))
      {
        pol2.back().x += nv.x;
        pol2.back().y += nv.y;
      }
      else
      {
        pol2.push_back(nv);
        v = nv;
      }
    }
    
    //find similar vectors - optimize polygon
    pol3.clear();
    int oi = 0;
    pol3.push_back(pol2[oi]);
    
    for (int i = 1; i < (int) pol2.size(); i++)
    {
      nv = pol2[i];

      x = 0;
      y = 0;
      bool fit = true;
      for (int j = oi; j <= i; j++)
      {
        //check for distance 
        x += pol2[j].x;
        y += pol2[j].y;
        
        double a = atan2((double) y, (double) x) - atan2((double) pol3.back().y, (double) pol3.back().x);
        if (a > PI)
          a -= 2 * PI;
        if (a < PI)
          a += 2 * PI;
        if (fabs(sin(a) * sqrt((double) (x * x) + (double) (y * y))) > (double) precision)
        {
          fit = false;
          break;
        }
      }
      
      if (fit)
      {
        pol3.back().x += nv.x;
        pol3.back().y += nv.y;
      }
      else
      {
        pol3.push_back(nv);
        oi = i;
      }
    }
 
    // fill result data
    IntegerMatrix mat(pol3.size(), 2);
    x = ph.x + 1; // convert from 0 to 1 based array
    y = ph.y + 1;
    mat(0, 0) = x;
    mat(0, 1) = y;
    for (int i = 0 ; i < (int) pol3.size() - 1; i ++)
    {
      x += pol3[i].x;
      y += pol3[i].y;
      mat(i + 1, 0) = x;
      mat(i + 1, 1) = y;
    }
    
    polygons(polyidx++) = mat;
  }
  
  
/*test function
  IntegerMatrix mat(lenx + 1, leny + 1);
  for (int x = lenx; x >= 0; x--)
    for (int y = leny; y >= 0; y--)
      mat(x, y) = fence[x][y];
  
  polygons(0) = mat;
  */

  delete[] fence;
  delete[] polyxy; 
  
  return polygons;
} //-------------------------------------------------------------------------------


struct polygonEx_t
{
  int headX;
  int headY;
  double fill;
  double centerX;
  double centerY;
  double weight;
};

#define WELD_CELL(p) { polyxy[COORDS(x, y)] = p; \
  polyheads[p].weight ++; \
  polyheads[p].centerX += (((double) x) - polyheads[p].centerX) / polyheads[p].weight; \
  polyheads[p].centerY += (((double) y) - polyheads[p].centerY) / polyheads[p].weight; }
  
//' Recognizes zones in the raster matrix and return a list of their analysis.
//' 
//' Takes a raster matrix and recognize zones exactly as the raster2vector function does. However does not return
//' polygon vectors, but weight and center of mass for each polygone (or so called zone).
//'
//' @param raster input matrix of numeric values.
//' @param from lower (greater than) margin for recognizable zone values.
//' @param to upper (less than or equal) margin.
//' @param step values gradient.
//' @return matrix of zones by four values: fill (ceiling values), weight (in cells), x and y of the mass center.
//' @examples library(datasets)
//' rasterZoneAnalyzer(volcano, 120, 200, 20)
//' @export
// [[Rcpp::export]]
NumericMatrix rasterZoneAnalyzer(NumericMatrix raster, double from = 0, double to = 1, double step = 0.1)
{
  int lenx0 = raster.nrow();
  int leny0 = raster.ncol();
  
  std::vector<polygonEx_t> polyheads;
  int lenx = lenx0 + 1;
  int leny = leny0 + 1;
  int* polyxy = new int[lenx * leny];
  int polylen = 0;
  
  for (int x = lenx0; x >= 0; x--)
    for (int y = leny0; y >= 0; y--)
    {
      // in all cases set unknown polygonID
      polyxy[COORDS(x, y)] = POLY_NA;
      
      // margin fences
      if (x == lenx0 || y == leny0)
        continue;
      
      // get fills for current cell
      double fill = raster[x + y * lenx0];
      
      // background
      if (fill <= from || fill > to)
        continue;
      
      //get neighbour polygons
      int polyidx = polyxy[COORDS(x + 1, y)];
      int polyidy = polyxy[COORDS(x, y + 1)];
      
      // get fills for this and neighbour cells
      fill = ceilf(fill / step) * step;
      double filldx = polyidx == POLY_NA ? FILL_NA : polyheads[polyidx].fill;
      double filldy = polyidy == POLY_NA ? FILL_NA : polyheads[polyidy].fill;
      
      if (fill != filldx && fill != filldy)
      {
        // new fill, no neighbours like that
        // create a new polygon head
        polygonEx_t pol;
        pol.headX = x;
        pol.headY = y;
        pol.fill = fill;
        pol.centerX = (double) x;
        pol.centerY = (double) y;
        pol.weight = (double) 1;
        polyheads.push_back(pol);
        
        polyxy[COORDS(x, y)] = polyheads.size() - 1;
        polylen ++;
      }
      else if (fill == filldx && fill != filldy)
      {
        // right neighbour
        WELD_CELL(polyidx);
      }
      else if (fill != filldx && fill == filldy)
      {
        // upper neighbour
        WELD_CELL(polyidy);
      }
      else
      {
        // both neighbours - right and upper
        WELD_CELL(polyidx); //weld to right
        
        if (polyidx != polyidy)
        {
          // welding different polygons with the same fill
          polyheads[polyidx].weight += polyheads[polyidy].weight;
          polyheads[polyidx].centerX += (polyheads[polyidy].centerX - polyheads[polyidx].centerX) * polyheads[polyidy].weight / polyheads[polyidx].weight;
          polyheads[polyidx].centerY += (polyheads[polyidy].centerY - polyheads[polyidx].centerY) * polyheads[polyidy].weight / polyheads[polyidx].weight;

          polyheads[polyidy].fill = FILL_NA; // mark as removed
          polylen--;
          
          for (int i = leny0; i > y; i--) // mark on the matrix
            if (polyxy[COORDS(x, i)] == polyidy)
              polyxy[COORDS(x, i)] = polyidx;
        }
      }
    }
    
    
  NumericMatrix zones(polylen, 4);

  int zoneidx = 0;
  for (int p = polyheads.size() - 1; p >= 0; p--)
  {
    polygonEx_t ph = polyheads[p];
    
    if (ph.fill == FILL_NA)
      continue; // omit joined pieces

    // fill result data
    zones(zoneidx, 0) = ph.fill;
    zones(zoneidx, 1) = ph.weight;
    zones(zoneidx, 2) = ph.centerX + 1; // convert from 0 to 1 based array
    zones(zoneidx, 3) = ph.centerY + 1;
    zoneidx ++;
  }

  delete[] polyxy; 
  return zones;
}
