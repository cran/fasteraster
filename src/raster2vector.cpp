/*
 raster2vector function to convert raster image matrix into list of polygons
 Copyright (C) 2016  Andy Bosyi
 
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
#include <math.h>

using namespace Rcpp;

#define FILL_NA DBL_MAX
#define POLY_NA -1
#define FENCE_NONE 0
#define FENCE_X 1
#define FENCE_Y 2
#define FENCE_XY 3
#define FENCE_UNSET_DX fence[x + 1][y] &= FENCE_Y
#define FENCE_UNSET_DY fence[x][y + 1] &= FENCE_X
#define FENCE_OPEN_RIGHT !(fence[x + 1][y] & FENCE_X)
#define FENCE_OPEN_DOWN !(fence[x][y] & FENCE_Y)
#define FENCE_OPEN_LEFT !(fence[x][y] & FENCE_X)
#define FENCE_OPEN_UP !(fence[x][y + 1] & FENCE_Y)

struct polygon_t
{
  int x;
  int y;
  double f;
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
//' @examples raster2vector(volcano, 120, 200, 20)
//' @export
// [[Rcpp::export]]
List raster2vector(NumericMatrix raster, double from = 0, double to = 1, double step = 0.1, int precision = 1, bool exclave = true)
{
  int lenx = raster.nrow();
  int leny = raster.ncol();
  
  std::vector<polygon_t> polyheads;
  int fence[lenx + 1][leny + 1];
  int polyxy[lenx + 1][leny + 1];

  for (int x = lenx; x >= 0; x--)
    for (int y = leny; y >= 0; y--)
    {
      // anyway set fence to the lower left corner
      fence[x][y] = FENCE_XY;

      // get fills for current cell
      double f = raster(x, y);
            
      // margin fences or background
      if (x == lenx || y == leny || f <= from || f > to)
      {
        polyxy[x][y] = POLY_NA;
        continue;
      }

      //get neighbour polygons
      int pdx = polyxy[x + 1][y];
      int pdy = polyxy[x][y + 1];
      
      // get fills for this and neighbour cells
      f = ceilf(f / step) * step;
      double fdx = pdx == POLY_NA ? FILL_NA : polyheads[pdx].f;
      double fdy = pdy == POLY_NA ? FILL_NA : polyheads[pdy].f;

      if (f != fdx && f != fdy)
      {
        // new fill, no neighbours like that
        // create new polygon head
        polygon_t p;
        p.x = x;
        p.y = y;
        p.f = f;
        polyheads.push_back(p);
        polyxy[x][y] = polyheads.size() - 1;
      }
      else if (f == fdx && f != fdy)
      {
        // right neighbour
        FENCE_UNSET_DX;
        polyxy[x][y] = polyxy[x + 1][y];        
      }
      else if (f != fdx && f == fdy)
      {
        // upper neighbour
        FENCE_UNSET_DY;
        polyxy[x][y] = polyxy[x][y + 1];
      }
      else
      {
        // both neighbours - right and upper
        
        if (pdx != pdy)
        {
          // welding different polygons with the same fill
          polyheads[pdy].f = FILL_NA; // mark as removed
          for (int i = leny - 1; i > y; i--)
            if (polyxy[x][i] == pdy)
              polyxy[x][i] = pdx;
          FENCE_UNSET_DX;
          FENCE_UNSET_DY;
          polyxy[x][y] = pdx;
        }
        else
        {
          // same fills, same polygons
          
          if (fence[x + 1][y + 1] & FENCE_Y && !exclave)
          {
            // there is an internal polygon inside
            
            // keep the horizontal fence and weld to right
            FENCE_UNSET_DX;
            polyxy[x][y] = pdx;
          }
          else
          {
            // weld to both
            FENCE_UNSET_DX;
            FENCE_UNSET_DY;
            polyxy[x][y] = pdx;
          }
        }
        
      }
    }
  
  int polylen = 0;
  for (int p = polyheads.size() - 1; p >= 0; p--)
    if (polyheads[p].f != FILL_NA)
      polylen++;
  
  List polygons(polylen);
  std::vector<vector_t> pol;
  std::vector<vector_t> pol1;
  std::vector<vector_t> pol2;
  std::vector<vector_t> pol3;
  
  int polyidx = 0;
  for (int p = polyheads.size() - 1; p >= 0; p--)
  {
    if (polyheads[p].f == FILL_NA)
      continue; // omit joined pieces
   
    pol.clear();
      
    polygon_t ph = polyheads[p];
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

      if (pol.size() > 0 && x == ph.x && y == ph.y & pol[0].x == nv.x && pol[0].y == nv.y)
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
    for (int i = 0; i < pol.size(); i++)
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
    for (int i = 0; i < pol1.size(); i++)
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
    
    for (int i = 1; i < pol2.size(); i++)
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
        
        double a = atan2(y, x) - atan2(pol3.back().y, pol3.back().x);
        if (a > PI)
          a -= 2 * PI;
        if (a < PI)
          a += 2 * PI;
        if (fabs(sin(a) * sqrt(x * x + y * y)) > precision)
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
    x = ph.x;
    y = ph.y;
    mat(0, 0) = x;
    mat(0, 1) = y;
    for (int i = 0 ; i < pol3.size() - 1; i ++)
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

  return polygons;
}
