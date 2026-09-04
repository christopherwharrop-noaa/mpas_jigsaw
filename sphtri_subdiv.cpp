#include <cstdint>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <limits>
#include <fstream>
#include <algorithm>

#ifndef MAX_EDGES_ON_PT
#define MAX_EDGES_ON_PT 10
#endif

#ifdef DEBUG
#define DEBUG_OUT( s ) std::cout << s
#else
#define DEBUG_OUT( s )
#endif

////////////////////////////////////////////////////////////////////////////////
/// string helper funcs
std::string
trim( std::string input )
{
  std::string::iterator first = std::find_if_not(
                                                input.begin(),
                                                input.end(),
                                                []( char c ) { return std::isspace( static_cast< int >( c ) ); }
                                                );


  if ( first == input.end() )
  {
    // entirely space, why are you here?
    return std::string();
  }

  size_t pos = first - input.begin();
  size_t len = std::find_if_not(
                                input.rbegin(),
                                input.rend(),
                                []( char c ) { return std::isspace( static_cast< int >( c ) ); }
                                ) - input.rbegin();

  return input.substr( pos, input.size() - len - pos );
}

std::vector<std::string>
tokenize( std::string input, std::string delim, bool trimTokens=false )
{
  size_t pos = 0;
  size_t prev = 0;
  std::vector< std::string > tokens;
  std::string token;

  while ( ( pos = input.find( delim, prev ) ) != std::string::npos )
  {
    token = input.substr(prev, pos - prev);
    if ( trimTokens )
    {
      token = trim( token );
    }
    tokens.push_back( token );
    prev = pos + delim.length();
  }
  // get last token
  token = input.substr(prev);
  if ( trimTokens )
  {
    token = trim( token );
  }
  tokens.push_back( token );
  return tokens;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Simple class to read/write and contain JIGSAW mesh info
////////////////////////////////////////////////////////////////////////////////
class JMesh
{
  public:
    JMesh( )
    : dims_( 0 )
    , meshID_( 9999 )
    , meshType_( "" )
    {}

    void read( std::string filename )
    {
      std::cout << "Reading mesh file " << filename << std::endl;
      std::ifstream file( filename );
      std::string line;

      if ( !file.is_open() )
      {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        throw 1;
      }
      
      size_t pos;
      uint32_t npoints        = 0;
      uint32_t npointsRead    = 0;
      uint32_t ntriangles     = 0;
      uint32_t ntrianglesRead = 0;
      bool pointRead = false;
      bool triRead = false;

      while ( std::getline( file, line ) )
      {
        if ( pointRead )
        {
          if ( npointsRead == npoints ) { pointRead = false; }
          else
          {
            // Line contains x;y;z;? info
            std::vector< std::string > tokens = tokenize( line, ";" );
            for ( uint8_t i = 0; i < dims_; i++ )
            {
              points_[npointsRead * dims_ + i] = std::stod( tokens[i] );
            }
            npointsRead++;
            continue;
          }
        }
        else if ( triRead )
        {
          if ( ntrianglesRead == ntriangles ) { triRead = false; }
          else
          {
            // Line contains p1;p2;p3;? info
            std::vector< std::string > tokens = tokenize( line, ";" );
            for ( uint8_t i = 0; i < dims_; i++ )
            {
              triangles_[ntrianglesRead * dims_ + i] = std::stod( tokens[i] );
            }
            ntrianglesRead++;
            continue;
          }
        }

        if ( npoints == 0 )
        {
          pos = line.find( "POINT=" );
          if ( pos != std::string::npos )
          {
            npoints = std::stoi( line.substr( pos + 6 ) );
            points_.resize( npoints * 3 );
            pointRead = true;
            std::cout << "  Points   : " << npoints << std::endl;
            continue;
          }
        }

        if ( ntriangles == 0 )
        {
          pos = line.find( "TRIA3=" );
          if ( pos != std::string::npos )
          {
            ntriangles = std::stoi( line.substr( pos + 6 ) );
            triangles_.resize( ntriangles * 3 );
            triRead = true;
            std::cout << "  Triangles: " << ntriangles << std::endl;
            continue;
          }
        }

        // Try finding some mesh info
        if ( dims_ == 0 )
        {
          pos = line.find( "NDIMS=" );
          if ( pos != std::string::npos )
          {
            dims_ = std::stoi( line.substr( pos + 6 ) );
            std::cout << "  Mesh dims: " << dims_ << std::endl;
          }
        }

        if ( meshID_ == 9999 )
        {
          pos = line.find( "MSHID=" );
          if ( pos != std::string::npos )
          {
            std::vector< std::string > tokens = tokenize( line.substr( pos + 6 ), ";" );
            meshID_ = std::stoi( tokens[0] );
            meshType_ = tokens[1];
            std::cout << "  Mesh ID  : " << meshID_ << std::endl;
            std::cout << "  Typename : " << meshType_ << std::endl;
          }
        }

        // Unknown line
      }

      std::cout << "  Finished reading mesh file" << std::endl;

      file.close();
    }

    void write( std::string filename )
    {
      std::cout << "Writing mesh file " << filename << std::endl;

      std::ofstream file( filename );
      file << std::setprecision(std::numeric_limits<double>::max_digits10);

      file << "# " << filename << "; created by SUBDIVISION TOOL\n"
           << "MSHID=" << meshID_ << ";" << meshType_ << "\n"
           << "NDIMS=" << dims_ << "\n";

      // Points
      file << "POINT=" << point_count() << "\n";
      for ( uint32_t i = 0; i < point_count(); i++ )
      {
        for ( uint8_t j = 0; j < dims_; j++ )
        {
          file << points_[i * dims_ + j] << ";";
        }
        // JIGSAW has a number at the end, don't know what it means...
        file << "0\n";
      }

      // Triangles
      file << "TRIA" << dims_ << "=" << triangle_count() << "\n";
      for ( uint32_t i = 0; i < triangle_count(); i++ )
      {
        for ( uint8_t j = 0; j < dims_; j++ )
        {
          file << triangles_[i * dims_ + j] << ";";
        }
        // JIGSAW has a number at the end, don't know what it means...
        file << "0\n";
      }

      file.close();
    }

    uint32_t triangle_count( ) { return triangles_.size() / dims_; }
    uint32_t point_count( ) { return points_.size() / dims_; }

  // We need to manipulate these fields a bit
  // private:
    uint32_t                dims_;
    uint32_t                meshID_;
    std::string             meshType_;
    std::vector< double   > points_;
    std::vector< uint32_t > triangles_;
};


namespace math
{

double
dot(
    double *a, ///< 3D point
    double *b  ///< 3D point
    )
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double
norm(
      double *a ///< 3D point
    )
{
  return std::sqrt( dot( a, a ) );
}

void
normalize(
            double *a ///< 3D point
          )
{
  double mag = norm( a );
  for ( uint8_t i = 0; i < 3; i++ ) { a[i] /= mag; }
}

////////////////////////////////////////////////////////////////////////////////
/// Small-angle tolerant calculation of angle between two 3D points in radians
///
/// https://people.eecs.berkeley.edu/~wkahan/Mindless.pdf Section 12
/// This is a slightly slower but more accurate and tolerant to singularities
/// take on the classical acos( dot( a, b ) )
/// Some optimization is made assuming a and b are unit vectors
////////////////////////////////////////////////////////////////////////////////
double
geodesic_dist(
              double *a, ///< 3D unit vector
              double *b  ///< 3D unit vector
              )
{
  double tmpA[ 3 ] = { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
  double tmpB[ 3 ] = { a[0] + b[0], a[1] + b[1], a[2] + b[2] };

  return 2.0 * std::atan2( norm( tmpA ), norm( tmpB ) );
}

////////////////////////////////////////////////////////////////////////////////
/// Spherical linear interpolation between two points a->b where an input value
/// of 0 denotes fully at point a, 1 denotes point b, and any value between is
/// some linear proportion along the arc length of the spherical geodesic.
////////////////////////////////////////////////////////////////////////////////
void
slerp( 
      double   *a,   ///< 3D unit vector
      double   *b,   ///< 3D unit vector
      double   *out, ///< 3D unit vector
      double    at   ///< Linear parameter of where we should interp to
      )
{
  double omega = geodesic_dist( a, b );
  double coefA = std::sin( omega * ( 1.0 - at ) ) / std::sin( omega );
  double coefB = std::sin( omega * ( at ) ) / std::sin( omega );

  out[0] = a[0] * coefA + b[0] * coefB;
  out[1] = a[1] * coefA + b[1] * coefB;
  out[2] = a[2] * coefA + b[2] * coefB;

  normalize( out );
}

////////////////////////////////////////////////////////////////////////////////
/// Partial sum of the first n natural numbers
////////////////////////////////////////////////////////////////////////////////
uint32_t
partial_sum( uint32_t n )
{
  return n * ( n + 1 ) / 2;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert the barycentric weights to a corresponding point inside the spherical
/// triangle. A normal planar barycentric calculation yields distortions close to
/// the vertices of the original triangle for large triangles, so slerp is used
/// instead.
///
/// This functions by taking the weights lambda<a,b,c> and constructing a geodesic
/// path from point A to the opposing edge BC that naturally intersects the point
/// to be found.
///
/// The point on BC is determined by taking lambda_a and distributing
/// the weight proportionally to lambda_b and lambda_c based on the ratio of
/// their sum - that is if lambda_c = 0, regardless of lambda_b, the entire weight
/// of lambda_a is added to lambda_b. The value lambda_b is then used as the
/// linear weight along BC using slerp to generate the point opposite point A.
/// This point P now functions as the start point for a slerp calculation from
/// point P -> point A. The linear weight is equal to lambda_a.
///
/// As barycentric weights inform proportionally how "close" a point is to any
/// given point, this geodesic linear approximation using slerp translates well.
////////////////////////////////////////////////////////////////////////////////
void
barycentric_to_spherical(
                          double *tri,     ///< 3D tri of [a,b,c] where each point is <xyz>
                          float  *weights, ///< 3 barycentric weights for a,b,c
                          double *out      ///< 3D point xyz
                          )
{
  // About 2.5x slower but more accurate
  // Trivial case of [1, 0, 0] weights to avoid divide by zero
  if ( weights[0] == 1.0 ) { std::memcpy( out, &tri[0], sizeof(tri[0]) * 3 ); }

  float bc_weight = weights[1] + weights[0] * ( weights[1] / ( weights[1] + weights[2] ) );
  double bc_pt[3];
  slerp( &tri[6], &tri[3], bc_pt, bc_weight );
  slerp( bc_pt, &tri[0], out, weights[0] );

  // // This is just a sanity planar projected
  // out[0] = 0.0;
  // out[1] = 0.0;
  // out[2] = 0.0;
  // for ( uint8_t k = 0; k < 3; k++ )
  // {
  //   out[0] += weights[k] * tri[k * 3 + 0];
  //   out[1] += weights[k] * tri[k * 3 + 1];
  //   out[2] += weights[k] * tri[k * 3 + 2];
  // }
  // normalize( out );
}

} // namespace math


namespace tri_subdiv
{

struct PointInfo
{
  uint32_t ptID;
  uint32_t ptLayer;
  uint32_t ptIdx;
  int8_t   ptType;
  uint32_t infoIdx; ///< For original points the corresponding original index
                    ///< For exterior points the corresponding original edge index
                    ///< For interior points the corresponding index when only counting these pts
};

////////////////////////////////////////////////////////////////////////////////
/// Subdivide a triangle by an integer number of layers such that each
/// subsequent layer adds one more subdivision point to each side of the main
/// triangle.
///
/// For example, a triangle using 3 layers results in the following subdivision:
///              * C
///       S #2  / \          // Note that along any side there are 3 "segments"
///      R     /8,4\         // of triangles. There will always be layers + 1
///     E     *-----*        // points along each side of the main triangle
///    Y #1  / \7,3/ \       //
///   A     /3,2\ /6,2\      // Subdivided triangles are annotaed with N,M
///  L     *-----*-----*     // where N = index into total subdivided triangles
///   #0  / \2,1/ \5,1/ \    // and M = index into local subdivided layer
///      /0,0\ /1,0\ /4,0\   //
///     *-----*-----*-----*  // The input triangle is defined by ABC (ccw)
///     A  ^               B
///         \  index into layers
///
/// Every return triangle will be A'B'C' (ccw) whose geometry is, nominally,
///        * C'     C'*-----* B'
///       / \    or    \ Y /
///      / X \          \ /
///     *-----*          * A'
///     A'    B'
///
/// Where the main triangle is composed of either X or Y orientation style
/// triangles.
///
/// The fully subdivided triangle will be indexed by a unique relative set of
/// points as shown below:
///           #3
///             9* C
///       S     / \          // The point regime will have one extra layer than
///      R #2  /   \         // the subdivision
///     E    5*----8*        //
///    Y     / \   / \       //
///   A #1  /   \ /   \      //
///  L    2*----4*----7*     //
///       / \   / \   / \    //
///  #0  /   \ /   \ /   \   //
///    0*----1*----3*----6*  //
///     A  ^               B
///         \  index into layers
///
/// The return trianle will be the point ids in this subdivision
///
////////////////////////////////////////////////////////////////////////////////
void
subdiv_tri_ptids(
                  uint32_t  triID,    ///< relative subdivided triangle id
                  uint32_t  *pointIDs ///< relative subdivided point ids for tri
                  )
{
  uint32_t i = static_cast< uint32_t >( std::sqrt( triID ) );
  uint32_t j = triID - i * i;

  // Find pt id of A' by finding decomposed ptLayer and ptIdx from i/j
  // A' will always be on the interior layer
  uint32_t ptLayerA = i;
  uint32_t ptIdxA = static_cast< uint32_t >( j / 2 );
  pointIDs[0] = math::partial_sum( ptLayerA ) + ptIdxA;

  if ( j % 2 == 0 )
  {
    // B' is same ptIdx but one layer over and C' is one index greater in B' layer
    pointIDs[1] = math::partial_sum( ptLayerA + 1 ) + ptIdxA;
    pointIDs[2] = pointIDs[1] + 1;
  }
  else
  {
    // C' is one index greater than A' and B' is one greater layer and index 
    pointIDs[2] = pointIDs[0] + 1;
    pointIDs[1] = math::partial_sum( ptLayerA + 1 ) + ptIdxA + 1;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the barycentric weights of the subdivided point ID based on the
/// subdivision layout noted in subdiv_tri_ptids()
////////////////////////////////////////////////////////////////////////////////
void
subdiv_ptid_to_barycentric(
                            uint32_t pointID,
                            uint32_t numTriLayers,
                            float    *weights
                            )
{
  uint32_t ptLayer = std::floor( ( std::sqrt( 8.0f * pointID + 1.0f ) - 1.0f ) / 2.0f );
  uint32_t ptIdx   = pointID - math::partial_sum( ptLayer );

  // Weight for A in base triangle
  weights[0] = static_cast< float >( numTriLayers - ptLayer ) / numTriLayers;
  // Weight for B in base triangle
  weights[1] = static_cast< float >( ptLayer - ptIdx ) / numTriLayers;
  // Weight for C in base triangle
  weights[2] = static_cast< float >( ptIdx ) / numTriLayers;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the info of point this ID is in a subdivision.
///
/// https://en.wikipedia.org/wiki/Triangular_number
/// https://math.stackexchange.com/a/1417583
///
/// The layer is the greatest value x yielding a triangular number less than the ID
/// The index in the layer is the ID minus same triangular number of x
///
/// The ptType will be:
///   -1 = Interior
///    0 = Original
///    1 = Exterior
///
/// The infoIdx will be filled according to the type:
///  * For original points the value will be original point index (0, 1, or 2)
///
///  * For exterior points the value will be the index of the edge in the original
///    triangle that precedes the edge this point lies on (0, 1, or 2)
///    For triangle ABC, edge AB = 0, BC = 1, CA = 2
///           C
///        2 / \ 1
///         A---B
///           0
///
///  * For interior points the value will be the point id as it corresponds to
///    the unique incremental ordering of only interior points. This follows the
///    same layering and index incrementation as the fully subdivided triangle
///    but only counts interior points.
///
////////////////////////////////////////////////////////////////////////////////
PointInfo
subdiv_ptid_info(
                  uint32_t pointID,     ///< relative point id
                  uint32_t numTriLayers
                  )
{
  PointInfo info;
  info.ptLayer = std::floor( ( std::sqrt( 8.0f * pointID + 1.0f ) - 1.0f ) / 2.0f );
  info.ptIdx   = pointID - math::partial_sum( info.ptLayer );

  uint32_t maxID = math::partial_sum( numTriLayers + 1 ) - 1;

  info.ptID = pointID;
  // Original points
  if ( pointID == 0 || pointID == maxID || pointID == (maxID - numTriLayers) )
  {
    info.ptType = 0;
    if      ( pointID == 0 )     { info.infoIdx = 0; }
    else if ( pointID == maxID ) { info.infoIdx = 2; }
    else                         { info.infoIdx = 1; }
  }
  // Obvious exterior points (lies on: AB or BC or CA)
  else if ( info.ptIdx == 0 || ( pointID > (maxID - numTriLayers) ) || info.ptIdx == info.ptLayer )
  {
    info.ptType = 1;
    if      ( info.ptIdx == 0 )            { info.infoIdx = 0; }
    else if ( info.ptIdx == info.ptLayer ) { info.infoIdx = 2; }
    else                         { info.infoIdx = 1; }
  }
  // Must be an interior point (not original and not exterior)
  else
  {
    info.ptType = -1;
    info.infoIdx = math::partial_sum( info.ptLayer - 2 ) + info.ptIdx - 1;
  }
  return info;
}


////////////////////////////////////////////////////////////////////////////////
/// Generate the relative weights and point IDs of a triangle subdivided into
/// this many layers
///
/// As this information is generally the same for each triangle this can be used
/// as a template to apply subdivision to any triangle
////////////////////////////////////////////////////////////////////////////////
void
gen_subdivide_triangle(
                        uint32_t numTriLayers,
                        float    *weights,
                        uint32_t *relPtIDs
                        )
{
  // Points
  uint32_t total = math::partial_sum( numTriLayers + 1 );
  for ( uint32_t i = 0; i < total; i++ )
  {
    subdiv_ptid_to_barycentric( i, numTriLayers, &weights[i * 3] );
  }

  // Triangles
  total = numTriLayers * numTriLayers;
  for ( uint32_t i = 0; i < total; i++ )
  {
    subdiv_tri_ptids( i, &relPtIDs[i * 3] );
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Apply the generalized triangle subdivision to a specific triangle ABC with
/// point IDs <p1,p2,p3> and write out the new triangles using new point IDs
/// as well as any new points designated as owned by this original triangle.
///
/// The triangle output is a series of triangles written out as:
/// [tOrig0_0, tOrig0_1, ... tOrig0_M, tOrig1_0, ... tOrig1_M, ..., tOrigN_M]
///   * N = num original tri
///   * M = numTriLayers^2
///   * The total global triangles is equal to numTriLayers^2 * num original tri
///   * triangle sets are grouped by M for each original tri
/// Thus, to find the global ID of a triangle, one only needs the original tri ID
/// and the relative tri ID within the subdivision. Each entry is a set of 3
/// global point IDs <p'1,p'2,p'3>
///
/// The point output is a series of points written out as:
/// [ (pOrig0, ... pOrigX), (pInt0, ... pIntY), (pExt0, ... pExtZ) ]
///   * X = num original points
///   * Y = num interior points per subdiv * num original tri
///   * Z = num edges * (numTriLayers - 1)
///   * Points are grouped logically by [original points, interior, exterior]
///     The use of () in the above diagram is illustrative and is not part of the
///     structure
/// Thus, to find the global ID of a point, one needs to know if it is an original
/// point, interior, or exterior. 
///   * For original points, use the same point ID in tri pt ID
///   * For interior points, one needs the original tri ID, start of (pInt0), and
///     relative interior ID.
///   * For exterior points, one needs the edge ID, start of (pExt0), and relative
///     ID along the edge.
/// Writing out of point cartesian values only occurs for new points. These are
/// (all the interior points) and (exterior points who lie on an edge of original
/// tri ABC where the point IDs p0 < p1, e.g. A < B for edge AB). Each entry is
/// a set of coordinates <x,y,z>.
///
/// The generalized triangle subdivision informs the relative tri of relative
/// point IDs, and for relative point IDs to the relative barycentric coords
///
/// This information in conjuction with determining respective global IDs allows
/// fully writing out a subdivided triangle and any owned points.
////////////////////////////////////////////////////////////////////////////////
void
apply_subdivision(
                  /// Input
                  uint32_t  numTriLayers,   ///< total number of triangle layers
                  float    *weights,        ///< generalized weights for subdivision
                  uint32_t *relPtIDs,       ///< generalized point IDs for subdivision
                  uint32_t  triID,          ///< original triangle ID
                  double   *triXYZ,         ///< original triangle XYZ values
                  uint32_t *triPtID,        ///< original triangle point IDs of ABC
                  uint32_t  interiorOffset, ///< global offset of points for interior pts to be generated
                  uint32_t  exteriorOffset, ///< global offset of points for exterior pts to be generated
                  int32_t  *edgeMap,        ///< map of original point IDs to edge info
                  /// Output
                  uint32_t *triangles,      ///< Full set of subdivided triangle mesh by unique point ID
                  double   *points          ///< Set of unique interior points in subdivided triangle mesh
                  )
{
  // Points
  uint32_t total = math::partial_sum( numTriLayers + 1 );
  uint32_t nInterior = ( numTriLayers * numTriLayers - 3 * numTriLayers + 2 ) / 2;
  uint32_t globalPtMap[total];

  DEBUG_OUT( "FOR TRIANGLE " << triID << std::endl );

  for ( uint32_t i = 0; i < total; i++ )
  {
    DEBUG_OUT( " FOR RELATIVE POINT ID " << i << std::endl );
    uint32_t globalPtID = 0;
    // Write out ONLY points that need to be written out
    // Skip original points and shared exterior that are not on owning edge
    PointInfo info = subdiv_ptid_info( i, numTriLayers );
    if ( info.ptType == 0 )
    {
      // Inform global subdiv pt id as original point
      globalPtID = triPtID[info.infoIdx];
    }
    else if ( info.ptType == -1 )
    {
      // Inform global subdiv pt id from tri id as interior point
      globalPtID = interiorOffset + triID * nInterior + info.infoIdx;

      // Write point out
      DEBUG_OUT( "  Writing to point ID " << globalPtID << std::endl );
      math::barycentric_to_spherical( triXYZ, &weights[i *3], &points[globalPtID * 3] );
    }
    else if ( info.ptType == 1 )
    {
      int32_t  *edgeMapAtPt = &edgeMap[triPtID[info.infoIdx] * MAX_EDGES_ON_PT * 2];
      uint32_t  edgeID      = 0;
#ifdef DEBUG
      std::cout << "  Tri [" << triPtID[0] << "," 
                             << triPtID[1] << ", "
                             << triPtID[2] << "]" << std::endl;
      std::cout << "  Looking for point id " << triPtID[(info.infoIdx + 1) % 3]
                << " in &edgeMap[" <<  triPtID[info.infoIdx] << "]:" << std::endl;
      for ( uint8_t e = 0; e < MAX_EDGES_ON_PT; e++ )
      {
        if ( edgeMapAtPt[e * 2] == -1 ) { continue; }
        std::cout << "    point: " << edgeMapAtPt[e * 2]
                  << " edge ID: " << edgeMapAtPt[e * 2 + 1] << std::endl;
      }
#endif

      for ( uint8_t e = 0; e < MAX_EDGES_ON_PT; e++ )
      {
        if ( edgeMapAtPt[e * 2] == -1 ) { continue; }
        if ( edgeMapAtPt[e * 2] == triPtID[(info.infoIdx + 1) % 3] )
        {
          edgeID = edgeMapAtPt[e * 2 + 1];
          break;
        }
      }

      // Find the local index along the edge when traversing ccw
      uint32_t localEdgeIdx = 0;
      if      ( info.infoIdx == 0 ) { localEdgeIdx = info.ptLayer - 1; }
      else if ( info.infoIdx == 1 ) { localEdgeIdx = info.ptID - ( total - numTriLayers ); }
      else                          { localEdgeIdx = (numTriLayers - 1) - info.ptIdx; }

      // If the owned edge is registered as cw for us (i.e we do not own the edge
      // and are traversing the indices in reverse), flip the indexing
      bool owned = triPtID[info.infoIdx] < triPtID[(info.infoIdx + 1) % 3];
      if ( !owned ) { localEdgeIdx = ( numTriLayers - 2 ) - localEdgeIdx; }

#ifdef DEBUG
      std::cout << "    Local pt id " << info.ptID 
                << "[idx:" << info.ptIdx << "/layer:" << info.ptLayer << "]"
                << " lies on edge " << edgeID << "[" << info.infoIdx << "]"
                << " local edge idx " << localEdgeIdx << std::endl;
#endif

      // Inform global subdiv pt id from edge map as exterior point
      globalPtID = exteriorOffset + edgeID * ( numTriLayers - 1 ) + localEdgeIdx; 

      // Write point out
      if ( owned )
      {
        DEBUG_OUT( "      Writing to point ID " << globalPtID << std::endl );
        math::barycentric_to_spherical( triXYZ, &weights[i *3], &points[globalPtID * 3] );
      }
    }

    globalPtMap[i] = globalPtID;
  }

  // Triangles
  total = numTriLayers * numTriLayers;
  for ( uint32_t i = 0; i < total; i++ )
  {
    // Write out all triangles
    uint32_t globalTriID = triID * total + i;
    triangles[globalTriID * 3 + 0] = globalPtMap[relPtIDs[i * 3 + 0]];
    triangles[globalTriID * 3 + 1] = globalPtMap[relPtIDs[i * 3 + 1]];
    triangles[globalTriID * 3 + 2] = globalPtMap[relPtIDs[i * 3 + 2]];
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Generate a point to edge ID map and return the number of edges found
///
/// The map iterates through all triangles and catalogues edges on a vertex.
/// To guarantee uniqueness, an edge ID is only added when p1 id < p2 id when
/// traversing edges in the triangle. If true, the edge ID is added to the map
/// for p1 and p2. Triangles must be consistent in traversal order, that is either
/// clockwise or counterclockwise.
///
/// The map is laid out as:
/// [(pA eX, pB eY, ... )_0, (pA' eX', pB' eY', ... )_1 ]
///   * Each point ()_# has MAX_EDGES_ON_PT * 2 entries consisting of alternating
///     point id-edge id pairs, where # = point ID in question
///   * Within the ()_# grouping, p* notes the edge from point # - point p*
///     REGARDLESS of cw/ccw directionality
///   * The correspoinding edge ID within the pair p*-e* notes the unique edge ID
///     when edge point # - p* is correct
///
/// To find the edge ID of any two points one can use the map by first indexing
/// the group ()_# for the point, and then finding the p* ID within the pairs
/// corresponding to the next point.
////////////////////////////////////////////////////////////////////////////////
uint32_t
gen_edge_map(
              uint32_t *triangles,
              uint32_t  nTri,
              int32_t  *edgeMap
              )
{
  uint32_t nEdge = 0;
  std::map< uint32_t, uint8_t > currIdx;

  // Loop over all triangle [pt1, pt2, pt3]
  for ( uint32_t i = 0; i < nTri; i++ )
  {
    uint32_t *tri = &triangles[i * 3];
    // Loop over all edges of this triangle
    for ( uint8_t e = 0; e < 3; e++ )
    {
      uint8_t next = (e + 1) % 3;
      // This is a unique edge
      if ( tri[e] < tri[next] )
      {
        // Fill out map for originating point
        edgeMap[tri[e] * MAX_EDGES_ON_PT * 2 + currIdx[tri[e]] * 2 + 0] = tri[next];
        edgeMap[tri[e] * MAX_EDGES_ON_PT * 2 + currIdx[tri[e]] * 2 + 1] = nEdge;
        currIdx[tri[e]]++;

        // Fill out map for leading point
        edgeMap[tri[next] * MAX_EDGES_ON_PT * 2 + currIdx[tri[next]] * 2 + 0] = tri[e];
        edgeMap[tri[next] * MAX_EDGES_ON_PT * 2 + currIdx[tri[next]] * 2 + 1] = nEdge;
        currIdx[tri[next]]++;

        nEdge++;
      }
    }
  }
  return nEdge;
}

} // namespace tri_subdiv


void
help( void )
{
  std::cout << "Usage: [cmd] <input mesh> <output mesh> <number of subdivisions>" << std::endl
            << "  Subdivide a spherical JIGSAW mesh using an integer number of" << std::endl
            << "  subdivisions along an edge, e.g. A subdivision of 2 results in:" << std::endl
            << "" << std::endl
            << "             * "            << std::endl
            << "            / \\ "          << std::endl
            << "           /   \\ "         << std::endl
            << "          *-----* "         << std::endl
            << "         / \\   / \\ "      << std::endl
            << "        /   \\ /   \\ "     << std::endl
            << "       *-----*-----* "      << std::endl
            << "      / \\   / \\   / \\ "  << std::endl
            << "     /   \\ /   \\ /   \\ " << std::endl
            << "    *-----*-----*-----* "   << std::endl
            << "  Where two cuts/points are applied along each edge" << std::endl
            << "  of every triangle." << std::endl
            << "" << std::endl
            << "  The change in triangle area is roughly equal to area / (n+1)^2" << std::endl;
}

int32_t main( int32_t argc, char const *argv[] )
{
  if ( argc != 4 )
  {
    std::cerr << "Error: Incorrect number of arguments!" << std::endl;
    help();
    return -1;
  }

  std::string filename( argv[1] );
  std::string output( argv[2] );
  uint32_t numTriLayers = std::stoi( std::string( argv[3] ) ) + 1;

  JMesh meshIn;
  meshIn.read( filename );

  if ( meshIn.dims_ != 3 )
  {
    std::cerr << "Error: subdivision only works on 3D spherical meshes right now" << std::endl;
    throw 1;
  }

  // Make into unit vectors
  double norm = math::norm( &meshIn.points_[0] );
  for ( uint32_t i = 0; i < meshIn.point_count(); i++ )
  {
    math::normalize( &meshIn.points_[i * 3] );
  }

  std::cout << "Subdividing mesh in " << filename
            << " by " << numTriLayers - 1
            << " subdivision(s)" << std::endl; 

  std::vector< int32_t > edgeMap( meshIn.point_count() * MAX_EDGES_ON_PT * 2, -1 );
  uint32_t nEdge = tri_subdiv::gen_edge_map(
                                            meshIn.triangles_.data(),
                                            meshIn.triangle_count(),
                                            edgeMap.data()
                                            );

#ifdef DEBUG
  for ( uint32_t i = 0; i < meshIn.point_count(); i++ )
  {
    std::cout << "Point " << i << " Edge Map : " << std::endl;
    int32_t *edgeMapAtPt = &edgeMap[ i * MAX_EDGES_ON_PT * 2 ];
    for ( uint8_t j = 0; j < MAX_EDGES_ON_PT; j++ )
    {
      if ( edgeMapAtPt[j * 2] == -1 ) { continue; }
      std::cout << "  p : " << edgeMapAtPt[j * 2 + 0] << " e : " << edgeMapAtPt[j * 2 + 1] << std::endl;
    } 
  }
#endif

  std::cout << "Found " << nEdge << " egdes" << std::endl;
  uint32_t interiorPtPerTri = ( numTriLayers * numTriLayers - 3 * numTriLayers + 2 ) / 2;
  uint32_t totalInterior    = interiorPtPerTri * meshIn.triangle_count();
  uint32_t totalPts         = meshIn.point_count()                // original points
                            + totalInterior                       // interior points
                            + ( (numTriLayers - 1) * nEdge );  // exterior points
  uint32_t totalTri         = meshIn.triangle_count() * numTriLayers * numTriLayers;

  std::cout << "Will subdivide mesh into:" << std::endl
            << "     Points: " << totalPts << std::endl
            << "  Triangles: " << totalTri << std::endl;


  // Prepare subdivision info
  std::vector< float > weights( (interiorPtPerTri + 3 * numTriLayers) * 3 );
  std::vector< uint32_t > relPtIDs( numTriLayers * numTriLayers * 3 );
  tri_subdiv::gen_subdivide_triangle( numTriLayers, weights.data(), relPtIDs.data() );

#ifdef DEBUG
  std::cout << "Relative Triangle subdivision info:" << std::endl;
  std::cout << "  Using barycentric weights:" << std::endl;
  for ( uint32_t i = 0; i < math::partial_sum( numTriLayers + 1 ); i++ )
  {
    float *weightsAt = &weights[i * 3];
    std::cout << "    Point ID " << i << "[" 
              << weightsAt[0] << ","
              << weightsAt[1] << ", "
              << weightsAt[2] << "]" << std::endl;
  }
  std::cout << "  Using relative points for triangles:" << std::endl;
  for ( uint32_t i = 0; i < numTriLayers * numTriLayers; i++ )
  {
    uint32_t *relPtIDAt = &relPtIDs[i * 3];
    std::cout << "    Triangle ID " << i << "[" 
              << relPtIDAt[0] << ","
              << relPtIDAt[1] << ", "
              << relPtIDAt[2] << "]" << std::endl;
  }
#endif

  JMesh meshOut;
  meshOut.points_.resize( totalPts * 3 );
  meshOut.triangles_.resize( totalTri * 3 );
  // Get original points
  std::memcpy(
              meshOut.points_.data(),
              meshIn.points_.data(),
              sizeof(meshIn.points_[0]) * meshIn.points_.size()
              );

  std::cout << "Starting subdivision..." << std::endl;

  for ( uint32_t i = 0; i < meshIn.triangle_count(); i++ )
  {
    double   tri[9];
    // Copy points ABC into local buffer
    std::memcpy( &tri[0], &meshIn.points_[meshIn.triangles_[i * 3 + 0] * 3], sizeof(tri[0]) * 3 );
    std::memcpy( &tri[3], &meshIn.points_[meshIn.triangles_[i * 3 + 1] * 3], sizeof(tri[0]) * 3 );
    std::memcpy( &tri[6], &meshIn.points_[meshIn.triangles_[i * 3 + 2] * 3], sizeof(tri[0]) * 3 );

    // Subdivide
    tri_subdiv::apply_subdivision(
                                  numTriLayers,
                                  weights.data(),
                                  relPtIDs.data(),
                                  i,
                                  tri,
                                  &meshIn.triangles_[i * 3],
                                  meshIn.point_count(),
                                  totalInterior + meshIn.point_count(),
                                  edgeMap.data(),
                                  meshOut.triangles_.data(),
                                  meshOut.points_.data()
                                  );
  }
  std::cout << "Subdivision complete" << std::endl;

  for ( double &p : meshOut.points_ )
  {
    p *= norm;
  }
  meshOut.dims_ = meshIn.dims_;
  meshOut.meshID_ = meshIn.meshID_;
  meshOut.meshType_ = meshIn.meshType_;
  meshOut.write( output );
  std::cout << "SUCCESS" << std::endl;
  return 0;
}
