#ifndef __CARBUNCLE_LEVELSET_READER_HPP__
#define __CARBUNCLE_LEVELSET_READER_HPP__

#include <string>
#include <fluids/levelset_fluid3.h>

namespace Carbuncle {

class LevelSetReader {
public:
  LevelSetReader(std::string filename);

  /**
   * @brief this method reads a Signal Distance Field file, as described in
   * https://github.com/christopherbatty/SDFGen:
   *
   * The input file format is:<ni> <nj> <nk>
   *  <origin_x> <origin_y> <origin_z>
   *  <dx>
   *  <value_1> <value_2> <value_3> [...]
   *
   * * (ni,nj,nk) are the integer dimensions of the resulting distance field.
   * * (origin_x,origin_y,origin_z) is the 3D position of the grid origin.
   * * <dx> is the grid spacing.
   * * <value_n> are the signed distance data values, in ascending order of i,
   then j, then k. The output filename will match that of the input, with the
   OBJ suffix replaced with SDF.
   * @param levelset
   */
  void read(Leviathan::LevelSetFluid3 &levelset);

private:
  std::string _filename;
};

} // namespace Carbuncle

#endif