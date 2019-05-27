#include <structures/levelset2.h>
#include <structures/levelset3.h>
#include <string>

class LevelSetFluid2 : public Ramuh::LevelSet2 {
public:
  LevelSetFluid2();

  void run();

  void loadConfiguration(std::string filename);

  void writeVelocityField();

  std::string &getFolderName();

  std::string &getDataFolderName();

  int getVelocityOrder();

  int getLevelsetOrder();

  int getFramesNumber();

private:
  int _resolution, _nFrames;
  int _velocityAdvectionOrder;
  int _levelsetAdvectionOrder;
  std::string _objFolderName, _dataFolderName;
  bool _isPressure2nd;
};

class LevelSetFluid3 : public Ramuh::LevelSet3 {
public:
  LevelSetFluid3();

  void run();

  void loadConfiguration(std::string filename);

  void writeVelocityField();

  std::string &getFolderName();

  std::string &getDataFolderName();

  int getVelocityOrder();

  int getLevelsetOrder();

  int getFramesNumber();

private:
  int _resolution, _nFrames;
  int _velocityAdvectionOrder;
  int _levelsetAdvectionOrder;
  std::string _objFolderName, _dataFolderName;
  bool _isPressure2nd;
};
