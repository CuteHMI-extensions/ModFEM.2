
#include "BinaryFileReader.h"
#include "../MeshModule/hHybridMesh.h"

namespace MeshRead {
  
  BinaryFileReader::BinaryFileReader()
	:_file_name("not initialized!")
  {
  }

  BinaryFileReader::BinaryFileReader(const std::string & fileName)
	: _file_name(fileName)
  {
  }

  BinaryFileReader::~BinaryFileReader()
  {
	Free();
  }
  
bool BinaryFileReader::Init() {
    _file.open(_file_name.c_str(), std::ios::in | std::ios::binary);
    return _file.is_open();
}

bool BinaryFileReader::Init(const std::string & s) {
    _file_name = s;
    _file.open(_file_name.c_str(), std::ios::in | std::ios::binary);
    return _file.is_open();
}

void BinaryFileReader::Free() 
{
    _file.close();
}

  std::string BinaryFileReader::Name() const
{
  return std::string("BinaryFileReader")+_file_name;
}

bool BinaryFileReader::doRead(hHybridMesh * m)
{
  if(_file.good()) {
	//char separator('E');
	//int meshId(0);
	//_file >> meshId;
	_file >> m->name_ ; //>> separator;
	  //assert(separator == ' ');
	m->vertices_.read(_file);

	//_file >> separator;
	//assert(separator == ' ');
	m->edges_.read(_file);

	//_file >> separator;
	//assert(separator == ' ');
	m->faces_.read(_file);

	//_file >> separator;
	//assert(separator == ' ');
	m->elements_.read(_file);

	_file >> m->gen_ >> m->maxGen_ >> m->maxGenDiff_;

    //hObj::myMesh = m;
	m->test();
  }
  return true;
}

int BinaryFileReader::GetCoordinatesDimension() const 
{
    return 3;
}

int BinaryFileReader::GetVerticesCount()
{
    return 0;
}

bool BinaryFileReader::GetNextVertex(double *)
{
    return false;
}

void BinaryFileReader::GetElementCount(int*)
{
}

bool BinaryFileReader::GetNextElement(Tind*, Tind*,Tind*,Tind&,Tind&,int8_t&,int8_t&)
{
    return false;
}

bool BinaryFileReader::GetBoundaryConditions(double **, int&)
{
    return false;
}
  
  bool BinaryFileReader::needsSetup() const
  {
	return false;
  }
  
};
