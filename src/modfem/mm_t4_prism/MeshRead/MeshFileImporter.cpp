
#include "MeshFileImporter.h"
#include <cassert>
#include <iostream>

#include "../MeshModule/hHybridMesh.h"

namespace MeshRead
{

MeshFileImporter::MeshFileImporter() : _vert_count(0),_vert_readed(0),_elem_count(0), _elem_readed(0)
{
}


MeshFileImporter::MeshFileImporter(const std::string & file_name) : _vert_count(0),_vert_readed(0),_elem_count(0), _elem_readed(0)
{
	_file_name = file_name;
}

MeshFileImporter::~MeshFileImporter(void)
{
	Free();
}
  
  std::string MeshFileImporter::Name() const
  {
	return _file_name;
  }
  
bool	MeshFileImporter::Init(const std::string & file_name)
{
    _file_name = file_name;
    return MeshFileImporter::Init();
}

bool	MeshFileImporter::Init()
{
  Free();
	if(false == _file_name.empty())
	{
		_file.open(_file_name.c_str(),std::ios_base::in);
	}

	mmv_out_stream << "\nMeshModule::opening file: " << _file_name << " is_open=" << _file.is_open();

	return _file.is_open();
}

void	MeshFileImporter::Free()
{
    if(_file.is_open())
    {
        _file.close();
    }

    _vert_count = 0;
    _vert_readed = 0;
    _elem_count = 0;
    _elem_readed = 0;

}

int	MeshFileImporter::getPos()
{
    char c=_file.peek();
    if(c=='\n' || c=='\r') {
	_file.ignore(3,'\n');
	_file.putback(' ');
    }
    return static_cast<int>(_file.tellg());
}

void	MeshFileImporter::setPos(const int newPos) 
{
    _file.seekg(newPos);
}
  
/** returns number of mesh dimensions: 1D/2D/3D
*/
  int MeshFileImporter::GetCoordinatesDimension() const
  {return 3;}

/** /return no. of vertices in initial mesh.
*/
  int MeshFileImporter::GetVerticesCount()
  {
	return _vert_count;
  }

/**	/param  coords array for coordinates of vertex
/return true if there is next vertex and false if there isn't.
*/
bool MeshFileImporter::GetNextVertex(double coords[])
{ return false; };

/** /return no. of elements in inintial mesh
*/
  void MeshFileImporter::GetElementCount(int element_type[])
  {};

	/** /param  vertices array for current element vertices numbers
	/param  neighbours array for current element neighoubrs(other elements) numbers.
	/param	bc number of boundary condition at this element
	/return true if there is next element and false if this is last element.
	*/
  bool MeshFileImporter::GetNextElement(Tind vertices[], Tind neighbours[], Tind faces[], Tind & element_type , Tind & father, int8_t  & material, int8_t & ref)
  {return false;}
	/** /param bc array for boundary condition parameters
	bc[0] - Dirichlet BC
	bc[1] - Neumann BC
	bc[2] - Cauchy BC par1
	bc[3] - Cauchy BC par2
	*/
  bool MeshFileImporter::GetBoundaryConditions(double ** bc, int & bcCount)
  {return false;}

  
} // namespace

