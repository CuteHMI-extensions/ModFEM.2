#include "SimplePointReader.h"

namespace MeshRead {
  
  SimplePointReader::SimplePointReader() {}
  SimplePointReader::SimplePointReader(const std::string &  file_name)
  {
	_file_name=file_name;
  }
  
bool SimplePointReader::Init(const std::string & file_name)
{
  _file_name = file_name;
  return Init();

}

  
  bool SimplePointReader::Init()
  {
	if(MeshFileImporter::Init()){
  	char tmp[255]={0};
	_vert_count=0;
	_file.clear();
	do {
	  _file.getline(tmp,255);
	  ++_vert_count;
	}while(_file.good());
	_file.clear();
	_file.seekg(0);
	}
  return _file.is_open();	  
	
  }

  bool SimplePointReader::doRead(hHybridMesh * mesh)
  {
	return false;
  }
  
bool SimplePointReader::GetNextVertex(double coords[])
{
  double tmp;
  if(_vert_readed < _vert_count && _file.good()) {
	_file >> tmp //temporary - some index or something
		  >> coords[0] >> coords[1] >> coords[2];
  }
  return _file.good();
}
  

	        /** returns number of mesh dimensions: 1D/2D/3D
        */
    int SimplePointReader::GetCoordinatesDimension() const {
  return 3;
  }

        /** /return no. of vertices in initial mesh.
        */
    int SimplePointReader::GetVerticesCount()
  {
  	return _vert_count;
  }


        /** /return no. of elements in inintial mehs
        */
  void SimplePointReader::GetElementCount(int element_type[])
  {
  }

        /** /param  vertices array for current element vertices numbers
        /param  neighbours array for current element neighoubrs(other elements) numbers.
        /param	bc number of boundary condition at this element
        /return true if there is next element and false if this is last element.
        */
  bool SimplePointReader::GetNextElement(Tind vertices[], Tind neighbours[], Tind faces[], Tind & element_type , Tind & father, int8_t  & material, int8_t & ref)
  {
	return false;
  }

        /** /param bc array for boundary condition parameters
        bc[0] - Dirichlet BC
        bc[1] - Neumann BC
        bc[2] - Cauchy BC par1
        bc[3] - Cauchy BC par2
        */
  bool SimplePointReader::GetBoundaryConditions(double ** bc, int & bcCount)
  {
	return false;
  }

		//		bool LineParse(char const* str, std::vector<double>& v);
		//		bool VertexLineParse(char const* str, std::vector<double>& v);

		//int GetEdgesCount();
  //		void SimplePointReader::GetFacesCount(int type_count[]);
  //		bool SimplePointReader::GetNextFace(Tind edges[],Tind & face_type, int8_t & bc, Tind neigh[]);
		//bool GetNextEdge(Tind vertices[],Tind & edge_type);


  
}
