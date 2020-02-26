#ifndef MESH_FILE_IMPORTER_H
#define MESH_FILE_IMPORTER_H

#include "IMeshReader.h"

#include <fstream>
#include <string>

namespace MeshRead
{

class MeshFileImporter :
	public MeshRead::IMeshReader
{
public:
	MeshFileImporter(void);
	MeshFileImporter(const std::string & file_name);
	virtual ~MeshFileImporter(void);

	/** /return true if MeshBuilder initialized successfully, otherwise returns false.
	*/
	bool Init(const std::string & file_name);

	/** /return true if MeshBuilder initialized successfully, otherwise returns false.
	*/
	bool Init();
	
	/** /return name of mesh source
	 */
	std::string Name() const;
	
	/** Frees resources used by builder.
	*/
	void Free();

	bool	doRead(hHybridMesh * mesh)=0;

	/** returns number of mesh dimensions: 1D/2D/3D
	*/
	int GetCoordinatesDimension() const ;

	/** /return no. of vertices in initial mesh.
	*/
	int GetVerticesCount() ;

	/**	/param  coords array for coordinates of vertex
	/return true if there is next vertex and false if there isn't.
	*/
	bool GetNextVertex(double coords[]) ;

	/** /return no. of elements in inintial mesh
	*/
	void GetElementCount(int element_type[]) ;

	/** /param  vertices array for current element vertices numbers
	/param  neighbours array for current element neighoubrs(other elements) numbers.
	/param	bc number of boundary condition at this element
	/return true if there is next element and false if this is last element.
	*/
	bool GetNextElement(Tind vertices[], Tind neighbours[], Tind faces[], Tind & element_type , Tind & father, int8_t  & material, int8_t & ref) ;

	/** /param bc array for boundary condition parameters
	bc[0] - Dirichlet BC
	bc[1] - Neumann BC
	bc[2] - Cauchy BC par1
	bc[3] - Cauchy BC par2
	*/
	bool GetBoundaryConditions(double ** bc, int & bcCount) ;

protected:
	int	getPos();
	void	setPos(const int newPos);
	std::string		_file_name;
	std::fstream	_file;
	int 			_vert_count,
                    _vert_readed,
		    		_elem_count,
		    		_elem_readed,
					gridPos_,
					elemPos_,
					facePos_;

};

}// namespace

#endif //MESH_FILE_IMPORTER_H
