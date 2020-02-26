#ifndef DMPFILEIMPORTER_H
#define DMPFILEIMPORTER_H

#include "MeshFileImporter.h"

namespace MeshRead
{

class DmpFileImporter : public MeshRead::MeshFileImporter
{
    public:
        DmpFileImporter(const std::string & fileName);
        DmpFileImporter();
        ~DmpFileImporter();

        bool  Init();

		bool doRead(hHybridMesh * mesh);

        /** /return no. of vertices in initial mesh.
	*/
	int GetVerticesCount() ;

	/**	/param  coords array for coordinates of vertex
	/return true if there is next vertex and false if there isn't.
	*/
	bool GetNextVertex(double coords[]) ;

	/** /return no. of elements in inintial mesh
	*/
	void GetElementCount(int type_count[]) ;

	/** /param  vertices array for current element vertices numbers
	/param  neighbours array for current element neighoubrs(other elements) numbers.
	/param	bc number of boundary condition at this element
	/return true if there is next element and false if this is last element.
	*/
	bool GetNextElement(Tind vertices[], Tind neighbours[], Tind faces[], Tind & element_type , Tind & father, int8_t  & material, int8_t & ref);

	/** /param bc array for boundary condition parameters
	bc[0] - Dirichlet BC
	bc[1] - Neumann BC
	bc[2] - Cauchy BC par1
	bc[3] - Cauchy BC par2
	*/

	int GetCoordinatesDimension() const { return 3; }
    bool GetBoundaryConditions(double **, int &) { return false; }

    int GetEdgesCount();
    void GetFacesCount(int type_count[]);
    bool GetNextFace(Tind edges[],Tind & face_type, int8_t & bc, Tind neigh[]);
    bool GetNextEdge(Tind vertices[],Tind & edge_type);


    protected:
        int     _face_count,
                _face_readed,
                _edge_count,
                _edge_readed;
    private:
};

}

#endif // DMPFILEIMPORTER_H
