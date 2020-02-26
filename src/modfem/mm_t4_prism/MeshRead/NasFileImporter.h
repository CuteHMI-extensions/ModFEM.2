#ifndef NASFILEIMPORTER_H
#define NASFILEIMPORTER_H

#include <vector>

#include "MeshFileImporter.h"

namespace MeshRead
{

class NasFileImporter
    : public MeshFileImporter
{
    public:
        NasFileImporter();
        NasFileImporter(const std::string & file_name);

        /** /return true if MeshBuilder initialized successfully, otherwise returns false.
        */
        bool Init(const std::string & file_name);

        /** /return true if MeshBuilder initialized successfully, otherwise returns false.
        */
        bool Init();

		bool	doRead(hHybridMesh * mesh);

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

        /** /return no. of elements in inintial mehs
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

		bool LineParse(char const* str, std::vector<double>& v);
		bool VertexLineParse(char const* str, std::vector<double>& v);

		//int GetEdgesCount();
		void GetFacesCount(int type_count[]);
		bool GetNextFace(Tind edges[],Tind & face_type, int8_t & bc, Tind neigh[]);
		//bool GetNextEdge(Tind vertices[],Tind & edge_type);

    protected:
		        int     _face_count,
                _face_readed,
                _edge_count,
                _edge_readed,
				_quads,
				_prisms;


};

} // namespace MeshRead

#endif // NASFILEIMPORTER_H
