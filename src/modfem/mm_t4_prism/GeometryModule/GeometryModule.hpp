#ifndef GEOMETRYMODULE_HPP_INCLUDED
#define GEOMETRYMODULE_HPP_INCLUDED

#include "../Common.h"
#include "../MeshRead/IMeshReader.h"
#include "../Field.hpp"

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


class NoGeometryModule
{
public:
	static NoGeometryModule & Instance()
	{
		static NoGeometryModule _instance;
		return _instance;
	}

    ~NoGeometryModule()
    {
    }

    void    Free()
    {
    }

    int     size() const { return 0; }

    int    ReadPoints(MeshRead::IMeshReader & reader)
    {
		return 0;
    }

    void    GetPointAt(IN OUT Tval coords[3])
    {}

    bool    test()
    {
        return true;
    }
};

template< typename Tval = double>
class FirstGeometryModule
{
    typedef Memory::Field<3,Tval> TCoord;
public:

    static FirstGeometryModule & Instance()
	{
		static FirstGeometryModule _instance;
		return _instance;
	}

    ~FirstGeometryModule()
    {
        Free();
    }

    void    Free()
    {
	  safeDeleteArray(_points);
        _size = 0;
    }

    int     size() const { return _size; }

    int    ReadPoints(MeshRead::IMeshReader & reader)
    {
        if(reader.Init())
        {
            _size = reader.GetVerticesCount();
			std::cout << " -> Reading " << _size  << " points";
			
            delete[] _points;

            if(_size > 0)
            {
                _points = new TCoord[_size];

                bool ok(true);
				int i(0);
				for(; i < _size && ok; ++i)
                {
                    ok = reader.GetNextVertex(_points[i].getArray());
                }

                if(i < _size)
                {
                    throw "GeometryModule::ReadPoints: error while reading: not all points readed!";
                }

            }
            //TODO: check readed data correctness
            //TODO: presort readed data
        }
        return _size;
    }

    void    GetPointAt(IN OUT Tval coords[3])
    {
        findNearest(coords[0],coords[1],coords[2]);
    }

    bool    test()
    {
        Tval v[] ={1,1,1};
        Tval v2[] ={1,1,1};
        GetPointAt(v);
        return v != v2;
    }

private:
    int         _size;
    TCoord*     _points;

protected:
    FirstGeometryModule() : _size(0), _points(NULL) { }
    FirstGeometryModule(const FirstGeometryModule & other);

    void    findNearest(Tval & x, Tval & y, Tval & z)
    {
        if(_size > 0)
		  {
			//std::cout << "\n("<<x<<","<<y<<","<<z<<")";
            TCoord last(_points[0]);
            Tval  dist = ( (x-last[0])*(x-last[0])+(y-last[1])*(y-last[1])+(z-last[2])*(z-last[2]) );
            for(int i(1); i < _size; ++i)
            {
                // if distance is smaller, change coords
			  const Tval newdist= ( (x-_points[i][0])*(x-_points[i][0])
							  +(y-_points[i][1])*(y-_points[i][1])
							  +(z-_points[i][2])*(z-_points[i][2]) );
			  if(dist >  newdist) 
                {
                    last[0] = _points[i][0];
                    last[1] = _points[i][1];
                    last[2] = _points[i][2];
                    dist = newdist;
                }
            }

            x = last[0];
            y = last[1];
            z = last[2];
			//std::cout << " detailed to ("<<x<<","<<y<<","<<z<<")";
        }
    }

};

template< typename Tval = double>
class FileGeometryModule
{
    typedef Memory::Field<3,Tval> TCoord;
public:

    static FileGeometryModule & Instance()
	{
		static FileGeometryModule _instance;
		return _instance;
	}

    ~FileGeometryModule()
    {
        Free();
    }

    void    Free()
    {
        //delete [] _points;
        _size = 0;
    }

    int     size() const { return _size; }

    int    ReadPoints(MeshRead::IMeshReader & reader)
    {
        _reader = & reader;

        return _reader == NULL ? 0  : 1;
    }

    void    GetPointAt(IN OUT Tval coords[3])
    {
        findNearest(coords[0],coords[1],coords[2]);
    }

    bool    test()
    {
        Tval v[] ={1,1,1};
        Tval v2[] ={1,1,1};
        GetPointAt(v);
        return v != v2;
    }

private:
    int         _size;
    //TCoord*     _points;
    MeshRead::IMeshReader * _reader;

protected:
    FileGeometryModule() : _size(0) { }
    FileGeometryModule(const FileGeometryModule & other);

    void    findNearest(Tval & x, Tval & y, Tval & z)
    {
        if(_reader != NULL)
        {

            if(_reader->Init())
            {
                _size = _reader->GetVerticesCount();

                //delete[] _points;

//                if(_size > 0)
//                {
//                    //_points = new TCoord[_size];
//
//                    bool ok(true);
//                    for(int i(0); i < _size && ok; ++i)
//                    {
//                        ok = reader->GetNextVertex(_points[i].getArray());
//                    }
//
//                    if(!ok)
//                    {
//                        throw "GeometryModule::ReadPoints: error while reading.";
//                    }


            TCoord last;
            _reader->GetNextVertex(last.getArray());
            Tval  dist = ( (x-last[0])*(x-last[0])+(y-last[1])*(y-last[1])+(z-last[2])*(z-last[2]) );
            TCoord current;
            for(int i(1); i < _size; ++i)
            {
                _reader->GetNextVertex(current.getArray());
                // if distance is smaller, change coords
                Tval newdist= ( (x-current[0])*(x-current[0])+(y-current[1])*(y-current[1])+(z-current[2])*(z-current[2]) );
                if(dist >  newdist)
                {
                    last[0] = current[0];
                    last[1] = current[1];
                    last[2] = current[2];
                    dist = newdist;
                }
            }

            _reader->Free();

            x = last[0];
            y = last[1];
            z = last[2];
            }
                //TODO: check readed data correctness
                //TODO: presort readed data
            }

    }

};
/**  @} */
#endif // GEOMETRYMODULE_HPP_INCLUDED
