#ifndef ARRAYPOOL_HPP_INCLUDED
#define ARRAYPOOL_HPP_INCLUDED

#include <iostream>
//#include <map>

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


#include "../Common.h"
//#include "VtsSqId.hpp"

template< class  BaseT>
class ArrayPool
{
public:
  //typedef std::map< VtsSqId< TnVerts >, uTind> hashPosByVts;
  //typedef hashPosByVts::value_type vtsPos;
  //typedef hashPosByVts::iterator vts2posIter;
  //typedef hashPosByVts::const_iterator vts2posCIter;
  
    ArrayPool(): pool_(NULL), size_(0), uniqueID_(0), sizeChange_(0), capacity_(0), nDeleted_(0)
	{
		clear();
	}
  
  
  ~ArrayPool()
  {
	safeDeleteArray(pool_);
  }
  
	BaseT & first(){ return *reinterpret_cast<BaseT*>(pool_); }
	const BaseT & first() const { return *reinterpret_cast<BaseT*>(pool_); }
	const BaseT*  end() const { return reinterpret_cast<const BaseT*>(pool_+size_+1); } // first after last
	uTind last() const { return static_cast<uTind>(size_); }
    size_t size() const { return size_-nDeleted_; }
	size_t capacity() const { return capacity_; }
	bool   empty() const { return size_==0; }
	void	clear()
	{
		safeDeleteArray(pool_);
		size_=0;
		size_=0;
		capacity_=0;
        nDeleted_=0;
	}
	BaseT & operator[](const uTind pos) { return reinterpret_cast<BaseT*>(pool_)[pos-FIRST]; }
	const BaseT & operator[](const uTind pos) const { return reinterpret_cast<const BaseT*>(pool_)[pos-FIRST]; }

	BaseT & at(const uTind pos) { assert(pos > 0 && pos <= size_); return reinterpret_cast<BaseT*>(pool_)[pos-FIRST];}
	const BaseT & at(const uTind pos) const { assert(pos > 0 && pos <= size_); return reinterpret_cast<const BaseT*>(pool_)[pos-FIRST];}

	BaseT&	getById(const ID id)
	{
		assert( at(BaseT::posFromId(id)).type_== BaseT::typeFromId(id));
		assert( at(BaseT::posFromId(id)).id_== id);
		return at(BaseT::posFromId(id));
	}

	const BaseT & getById(const ID id) const
	{
		assert( at(BaseT::posFromId(id)).type_== BaseT::typeFromId(id));
		assert( at(BaseT::posFromId(id)).id_== id);
		return at(BaseT::posFromId(id));
	}

  BaseT*  newObj(const uTind verts[])
  {
 		assert( size_+1 <= capacity_);
		assert(usedMemory() < totalMemory());
		assert(uniqueID_ == size_);	// checking if all objcts are created by newObj(..)
		
		BaseT & obj(* new (pool_+(sizeof(BaseT)*size_++)) BaseT(verts, uniqueID() ) );

		
		//obj.myMesh->registerHObj(obj.type_,verts,obj.pos_);
		//vts2pos_.insert(vtsPos(verts,obj.pos_));
		
		return  &obj;
	}

	size_t	totalMemory() const { return sizeof(BaseT)*capacity_; }
	size_t	usedMemory() const {return sizeof(BaseT)*size_;}

	void	reserve(const size_t newcapacity, bool force=false)
	{
		if(force || newcapacity > capacity_)
		{
			BYTE * newPool=new BYTE[newcapacity*sizeof(BaseT)];
			memset(newPool,0,newcapacity*sizeof(BaseT));
			if(size_ > 0)
			{
				memcpy(newPool,pool_,std::min(size_,newcapacity)*sizeof(BaseT));
				safeDeleteArray(pool_);
			}
			pool_ = newPool;
			capacity_ = newcapacity;
		}
	}

	void    requestChange(const Tind sizeChange)
	{
		sizeChange_+=sizeChange;
	}

	uTind memoryNeeded() const  // in bytes
	{
		return static_cast<uTind>(usedMemory() + sizeChange_*sizeof(BaseT));
	}

	void    adjust()
	{
		/// Pre: +all markings have beed done (correctly).
		///     + pool is not empty
        //if (memoryNeeded() > usedMemory())
		{
            reserve(size_+sizeChange_,true);

            for(BaseT *ptr(&first()); ptr != end(); ++ptr) {
                if(ptr->nMyClassSons_ == BaseT::delMark) {
                    memset(ptr,0,sizeof(BaseT));
                    ++nDeleted_;

                }
            }
			sizeChange_=0;
		}
	}

    typedef BaseT* Iterator;
    typedef const Iterator const_Iterator;


  
  
    std::ostream& write(std::ostream& stream) const
{
	if(stream.good())
	{
		stream << capacity_ << " ";     // in objects(indexes); end size of pool
		stream << size_ << " ";     // in obj; last used from pool,last existing (used) element
		stream << uniqueID_ << " ";  //
		stream << sizeChange_ << " " ;
		const std::streamsize size(size_*sizeof(BaseT));
		stream.write(reinterpret_cast<const char*>(pool_),size);
	}
	assert(stream.good());
	return stream;
}
  
  std::istream& read(std::istream& stream)
{
	if(stream.good())
	  {
		stream >> capacity_;     // in objects(indexes); end size of pool
		assert(stream.good());
		//std::cout << "\ncapacity=" << capacity_  ;
		reserve(capacity_,true);
		stream >> size_;     // in obj; last used from pool,last existing (used) element
		assert(stream.good());
		//std::cout << "\nsize=" << size_ ;
		stream >> uniqueID_;  //
		assert(stream.good());
		//std::cout << "\nuniqueId=" << uniqueID_ ;
		stream >> sizeChange_;
		assert(stream.good());
		//std::cout << "\nsizeChange=" << sizeChange_;
		char tmp('E');
		stream.get(tmp);
		assert(tmp == ' ');
		const std::streamsize size(size_*sizeof(BaseT));
		stream.read(reinterpret_cast<char*>(pool_),size);
		assert(stream.gcount() == static_cast<std::streamsize>(size_*sizeof(BaseT)));
	  }
	assert(stream.good());
	return stream;
}

  
private:
	uTind    uniqueID()
	{
		return static_cast<uTind>(++uniqueID_);
	}

    BYTE * pool_;
    size_t size_,uniqueID_,sizeChange_,capacity_,nDeleted_;
  //hashPosByVts vts2pos_;
};

/**  @} */

#endif // ARRAYPOOL_HPP_INCLUDED
