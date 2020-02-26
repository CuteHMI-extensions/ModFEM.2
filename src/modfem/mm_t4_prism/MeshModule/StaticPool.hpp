#ifndef STATICPOOL_HPP_INCLUDED
#define STATICPOOL_HPP_INCLUDED

#include <vector>

#include <cassert>
#include <string.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stack>




/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


#include "../Common.h"
#include "VtsSqId.hpp"

class hHybridMesh;

template< class BaseT, int TMaxVerts = BaseT::nVerts>
class   StaticPool
{
public:
  // up to this moment we have max 6 vertices in entity

  typedef VtsSqId<TMaxVerts> VtsId;
  typedef unordered_map< VtsId, uTind> HashPosByVts;
  typedef typename HashPosByVts::value_type VtsPos;
  typedef typename HashPosByVts::iterator Vts2posIter;
  typedef typename HashPosByVts::const_iterator Vts2posCIter;
 
  typedef BaseT* BaseTPtr;
private:
	static const int MaxVerts = TMaxVerts;
    uTind hashSize_;
	uTind capacity_;     // in objects(indexes); end size of pool
	uTind size_;     // in obj; last used from pool,last existing (used) element
	uTlong capacityB_;   // in bytes; real end size of pool
	BYTE* poolB_;     // pool (bytes)
	BYTE* endB_;      // in bytes; pointer to next after size_ element in pool
	BaseTPtr *hash_;   //hash_ table for random access
	uTind uniqueID_;  //
  Tlong memoryChangeB_;
	Tind  capacityChange_;
    std::stack<uTind>    emptyPos_;

    hHybridMesh * m;

  uTind    uniqueID()
	{
        uTind new_id = ++uniqueID_;
        updateHashSize() ;
        return new_id;
	}

public: 
  uTind dividedObjs_;

  HashPosByVts vts2pos_; 

    StaticPool(hHybridMesh * mesh): m(mesh), hashSize_(0),capacity_(0),size_(0),capacityB_(0),poolB_(NULL),endB_(NULL),hash_(NULL),uniqueID_(0),memoryChangeB_(0),capacityChange_(0), dividedObjs_(0)
	{
		clear();
	}
  
  ~StaticPool()
  {
	clear();
  }
  
  //	void initialize(const uTind newCapacity,const uTind newCapacityB = 0)
  //	{
  //	reserve(newCapacity,newCapacityB);
  //}

	//const uTind end()   const { return size_+1;}   // first after last element
	uTind last()  const {
        return uniqueID_;
	}
	BaseT* end() const {
        return reinterpret_cast<BaseTPtr>(endB_);
	}
	uTind size()  const {
		return size_;     // size in obj
	}
	uTind capacity() const {
		return capacity_;	// max size in obj
	}

	bool empty()  const {
		return (size_ == 0);    // true if pool is empty
	}

	void clear()
	{
	  vts2pos_.clear();
      while( emptyPos_.size() > 0) {
          emptyPos_.pop();

      }

	  uniqueID_ = 0;
		capacity_ = 0;
        updateHashSize();
		size_ = 0;
		safeDeleteArray(hash_);
        hash_ = NULL;

		capacityB_ = 0;
        endB_ = NULL;
		safeDeleteArray(poolB_);
        poolB_ = NULL;

		dividedObjs_=0;
	}

	void reserve(const uTind newCapacity,uTind newCapacityB=0, bool forceResize = false)
	{
	  //assert(newCapacity >= 0);
 
		if (newCapacity > 0) {
			capacity_ = newCapacity;
            updateHashSize();
			if (newCapacityB == 0) 	{
				newCapacityB = sizeof(BaseT)*newCapacity;
			}
			//assert(newCapacityB >= 0);
			if (newCapacityB > capacityB_ || forceResize) {
				const uTind oldSizeB(static_cast<uTind>(endB_ - poolB_));
				assert(oldSizeB <= capacityB_);
				const uTind oldCapacityB(capacityB_);
                BYTE * _newPool(NULL);
				try
				{
					capacityB_ = newCapacityB;
					_newPool = new BYTE[capacityB_+sizeof(BaseT)];
					memset(_newPool,0,capacityB_+sizeof(BaseT));
				}
				catch (const std::bad_alloc & e)
				{
					std::string msg(e.what());
					msg.append(" StaticPool::resize( ");
					char	tmp[64]={};
					//itoa(newCapacity,tmp,64,10);
					msg.append(tmp);
					msg.append(", ");
					//itoa(newCapacityB,tmp,64,10);
					msg.append(tmp);
					msg.append(") not enough memory to resize pool.");
					std::bad_alloc	err;
					throw err;
				}

                if (poolB_ != NULL) {
					if (forceResize && newCapacityB < oldCapacityB) {
						memcpy(_newPool,poolB_,newCapacityB);
						endB_ = _newPool+newCapacityB;
					}
					else {
						memcpy(_newPool,poolB_,oldSizeB);
						endB_ = _newPool+oldSizeB;
					}
					safeDeleteArray(poolB_);
				}
				else { // pool was empty
					endB_=_newPool;
				}
				poolB_ = _newPool;
                _newPool = NULL;
				rebuildHash();
			}
		}
	}

	void    reallocate(const uTind newCapacity,uTlong newCapacityB=0, bool forceResize = false)
	{
	  //assert(newCapacity >= 0);
	  //assert(newCapacityB >= 0);
	  
		if (newCapacity > 0)
		{
		  //const uTind oldCapacityB(capacityB_),oldCapacity(capacity_);
			capacity_ = newCapacity;
            updateHashSize();
			if (newCapacityB == 0)
			{
				newCapacityB = sizeof(BaseT)*newCapacity;
			}
			
			assert(newCapacityB / capacity_ >= sizeof(BaseT)); // it should be at least space for BaseT
			//std::cout << "Resizing StaticPool from " << " to " <<;
			
			if (newCapacityB > capacityB_ || forceResize)
			{
				//const uTind oldCapacityB(capacityB_);
                BYTE * _newPool(NULL);
				try
				{
					capacityB_ = newCapacityB;
					_newPool = new BYTE[capacityB_+sizeof(BaseT)];
					memset(_newPool,0,capacityB_+sizeof(BaseT));

					//delete [] hash_;
					//hash_ = new BaseT*[capacity_+FIRST];
					//memset(hash_,0, sizeof(BaseT*) * (capacity_+FIRST));
					rebuildHash();
				}
				catch (const std::bad_alloc & e)
				{
					std::string msg(e.what());
					msg.append(" StaticPool::reallocate(newCapacity=");
					char	tmp[64]={};
					//itoa(newCapacity,tmp,64,10);
					msg.append(tmp);
					msg.append(", newCapacityB=");
					//itoa(newCapacityB,tmp,64,10);
					msg.append(tmp);
					msg.append(") not enough memory to resize pool.");
					std::bad_alloc	err;
					throw err;
				}

                if (poolB_ != NULL)
				{

					size_t  offset(0);
					BaseTPtr	oldEndIt(reinterpret_cast<BaseT*>(endB_));
					BYTE*	_oldPool(poolB_);
					poolB_ = _newPool;	//!!
					endB_ = poolB_;
					for(BaseTPtr  it(reinterpret_cast<BaseTPtr>(_oldPool)); it < oldEndIt;
						it = reinterpret_cast<BaseTPtr>(it->next()) )
					  {
						BaseTPtr moved_obj(reinterpret_cast<BaseTPtr>(endB_));
                        assert(moved_obj != NULL);
                        if(it->nMyClassSons_ != BaseT::delMark) {
						  
						  memcpy(moved_obj,it,it->mySize_); // actual SHALLOW copying
						  hash_[moved_obj->pos_]=reinterpret_cast<BaseTPtr>(moved_obj); // hash_ copied elem
						  offset+=moved_obj->mySize_; // move to next position
						  //std::cout << "\nMoving " << moved_obj->pos_;
						  switch (moved_obj->nMyClassSons_)
						  {
							case BaseT::fullRefMark:  // if elem is marked for refine
							case BaseT::partialRefMark:
                              offset += moved_obj->hBreak(m); // left space for child elements
							  //std::cout << " Breaking! ";
							  break;
							case BaseT::derefMark: // back to state before adaptation
								assert(!"Deref mark found while reallocationg!");
                                moved_obj->derefine(m);
							  break;
				  		  }// otherwise do nothing
                          endB_ = poolB_+offset;
                          moved_obj->updatePointers(); // !!!!!! MUST DO IT !!!!!
                          assert(it->equals(*moved_obj));
						}
                        else { // deleting
                            --size_;
                            hash_[it->pos_] = NULL; // removing entry from hash table
                            emptyPos_.push(it->pos_);
                            //vts2pos_.erase(VtsId(it->verts()));
                            //assert(!"Delete mark found while reallocationg!");
						}

					}
                    rebuildVts2Pos();
                    //assert(first().myMesh->test());
                    //

					//if (forceResize && newCapacityB < oldCapacityB)
					//{
					//   //memcpy(_newPool,poolB_,newCapacityB);
					//    endB_ = _newPool+newCapacityB;
					//}
					//else
					//{
					//   //memcpy(_newPool,poolB_,oldCapacityB);
					//    endB_ = _newPool+(endB_-poolB_);
					//}

					safeDeleteArray(_oldPool);
				}
				else // pool was empty
				{
					endB_=poolB_;
				}
			}
		}
        else {
            this->clear();
        }
		checkHash();
        if(size_ != 0) {
            assert(capacityB_ / size_ >= sizeof(BaseT));
        }
        if(capacity_ != 0){
            assert(capacityB_ / capacity_ >= sizeof(BaseT));
        }
        assert(endB_ >= poolB_);
		assert(endB_ <= poolB_+newCapacityB);
		assert(static_cast<uTind>(endB_-poolB_) <= capacityB_);
		assert(size_ <= capacity_);
	}
  
  void rebuildVts2Pos()
  {
#ifndef _MSC_VER
      if(MaxVerts > 0)
#endif
      {
        vts2pos_.clear();
        for(constIterator<> it(this); !it.done(); ++it) {
          uTind verts[VtsId::length];
          uTind i(0);
          for(;i < it->typeSpecyfic_.nVerts_;++i) {
            verts[i]=it->verts(i);
          }
          for(;i < static_cast<uTind>(VtsId::length); ++i) {
            verts[i]=UNKNOWN;
          }
          assert(i == static_cast<uTind>(VtsId::length));
          assert(vts2pos_.find(verts) == vts2pos_.end());

          vts2pos_.insert(VtsPos(verts,it->pos_));

          VtsId tmpId(verts);
          assert(vts2pos_.count(VtsId(verts))==1);
          Vts2posIter iter = vts2pos_.find(VtsId(verts));
          assert( iter != vts2pos_.end());
          assert( tmpId == iter->first );
          assert(iter->second  == it->pos_);
        }
      }
  }
  
  inline uTind hashSize() const { return hashSize_; }
  void updateHashSize()  { hashSize_ = (capacity_==0) ? 0 : FIRST + std::max(capacity_,uniqueID_); }

	void rebuildHash()
	{
		safeDeleteArray(hash_);
		if (capacity_ > 0)
		{
			try
			{
                this->hash_ = new BaseTPtr[hashSize()];
                memset(hash_, 0, sizeof(BaseTPtr) * (hashSize()));
			}
			catch (const std::bad_alloc &)
			{
				throw "StaticPool::rebuildHash: not enugh memory to allocate";
			}
			//uTind i(FIRST);

            assert(hash_[0] == NULL);
			for (Iterator<>  it(this);!it.done();++it)
			{
			    assert(it->pos_ >= FIRST);
                assert(it->pos_ <= hashSize());
                assert(hash_[it->pos_]== NULL);

				hash_[it->pos_]= &it;
				it->updatePointers();

				assert(hash_[it->pos_]->pos_ == it->pos_);
			}
		}
	checkHash();
	}

	void checkHash() const
	{

        for(int i=0; i < hashSize(); ++i)
        {
            assert( (hash_[i] == NULL)
                    || ( (reinterpret_cast<BYTE*>(hash_[i]) >= poolB_)
                      && (reinterpret_cast<BYTE*>(hash_[i]) < endB_) ) );
        }

	  for(constIterator<allObj> it(this); !it.done(); ++it)
	  {
		assert(getById(it->id_).id_== it->id_);
		//it->test();
	  }
	}

	BaseT & first()
    {
	  assert(!empty());
		return *reinterpret_cast<BaseT*>(poolB_);
	}

	const BaseT & first() const
    {
	  assert(!empty());
	  return *reinterpret_cast<const BaseT*>(poolB_);
    }

    inline BaseT & at(const uTind pos_)
	{
        assert(pos_<= hashSize());
		assert(pos_ >= FIRST);
        assert(hash_[pos_] != NULL);
		assert(hash_[pos_]->pos_== pos_);
		return *hash_[pos_];
	}

	template<class RetType>
    inline RetType & at(const uTind pos_)
	{
		assert(RetType::myType == hash_[pos_]->type_);
        assert(pos_<= hashSize());
		assert(pos_ >= FIRST);
        assert(hash_[pos_] != NULL);
		assert(hash_[pos_]->pos_== pos_);
		return *reinterpret_cast<RetType*>(hash_[pos_]);
	}

    inline const BaseT & at(const uTind pos_) const
	{
        assert(pos_<= hashSize());
		assert(hash_[pos_]->pos_== pos_);
		return *hash_[pos_];
	}

	template<class RetType>
    inline const RetType & at(const uTind pos_) const
	{
		assert(RetType::myType == hash_[pos_]->type_);
        assert(pos_<= hashSize());
		assert(hash_[pos_]->pos_== pos_);
		return *static_cast<RetType*>(hash_[pos_]);
	}

    inline BaseT & getById(const ID id)
  {
	    assert( at(BaseT::posFromId(id)).type_== BaseT::typeFromId(id));
		assert( at(BaseT::posFromId(id)).id_== id);
		return at(BaseT::posFromId(id));
	}

	template<class RetType>
    inline RetType & getById(const ID id)
	{

		assert( at(BaseT::posFromId(id)).type_== BaseT::typeFromId(id));
		assert( at(BaseT::posFromId(id)).id_== id);
		return at<RetType>(BaseT::posFromId(id));
	}

    inline const BaseT & getById(const ID id) const
	{
		assert( at(BaseT::posFromId(id)).type_== BaseT::typeFromId(id));
		assert( at(BaseT::posFromId(id)).id_== id);
		return at(BaseT::posFromId(id));
	}

	template<class RetType>
    inline const RetType & getById(const ID id) const
	{
		assert( at(BaseT::posFromId(id)).type_== BaseT::typeFromId(id));
		assert( at(BaseT::posFromId(id)).id_== id);
		return at<RetType>(BaseT::posFromId(id));
	}

	// operator to get by pos
    inline const BaseT& operator[](const Tind pos_) const
    {
        assert(pos_<= hashSize());
        assert(pos_ >= FIRST);
        return *hash_[pos_];
    }
    inline BaseT& operator[](const Tind pos_)
    {
        assert(pos_<= hashSize());
        assert(pos_ >= FIRST);
        return *hash_[pos_];
    }

	// operator to get by id
    inline const BaseT& operator()(const ID id) const { return getById(id); }
    inline BaseT& operator()(const ID id)	   { return getById(id); }
  
private:

	void checkIfInVts2Pos(const uTind verts[]) const 
	{
#ifndef _MSC_VER
        if(MaxVerts > 0)
#endif
        {
            if(vts2pos_.find(VtsId(verts)) !=  vts2pos_.end() )  {
                mf_log_err("Duplicate object detected (max verts %d)!",MaxVerts);
                mf_print_array(verts,MaxVerts,"%d");
            }

        }
	}

	void insertIntoVts2Pos(const uTind verts[],const int pos) 
	{
#ifndef _MSC_VER
        if(MaxVerts > 0)
#endif
        {
            //assert(hash_[size_]->myMesh != NULL);
            //hash_[size_]->myMesh->registerHObj(hash_[size_]->type_,verts,hash_[size_]->pos_);
            // replaced by:
            vts2pos_.insert( VtsPos(verts, static_cast<uTind>(pos) ) );
            assert(vts2pos_[verts] == pos);
        }
	}

public:



	template< class T>
    T*  newObj(hHybridMesh * myMesh,const uTind verts[],void * memPtr=NULL)
  {
    uTind objpos = -1;
    if(!emptyPos_.empty()) {
        objpos = emptyPos_.top();
        emptyPos_.pop();
    }
    else {
        objpos = (uniqueID());
    }

    void * ptr(memPtr==NULL ? endB_ : memPtr);
        assert( objpos < hashSize());
		assert(usedMemory() < totalMemory());
        assert(hash_[objpos]==NULL);
		assert(ptr >= poolB_);
		assert(ptr < poolB_+capacityB_);
		// check if entity is NOT already in pool...
		checkIfInVts2Pos(verts);
	
        hash_[objpos] =  reinterpret_cast<BaseT*>( new(ptr) T(myMesh,verts,objpos) );
        ++size_;
        if(memPtr == NULL) {
		  endB_+=sizeof(T);
		}

		insertIntoVts2Pos(verts,hash_[objpos]->pos_);
		
        assert(uniqueID_ >= objpos);	// checking if all objcts are created by newObj(..)
		assert(usedMemory() <= totalMemory());
        assert(hash_[objpos] != NULL);
        assert( size_ <= FIRST+capacity_);
        return  reinterpret_cast<T*>(hash_[objpos]);
	}

	uTlong totalMemory() const
	{
		return capacityB_;
	}

	uTlong usedMemory() const
	{
		assert(endB_ >= poolB_);
		return static_cast<uTlong>(endB_-poolB_);
	}

    void    requestChange(const Tind sizeChange,const Tlong memChange)
	{
		capacityChange_+=sizeChange;
		memoryChangeB_+=memChange;
		assert(static_cast<size_t>(memoryChangeB_ / capacityChange_) >= sizeof(BaseT));
	}

	uTlong memoryNeeded() const  // in bytes
	{
		return usedMemory() + memoryChangeB_;
	}

	void    adjust()
	{
		/// Pre: +all markings have beed done (correctly).
		///     + pool is not empty
        if (memoryChangeB_!=0 || capacityChange_!= 0)
		  {
			assert(memoryChangeB_ / capacityChange_ >= static_cast<Tlong>(sizeof(BaseT)));
			reallocate(size()+capacityChange_, memoryNeeded(),true);
			memoryChangeB_=0;
			capacityChange_=0;
		}
	}


	uTind	nonDividedObjs() const
	{
		return size() - dividedObjs_;
	}


  //void register(const uTind vertices[],const uTind pos)
  //{
  //	
  //}
  
  //// I/O operators
  
	std::ofstream & operator << (std::ofstream & stream)
	{
	  return this->write(stream);
	}

	std::ifstream & operator >> (std::ifstream & stream)
	{
	  return this->read(stream);
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	struct	allObj
	{
	  static bool	check(const BaseT*) {return false;}
	};
	struct	notBroken
	{
	  static bool	check(const BaseT* obj) { return obj->isBroken(); }
	};
	struct atBoundary
	{
	  static bool check(const BaseT* obj) { return obj->isAtBoundary(); }
	};
	
	template< class TCondition = allObj>
	class   Iterator
	{
	public:

		Iterator(StaticPool * container) : _container(container)
		{
			if (!_container->empty())
			{
				_obj = reinterpret_cast<BaseT*>(_container->poolB_);
				while(TCondition::check(_obj))
				{
				  this->operator++();
				}
			}
			else
			{
				_obj = reinterpret_cast<BaseT*>(_container->endB_);
			}
		}

		Iterator(StaticPool * container,const uTind pos_) : _container(container)
		{
			if (!container->empty() && pos_<= _container->last())
			{
				_obj = _container->at(pos_);
			}
			else
			{
				_obj = reinterpret_cast<BaseT*>(_container->endB_);
			}
		}

		template< class T>
		T*  as() {
			return reinterpret_cast<T*>(_obj);
		}

		BaseT*  operator&()  {
			return _obj;
		}
		BaseT&  operator*() {
			return *_obj;
		}
		BaseT*  operator->() {
			return _obj;
		}
		BaseT*  operator++()
		{
			assert(_obj->mySize_> 0);
			do{
				_obj = reinterpret_cast<BaseT*>( reinterpret_cast<BYTE*>(_obj) + _obj->mySize_ );
			}while(TCondition::check(_obj));
			assert(reinterpret_cast<BYTE*>(_obj) <= _container->endB_);
			return _obj;
		}
		BaseT*  operator++(int)
		{
			const BaseT* old(_obj);
			operator++();
			return old;
		}

		bool    done() const {
			return (reinterpret_cast<BYTE*>(_obj) >= _container->endB_);
		}
		bool    operator==(const Iterator & other) const  {
			return _container==other._container && _obj==other._obj;
		}
		bool    operator!=(const Iterator & other) const  {
			return !(*this==other);
		}

	private:
		StaticPool * _container;
		BaseT      *_obj;
	};

	template< class TCondition = allObj>
	class   constIterator
	{
	public:

		constIterator(const StaticPool * container) : _container(container)
		{
			if (!_container->empty())
			{
				_obj = reinterpret_cast<const BaseT*>(_container->poolB_);
				while(TCondition::check(_obj))
				{
				  this->operator++();
				}
			}
			else
			{
				_obj = reinterpret_cast<const BaseT*>(_container->endB_);
			}
		}

		constIterator(const StaticPool * container,const uTind pos_) : _container(container)
		{
			if (!container->empty() && pos_<= _container->last())
			{
				_obj = &_container->at(pos_);
			}
			else
			{
				_obj = reinterpret_cast<const BaseT*>(_container->endB_);
			}
		}

		template< class T>
		const T*  as() {
			return reinterpret_cast<const T*>(_obj);
		}

		const BaseT*  operator&() const {
			return _obj;
		}
		const BaseT&  operator*() const  {
			return *_obj;
		}
		const BaseT*  operator->()const  {
			return _obj;
		}
		const BaseT*  operator++()
		{
			assert(_obj->mySize_> 0);
			do{
				_obj = reinterpret_cast<const BaseT*>( reinterpret_cast<const BYTE*>(_obj) + _obj->mySize_ );
			}while(TCondition::check(_obj));
			assert( reinterpret_cast<const BYTE*>(_obj) <= _container->endB_);
			return _obj;
		}
		const BaseT*  operator++(int)
		{
			const BaseT* old(reinterpret_cast<const BaseT*>(_obj));
			operator++();
			return old;
		}

		bool    done() const {
			return (reinterpret_cast<const BYTE*>(_obj) >= _container->endB_);
		}
		bool    operator==(const constIterator & other) const  {
			return _container==other._container && _obj==other._obj;
		}
		bool    operator!=(const constIterator & other) const  {
			return !(*this==other);
		}

	private:
		const StaticPool * _container;
		const BaseT      * _obj;
	};

	constIterator<notBroken>    begin() const {
		return constIterator<notBroken>(this);
	}
	//constIterator    end()   const { return constIterator(this,size_);}

	Iterator<notBroken>    begin() {
		return Iterator<notBroken>(this);
	}
	//Iterator    end()   { return Iterator(this,size_);}

  
  std::ostream& write(std::ostream& stream) const
{
	if(stream.good())
	{
	  stream << capacity_ << " "; assert(stream.good());    // in objects(indexes); end size of pool
	  stream << size_ << " "     // in obj; last used from pool,last existing (used) element
			 << capacityB_ << " "   // in bytes; real end size of pool
			 << (endB_-poolB_) << " " // actual used memory in Bytes
			 << uniqueID_ << " "  //
			 << memoryChangeB_ << " "
			 << capacityChange_ << " "; assert(stream.good());
	  stream.write(reinterpret_cast<const char*>(poolB_),(endB_-poolB_));
		//BYTE* poolB_;     // pool (bytes)
		//BYTE* endB_;      // in bytes; pointer to next after size_ element in pool
		//BaseT* *hash_;   //hash_ table for random access
	}
	return stream;
}
  
  std::istream& read(std::istream& stream)
  {
	if(stream.good())
	  {
		int usedMemB(0);
	  stream >> capacity_;   assert(stream.good());// in objects(indexes); end size of pool
	  stream >> size_ ;    assert(stream.good());// in obj; last used from pool,last existing (used) element
	  stream >> capacityB_;   assert(stream.good());// in bytes; real end size of pool
	  stream >> usedMemB; assert(stream.good());
	  stream >> uniqueID_;  assert(stream.good());//
	  stream >> memoryChangeB_; assert(stream.good());
	  stream >> capacityChange_; assert(stream.good());
	  assert(stream.good());
	  assert(static_cast<uTind>(usedMemB) <= capacityB_);
	  std::cout << " capacity:" << capacity_
				<< " size:" << size_
				<< " usedMemB:" << usedMemB
				<< " uniqueId:"<< uniqueID_
				<< " memoryChangeB:" << memoryChangeB_
				<< " capacitChange:" << capacityChange_;

	  reserve(capacity_,capacityB_,true);
	  char tmp('E');
	  stream.get(tmp);
	  stream.read(reinterpret_cast<char *>(poolB_),usedMemB);
	  endB_ +=usedMemB;
	  rebuildHash();
	  rebuildVts2Pos();
	  //BYTE* poolB_;     // pool (bytes)
		//BYTE* endB_;      // in bytes; pointer to next after size_ element in pool
		//BaseT* *hash_;   //hash_ table for random access

	}
	return stream;
}
  
  
template<class T,int T1> friend std::ostream& operator<<(std::ostream& os,const StaticPool<T,T1> & pool);
template<class T,int T1> friend std::istream& operator>>(std::istream& os,const StaticPool<T,T1> & pool);
};

template <class T,int T1>
std::ostream & operator << (std::ostream & stream, const StaticPool<T, T1> & pool)
{
  return pool.write(stream);
}

template <class T, int T1>
std::istream & operator >> (std::istream & stream, StaticPool<T, T1> & pool)
{
  return pool.read(stream);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**  @} */
#endif // STATICPOOL_HPP_INCLUDED
