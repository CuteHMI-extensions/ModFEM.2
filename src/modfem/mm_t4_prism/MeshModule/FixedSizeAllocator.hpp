// ***************************************************************
//  FixedSizeAllocator   version:  1.0   date: 05/26/2008
//  -------------------------------------------------------------
//  Kazimierz Michalik
//  -------------------------------------------------------------
//  Copyright (C) 2008 - All Rights Reserved
// ***************************************************************
// 
// ***************************************************************
#ifndef ChunkList2_h__
#define ChunkList2_h__

#include <vector>

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


namespace Memory {
	using std::size_t;

	template<typename TField>
	class FixedSizeAllocator {
	protected:
		// dependent types
		struct Chunk;
		typedef typename std::vector<Chunk*>	vChunk;
		typedef typename vChunk::iterator		iChunk;
		typedef typename Chunk::BLOCK*			pBlock;
        typedef const typename Chunk::BLOCK*	cBlock;

	public:
		FixedSizeAllocator(const size_t size, const size_t num)
			: _blockSize(size), _blockNum(num), _firstFree(NULL) {
			_chunks.reserve(64);
		}

		~FixedSizeAllocator() {	clear(); }

		void	init(const size_t size, const size_t num) {
			clear();
			_blockSize	= size;
			_blockNum	= num;
			_firstFree	= 0;
		}

		TField*	alloc() {
			if(_firstFree == NULL)	// chunk is full or it's empty
				addChunk();
			_ptr = reinterpret_cast<TField*>(_firstFree);
			_firstFree = _firstFree->_next;
			return new (_ptr) TField();
		}

		void	free(register TField *const ptr) {
			if(ptr) {
				reinterpret_cast<pBlock>(ptr)->_next = _firstFree;
				_firstFree = reinterpret_cast<pBlock>(ptr);
			}
		}

	protected:
		FixedSizeAllocator(const FixedSizeAllocator &);

		void	clear() {
			for(register iChunk it(_chunks.begin()); it != _chunks.end(); ++it) 
				delete *it;
			_chunks.clear();
		}

		void	addChunk() {
			register Chunk *const _lastChunk = new Chunk(_blockSize, _blockNum);
			_chunks.push_back(_lastChunk);
			_firstFree = _lastChunk->_storage;
		}

		size_t		_blockSize;
		size_t		_blockNum;
		pBlock		_firstFree;
		vChunk		_chunks;
		TField*		_ptr;	
	};

	//////////////////////////////////////////////////////////////////////////
	template<typename TField>
	struct FixedSizeAllocator<TField>::Chunk {
		union BLOCK {
			unsigned char	_val[sizeof(TField)];
			BLOCK*			_next;
		} *_storage;
        size_t	_blockSize;
		size_t	_blockNum;

		Chunk(const size_t size, const size_t num)
			: _storage(new BLOCK[size * num]), _blockSize(size), _blockNum(num) {
			assert(sizeof(BLOCK) == sizeof(TField));
			memset(_storage, 0, size * num * sizeof(BLOCK));
			reset();
		}

		~Chunk() { delete [] _storage; }

		void reset() {
			register BLOCK* ptr(_storage);
			for(register BLOCK *const end(_storage + _blockSize * (_blockNum - 1)); ptr < end; ptr = ptr->_next)
				ptr->_next = ptr + _blockSize;
			ptr->_next = NULL;
		}

		bool	inside(register const TField *const ptr) const {
			register const TField *const begStorage = reinterpret_cast<TField*>(_storage);
			return ptr >= begStorage && ptr < begStorage + _blockNum;
		}
	};
}
/**  @} */
#endif // ChunkList2_h__
