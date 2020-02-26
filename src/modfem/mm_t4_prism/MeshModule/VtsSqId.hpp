#ifndef VTS_SQ_ID_HPP
#define VTS_SQ_ID_HPP

#include "../Common.h"

#include <algorithm>
#include <ostream>

#include "SuperFastHash.h"
#include <boost/config.hpp> //< for BOOST_NO_CXX11_HDR_UNORDERED_MAP

/// Vertices Sequence Identyfier
/// assure that vertices are correcly ordered (sorted)
/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


#ifdef _MSC_VER

template<int TLength>
class VtsSqId;

template<>
class VtsSqId<0> {};

#endif

template<int TLength>
class VtsSqId
{
  uTind v[TLength];
public:
  static const int length=TLength;
  VtsSqId(const uTind vts[])
  {
	memcpy(v,vts,sizeof(uTind)*TLength);
	std::sort(v,v+TLength);
  }
  
  VtsSqId(const VtsSqId & other)
  {
	memcpy(v,other.v,sizeof(uTind)*TLength);
  }

  /// OPERATORS  
  VtsSqId& operator=(const VtsSqId & other)
  {
	memcpy(v,other.v, sizeof(uTind)*TLength);
	return *this;
  }

  inline bool operator==(const VtsSqId & other) const
  {
	bool eq(true);
	for(int i(0); (i < TLength) && eq; ++i) {
	  if(v[i] != other.v[i]){
		eq=false;
	  }
	}
	return eq;
  }

 inline  bool operator!=(const VtsSqId & other) const
  {
	return !(*this==other);
  }

  inline bool operator<(const VtsSqId & other) const
  {
	register int i(0);
	while((v[i]==other.v[i]) && (i < TLength-1) ) {
	  ++i;
	}
	return v[i] < other.v[i];
  }

  inline bool operator>=(const VtsSqId & other) const
  {
	return ! ( this->operator<(other) );
  }

  inline bool operator>(const VtsSqId & other) const
  {
	register int i(0);
	while((v[i]==other.v[i]) && (i < TLength-1) ) {
	  ++i;
	}
	return v[i] > other.v[i];
  }
  
  inline bool operator<=(const VtsSqId & other) const
  {
	return ! ( this->operator>(other) );
  }

  
  
  template<int T>
  friend std::ostream &  operator<<(std::ostream & os,const VtsSqId<T> & id);
  
};
template<int T>
std::ostream &  operator<<(std::ostream & os,const VtsSqId<T> & id)
{
  for(int i(0); i < VtsSqId<T>::length; ++i) {
	os << id.v[i] << " ";
  }
  return os;
}

/**  @} */



#ifdef BOOST_NO_CXX11_HDR_UNORDERED_MAP
#include <boost/unordered/unordered_map.hpp>

namespace boost {

template<int T>
struct hash<VtsSqId <T> >
{
    typedef VtsSqId <T> argument_type;
    typedef size_t result_type;

    result_type operator () (const argument_type& x) const {
        return SuperFastHash( reinterpret_cast<const char*>(&x),sizeof(x));
    }
};

}

using boost::unordered::unordered_map;
#else
#include <unordered_map>
/// Hashes for std::unordered_map
namespace std {

template<int T>
struct hash<VtsSqId <T>>
{
    typedef VtsSqId <T> argument_type;
    typedef size_t result_type;

    result_type operator () (const argument_type& x) const {
        return SuperFastHash( reinterpret_cast<const char*>(&x),sizeof(x));
    }
};

template<int T>
    struct __is_fast_hash< hash< VtsSqId <T> > > : public std::false_type
    { };

}

using std::unordered_map;

#endif

#endif //guard header
