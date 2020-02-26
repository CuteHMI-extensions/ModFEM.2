#ifndef _FIELD_HPP_
#define _FIELD_HPP_

#include <cassert>
#include <string.h>
#include <iostream>

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */

namespace Memory {

	template <int DIM, typename T>
	class Field;

	template <int DIM, typename T>
	std::ostream & operator<< (std::ostream & os, const Field<DIM,T> & f);

	//////////////////////////////////////////////////////////////////////////
	template <int DIM, typename T>
	class Field {
	public:
		static const int size = DIM;
				 Field()				{ init(0); }
		explicit Field(const int c)		{ init(c); }
				 Field(const Field& f)	{ memcpy(this, &f, DIM * sizeof(T)); }

		void init(const int c)	{ memset(this,  c, DIM * sizeof(T)); }

        inline       T& operator[](const int i)			{ assert(i >= 0 && i < DIM); return arr[i]; }
        inline const T& operator[](const int i)	const	{ assert(i >= 0 && i < DIM); return arr[i]; }

		const Field& operator=(const T val) {
			for(register int i = 0; i < DIM; ++i)
				arr[i] = val;
			return *this;
		}
		const Field& operator=(const Field & f)	{
			if(this != &f)
				memcpy(this, &f, DIM * sizeof(T));
			return *this;
		}

		operator const T * () const { return arr; } //conversion operator
		operator T*		   ()		{ return arr; }

		void	addAX(register const Field& f, register const T a = 1.0) {
			for(register int i = 0; i < DIM; ++i)
				arr[i] += f.arr[i] * a;
		}

		T*		getArray() { return arr; }

		friend std::ostream& operator<< <>(std::ostream& os, const Field& f);

	private:
		T	arr[DIM];
	};

	//////////////////////////////////////////////////////////////////////////
	template <typename T>
	class Field<1, T>
	{
	public:
		Field()					: v(T(0))	{}
		Field(const T & c)		: v(c)		{}
		Field(const Field& f)	: v(f.v)	{}

		void	init(const int c)	{ v = c; }

		      T& operator[](const int)			{ return v; }
		const T& operator[](const int)	const	{ return v; }

		const Field& operator=(const T val)		{ v = val; return *this; }
		const Field& operator=(const Field& f)	{ v = f.v; return *this; }

		void	addAX(register const Field& f, register const T a = 1.0) { v += a * f.v; }

		friend std::ostream& operator<< <>(std::ostream& os, const Field& f);

	private:
		T	v;
	};

	//////////////////////////////////////////////////////////////////////////
	
	template <typename T>
	class Field<0, T>
	{
	public:
		void	init(const int) {};
	};

	template <int DIM>
	class Field<DIM, void>
	{
	public:
		void	init(const int) {};
	};

	//////////////////////////////////////////////////////////////////////////
	
	template <int DIM, typename T>
	std::ostream& operator<<(std::ostream& os, const Field<DIM, T>& f) {
		os << "(";
		for(register int i = 0; i < DIM; ++i)
			os << f.arr[i] << ", ";
		os << "\b\b)";
		return os;
	}

	template <typename T>
	std::ostream& operator<<(std::ostream& os, const Field<1, T>& f) {
		os << "(" << f.v << ")";
		return os;
	}

	//////////////////////////////////////////////////////////////////////////
	// coordinates
	//typedef Field<1, Tval>	Coord1D;
	//typedef Field<2, Tval>	Coord2D;
	//typedef Field<3, Tval>	Coord3D;

	//////////////////////////////////////////////////////////////////////////
	// degrees of freedom
	//typedef	Field<1, Tval>	Dofs1D;
	//typedef	Field<2, Tval>	Dofs2D;
	//typedef	Field<3, Tval>	Dofs3D;
	//typedef	Field<4, Tval>	Dofs4D;

	template <typename T>
	struct Type2Ref {
		typedef  T& refType;
		typedef const T& constRefType;
	};

	template <>
	struct Type2Ref<void> {
		typedef  void refType;
		typedef  void constRefType;
	};
}


/** @}*/

#endif
