/** @file Array_view.hpp
 
	Copyright (c) 2016-2017 Santeri Puranen.

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Affero General Public License as
	published by the Free Software Foundation, either version 3 of the
	License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Affero General Public License for more details.

	You should have received a copy of the GNU Affero General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.

	@author Santeri Puranen
	$Id: $
*/
#ifndef APEGRUNT_ARRAY_VIEW_HPP
#define APEGRUNT_ARRAY_VIEW_HPP

#include <iterator>
#include <type_traits>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
//#include <boost/range/iterator_range_core.hpp>

#include <boost/accumulators/accumulators.hpp>

#include "misc/Array_view_forward.h"
#include "misc/Array_view_iterator_forward.h"

namespace apegrunt {

template< typename T, std::size_t Extent, std::size_t Rank >
class Array_view
{
public:
	using value_type = T;
	using real_t = typename scalar_type<value_type>::type;

	using my_type = Array_view<T,Extent>;

	using iterator = Array_view_iterator<value_type>;
	using const_iterator = const Array_view_iterator<value_type>;

	Array_view( real_t * const data, std::size_t extent=Extent ) : m_data(data), m_extent(extent), m_nT(std::max(extent,recursive_extent<my_type>::value)*recursive_extent<T>::value)
	{
/*
		std::cout << "Array_view Extent=" << Extent << " extent=" << extent << " Rank=" << Rank << " m_size=" << m_size
				<< " (this=" << recursive_extent<my_type>::current_level
				<< ",sub=" << recursive_extent<my_type>::sub_level << ")" << std::endl;
*/
	}
	~Array_view() { }

	//enum { m_size=recursive_extent<my_type>::value };

	inline real_t* data() { return m_data; }
	inline const real_t* data() const { return m_data; }
	inline std::size_t size() const { return m_extent; }
//	inline std::size_t extent() const { return m_extent; }

	inline value_type operator[]( std::size_t i ) { return value_type( m_data+i*recursive_extent<T>::value ); }
	inline value_type operator[]( std::size_t i ) const { return value_type( m_data+i*recursive_extent<T>::value ); }

	inline iterator begin() { return iterator( this->data() ); }
	inline iterator end() { return iterator( this->data()+m_nT ); } // one past end

	inline const_iterator cbegin() const { return const_iterator( this->data() ); }
	inline const_iterator cend() const { return const_iterator( this->data()+m_nT ); } // one past end

private:
	real_t* const m_data;
	const std::size_t m_extent;
	const std::size_t m_nT; // #of T type elements pointed by m_data
//	const std::size_t m_size;
};

template< typename T, std::size_t Extent >
class Array_view<T,Extent,1> //: public Array_view_base< Array_view<T,Extent,1> > //typename std::enable_if< std::is_floating_point<T>::value >::type > : Array_view_base< Array_view<T,Extent> >
{
public:
	using value_type = T;
	using real_t = typename scalar_type<value_type>::type;

	using my_type = Array_view<T,Extent>;

	using iterator = Array_view_iterator<value_type>;
	using const_iterator = const Array_view_iterator<value_type>;

	Array_view( real_t * const data, std::size_t extent=Extent ) : m_data(data), m_extent(extent), m_nT(std::max(extent,recursive_extent<my_type>::value))
	{
/*
		std::cout << "Array_view Extent=" << Extent << " extent=" << extent << " Rank=" << Rank << " m_size=" << m_size
				<< " (this=" << recursive_extent<my_type>::current_level
				<< ",sub=" << recursive_extent<my_type>::sub_level << ")" << std::endl;
*/
	}
	~Array_view() { }

	//enum { m_nT=recursive_extent<my_type>::value };
	enum { Rank=1 };

	inline real_t* data() { return m_data; }
	inline const real_t* data() const { return m_data; }
	inline std::size_t size() const { return m_extent; }
//	inline std::size_t extent() const { return m_extent; }

	inline value_type& operator[]( std::size_t i ) { return *(m_data+i);	}
	inline const value_type& operator[]( std::size_t i ) const { return *(m_data+i);	}

	inline iterator begin() { return iterator( this->data() ); }
	inline iterator end() { return iterator( this->data()+m_nT ); } // one past end

	inline const_iterator cbegin() const { return const_iterator( this->data() ); }
	inline const_iterator cend() const { return const_iterator( this->data()+m_nT ); } // one past end

private:
	real_t* const m_data;
	const std::size_t m_extent;
	const std::size_t m_nT; // #of T type elements pointed by m_data
	//const std::size_t m_size;
};

template< typename T, std::size_t Extent >
std::ostream& operator<< ( std::ostream& os, const Array_view<T,Extent>& array_view )
{
	for( std::size_t i = 0; i < array_view.size(); ++i )
	{
		if( i != 0 ) { os << " "; }
		os << array_view[i];
	}
	return os;
}

template< typename T, std::size_t Extent, std::size_t Extent2 >
std::ostream& operator<< ( std::ostream& os, const Array_view< Array_view<T,Extent>,Extent2>& array_view )
{
	for( std::size_t i = 0; i < array_view.size(); ++i )
	{
		if( i != 0 ) { os << " "; }
		os << array_view[i];
	}
	return os;
}

template< std::size_t M, std::size_t N, typename Extractable >
class extractor
{
public:
	extractor( const Extractable& extractable ) : m_extractable( extractable ) { }

	// from C++14 forward
	const auto operator[]( std::size_t i ) const { return m_extractable[i]; }

	std::array< std::array<double,N>, M > operator()() const
	{
		std::array< std::array<double,N>, M > m;
		for( std::size_t i = 0; i < M; ++i )
		{
			for( std::size_t j = 0; j < M; ++j )
			{
				m[i][j] = m_extractable[i][j];
			}
		}
		return m;
	}

private:
	const Extractable& m_extractable;
};

template< std::size_t M, std::size_t N, typename Extractable >
extractor<M,N,Extractable> extract( const Extractable& extractable ) { return extractor<M,N,Extractable>(extractable); }

/*
template< std::size_t ExtractN, typename Extractable >
class extractor<0,ExtractN,Extractable>
{
public:
	extractor( const Extractable& extractable ) : m_extractable( extractable ) { }

	// from C++14 forward
	const auto&& operator[]( std::size_t i ) const { return m_extractable[i]; }

private:
	const Extractable& m_extractable;
};
*/

template< std::size_t M, std::size_t N, typename Extractable >
std::ostream& operator<< ( std::ostream& os, const extractor<M,N,Extractable>& extor )
{
	for( std::size_t i = 0; i < M; ++i )
	{
		//const auto& Nextor = extor[i];
		for( std::size_t j = 0; j < N; ++j )
		{
			if( i == 0 && j == 0 ) { os << extor[i][j]; }
			else
			{
				os << " " << extor[i][j];
			}
		}
	}
	return os;
}

template< std::size_t M, std::size_t N, typename Extractable, typename ...Args >
boost::accumulators::accumulator_set< Args... >& operator<< ( boost::accumulators::accumulator_set< Args... >& acc, const Extractable& extor )
{
	for( std::size_t i = 0; i < M; ++i )
	{
		for( std::size_t j = 0; j < N; ++j )
		{
			acc( extor[i][j] );
		}
	}
	return acc;
}

} // namespace apegrunt

#include "Array_view_iterator.hpp"

namespace apegrunt {

template< typename T, std::size_t Extent >
typename Array_view<T,Extent>::iterator begin( Array_view<T,Extent>& array_view ) { return array_view.begin(); }

template< typename T, std::size_t Extent >
typename Array_view<T,Extent>::iterator end( Array_view<T,Extent>& array_view ) { return array_view.end(); }

template< typename T, std::size_t Extent >
typename Array_view<T,Extent>::const_iterator begin( const Array_view<T,Extent>& array_view ) { return array_view.cbegin(); }

template< typename T, std::size_t Extent >
typename Array_view<T,Extent>::const_iterator end( const Array_view<T,Extent>& array_view ) { return array_view.cend(); }

template< typename T, std::size_t Extent >
typename Array_view<T,Extent>::const_iterator cbegin( const Array_view<T,Extent>& array_view ) { return array_view.cbegin(); }

template< typename T, std::size_t Extent >
typename Array_view<T,Extent>::const_iterator cend( const Array_view<T,Extent>& array_view ) { return array_view.cend(); }

} // namespace apegrunt

#include "Array_view_operations.hpp"

#endif // APEGRUNT_ARRAY_VIEW_HPP
