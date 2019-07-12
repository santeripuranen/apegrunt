/** @file IntegerSequence_interface.hpp
 
	Copyright (c) 2016-2019 Santeri Puranen. All rights reserved.
 
	By installing, copying or otherwise using the attached
	material ("product" or "software") you acknowledge and
	agree that the attached	material contains proprietary
	information of the copyright holder(s). Any use of the
	material is prohibited except as expressly agreed between
	the copyright holder(s) and the recipient.
 
	THIS PRODUCT ("SOFTWARE") IS PROVIDED "AS IS", WITHOUT WARRANTY
	OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
	THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
	PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
	COPYRIGHT HOLDER(S) BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY,
	WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
	IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
	THE SOFTWARE.

	@author Santeri Puranen
	$Id: $
*/

#ifndef APEGRUNT_INTEGERSEQUENCE_INTERFACE_HPP
#define APEGRUNT_INTEGERSEQUENCE_INTERFACE_HPP

#include <vector>

//#define BOOST_POOL_NO_MT
//#include <boost/pool/pool_alloc.hpp>

#include "IntegerSequence_forward.h"
#include "IntegerSequence_operations_forward.h"
#include "IntegerSequence_iterator.hpp"
#include "aligned_allocator.hpp"

namespace apegrunt {

template< typename RangeT >
class IntegerSequence
{
public:
	using range_t = RangeT;
	using my_type = IntegerSequence<range_t>;
	using value_type = typename range_t::index_t;
	using aligned = typename range_t::aligned;
	using const_iterator = apegrunt::iterator::IntegerSequence_const_iterator<range_t>;
	using iterator = const_iterator;

	IntegerSequence() : m_back_inserter( m_storage ), m_size(0) { }
	//IntegerSequence() : m_back_inserter( apegrunt::back_inserter(m_storage) ), m_size(0) { }
	IntegerSequence( const value_type& value ) : m_back_inserter( m_storage ), m_size(0) { this->push_back(value); }
	//
	IntegerSequence( const my_type& other )
	: m_storage( other.m_storage ),
	  m_back_inserter( other.m_back_inserter, m_storage ), // need to supply reference to the _copied_ m_storage, when copy-constructing back_inserter. The regular back_inserter copy constructor, would copy a reference to the old m_storage. This is an issue for example when reallocating std::vector<IntegerSequence>.
	  m_size(other.m_size)
	{ }
	IntegerSequence( my_type&& other ) : m_storage( std::move(other.m_storage) ), m_back_inserter( std::move(other.m_back_inserter) ), m_size( std::move(other.m_size) ) { }
	~IntegerSequence()// = default;
	{
		//boost::singleton_pool<boost::pool_allocator_tag, sizeof(range_t)>::purge_memory();
	}

	//IntegerSequence( my_type&& other ) noexcept : m_storage( std::move(other.m_storage) ), m_back_inserter( std::move(other.m_back_inserter) ), m_size(other.m_size) { }

	inline void push_back( const value_type& value ) { m_back_inserter = value; ++m_size; }

	inline std::size_t size() const { return m_size; }
	inline std::size_t bytesize() const { return sizeof(range_t) * m_storage.size(); }

	const_iterator cbegin() const { return const_iterator( m_storage.data(), this->one_past_end() ); }
	const_iterator cend() const { return const_iterator( this->one_past_end(), this->one_past_end() ); }

	iterator begin() const { return this->cbegin(); }
	iterator end() const { return this->cend(); }

private:
	/*
		using allocator_t = boost::pool_allocator<range_t,
		    boost::default_user_allocator_new_delete,
		    boost::details::pool::default_mutex,
		    0xFF, 0xFFFF>;
	*/
	// the aligned allocator is sometimes slow -- don't use it unless it's really needed (= aligned loads and SIMD needed to make downstream code significantly faster)
	//using allocator_t = memory::AlignedAllocator<range_t,alignof(range_t)>;
	//using container_t = std::vector< range_t, allocator_t >;
	using container_t = std::vector< range_t >;
	container_t m_storage;
	apegrunt::iterator::back_insert_iterator< container_t > m_back_inserter;
	std::size_t m_size;

	inline const range_t * const one_past_end() const { return &(m_storage.back())+1; } // this pointer should never be dereferenced

	friend std::ostream& operator<< <value_type> ( std::ostream& os, const my_type& container );
	friend my_type set_intersection<value_type,aligned::value>( const my_type& a, const my_type& b );
	friend double intersect_and_gather<value_type,double>( const my_type& a, const my_type& b, const std::vector<double>& weights );
};

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_INTERFACE_HPP
