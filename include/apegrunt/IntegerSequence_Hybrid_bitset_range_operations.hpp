/** @file IntegerSequence_Hybrid_bitset_range_operations.hpp
 
	Copyright (c) 2018-2019 Santeri Puranen. All rights reserved.
 
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

#ifndef APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_HPP
#define APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_HPP

#include <vector>
#include <ostream>

#include "IntegerSequence_forward.h"
#include "IntegerSequence_iterator.hpp"
#include "IntegerSequence_operations_forward.h"
#include "IntegerSequence_Hybrid_bitset_range_operations_forward.h"

namespace apegrunt {

// Specializations for IntegerSequence< Hybrid_bitset_range_element >

// This version is aligned at value_type::BITSTRING_SIZE chunks and does not
// have single index entries, but only bitstring and interval entries.
template< typename IndexT >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& container )
{
	os << std::hex; // use hex throughout

	for( const auto& pos: container.m_storage )	{ os << pos; }

	os << std::dec;

	return os;

}

// intersection
template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > set_intersection(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	using container_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	using value_type = Hybrid_bitset_range_element<IndexT,true>;
	using index_t = typename value_type::index_t;

	container_t intersection;

	auto aitr = cbegin(a.m_storage);
	const auto aenditr = cend(a.m_storage);

	auto bitr = cbegin(b.m_storage);
	const auto benditr = cend(b.m_storage);

	auto apos = (*aitr)();
	auto bpos = (*bitr)();

	while( aitr != aenditr && bitr != benditr )
	{
        if( apos < bpos )
        {
        	// increment apos
            if( aitr->is_range() && apos < (*aitr)() )
            {
            	apos += value_type::BITSTRING_SIZE;
            }
            else
            {
            	apos = ( ++aitr != aenditr ) ? (*aitr)() : 0;
            }
        }
        else
        {
            if( bpos < apos )
            {
            	// increment bpos
            	if( bitr->is_range() && bpos < bitr->range_end() )
            	{
            		bpos += value_type::BITSTRING_SIZE;
            	}
            	else { bpos = ( ++bitr != benditr ) ? (*bitr)() : 0; }
            }
            else // apos == bpos
            {
                if( aitr->is_range() )
                {
                	if( bitr->is_range() ) // both aitr and bitr and ranges
                	{
                		// fast-forward to the end of the shared range
                		auto range_end = std::min( aitr->range_end(), bitr->range_end() );
                		// gather weights
                		intersection.m_storage.emplace_back( apos,  range_end );
                    	// update positions
                		apos = ( range_end+1 < aitr->range_end() ) ? range_end+1 : ( ++aitr != aenditr ) ? (*aitr)() : 0;
                		bpos = ( range_end+1 < bitr->range_end() ) ? range_end+1 : ( ++bitr != benditr ) ? (*bitr)() : 0;
                 	}
                	else // aitr is a range, but bitr is not
                	{
                		// gather weights
                		intersection.m_storage.emplace_back( apos,  ~index_t(0) & bitr->get_bitfield() );
                		// update positions
                		apos = ( apos < aitr->range_end() ) ? (apos += value_type::BITSTRING_SIZE) : ( ++aitr != aenditr ) ? (*aitr)() : 0;
                		bpos = ++bitr != benditr ? (*bitr)() : 0;
                	}

                }
                else if( bitr->is_range() ) // bitr is a range, but aitr is not
                {
            		// gather weights
            		intersection.m_storage.emplace_back( apos,  ~index_t(0) & aitr->get_bitfield() );
            		// update positions
                	apos = ++aitr != aenditr ? (*aitr)() : 0;
            		bpos = ( bpos < bitr->range_end() ) ? (bpos += value_type::BITSTRING_SIZE) : ( ++bitr != benditr ) ? (*bitr)() : 0;
                }
                else // neither aitr nor bitr are ranges; equivalent to if( !aitr->is_range() && !bitr->is_range() )
                {
                	const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
                	if( bitfield )
                	{
                		// gather weights
                		intersection.m_storage.emplace_back( apos,  bitfield );
                 	}

                	// update positions
                	apos = ++aitr != aenditr ? (*aitr)() : 0;
                	bpos = ++bitr != benditr ? (*bitr)() : 0;
                }
            }

        }
    }
}

// intersect and gather for aligned container
template< typename IndexT, typename RealT >
RealT intersect_and_gather(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b,
		const std::vector<RealT>& weights
	)
{
	using container_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	using value_type = Hybrid_bitset_range_element<IndexT,true>;
	using index_t = typename value_type::index_t;
	using real_t = RealT;

	auto aitr = cbegin(a.m_storage);
	const auto aenditr = cend(a.m_storage);

	auto bitr = cbegin(b.m_storage);
	const auto benditr = cend(b.m_storage);

	auto apos = (*aitr)();
	auto bpos = (*bitr)();

	real_t weight(0);

	while( aitr != aenditr && bitr != benditr )
	{
        if( apos < bpos )
        {
        	// increment apos
            if( aitr->is_range() && apos < (*aitr)() )
            {
            	apos += value_type::BITSTRING_SIZE;
            }
            else
            {
            	apos = ( ++aitr != aenditr ) ? (*aitr)() : 0;
            }
        }
        else
        {
            if( bpos < apos )
            {
            	// increment bpos
            	if( bitr->is_range() && bpos < bitr->range_end() )
            	{
            		bpos += value_type::BITSTRING_SIZE;
            	}
            	else { bpos = ( ++bitr != benditr ) ? (*bitr)() : 0; }
            }
            else // apos == bpos
            {
                if( aitr->is_range() )
                {
                	if( bitr->is_range() ) // both aitr and bitr and ranges
                	{
                		// This appears to happen very rarely, but when it does we can
                		// fast-forward to the end of the shared range
                		auto range_end = std::min( aitr->range_end(), bitr->range_end() );

                		while( apos < range_end && bpos < range_end )
                		{
                    		// gather weights
                       		weight += value_type::gather( &weights[apos], ~index_t(0) );
                        	// update positions
                       		apos = ( apos < aitr->range_end() ) ? (apos += value_type::BITSTRING_SIZE) : ( ++aitr != aenditr ) ? (*aitr)() : 0;
                    		bpos = ( bpos < bitr->range_end() ) ? (bpos += value_type::BITSTRING_SIZE) : ( ++bitr != benditr ) ? (*bitr)() : 0;
                		}
                	}
                	else // aitr is a range, but bitr is not
                	{
                		// gather weights
                		weight += value_type::gather( &weights[apos], ~index_t(0) & bitr->get_bitfield() );
                    	// update positions
                		apos = ( apos < aitr->range_end() ) ? (apos += value_type::BITSTRING_SIZE) : ( ++aitr != aenditr ) ? (*aitr)() : 0;
                		bpos = ++bitr != benditr ? (*bitr)() : 0;
                	}

                }
                else if( bitr->is_range() ) // bitr is a range, but aitr is not
                {
            		// gather weights
            		weight += value_type::gather( &weights[apos], ~index_t(0) & aitr->get_bitfield() );
            		// update positions
                	apos = ++aitr != aenditr ? (*aitr)() : 0;
            		bpos = ( bpos < bitr->range_end() ) ? (bpos += value_type::BITSTRING_SIZE) : ( ++bitr != benditr ) ? (*bitr)() : 0;
                }
                else // neither aitr nor bitr are ranges; equivalent to if( !aitr->is_range() && !bitr->is_range() )
                {
                	const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
                	if( bitfield )
                	{
                		// gather weights
                		weight += value_type::gather( &weights[apos], bitfield );
                 	}

                	// update positions
                	apos = ++aitr != aenditr ? (*aitr)() : 0;
                	bpos = ++bitr != benditr ? (*bitr)() : 0;
                }
            }

        }
    }

	return weight;
}

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_HPP
