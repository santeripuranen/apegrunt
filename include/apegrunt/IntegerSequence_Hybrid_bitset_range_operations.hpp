/** @file IntegerSequence_Hybrid_bitset_range_operations.hpp
 
	Copyright (c) 2018-2021 Santeri Puranen.

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

#ifndef APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_HPP
#define APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_HPP

#include <vector>
#include <ostream>
#include <tuple>
#include <utility>

#include "IntegerSequence_forward.h"
#include "IntegerSequence_iterator.hpp"
#include "IntegerSequence_operations_forward.h"

namespace apegrunt {

// This version is aligned at value_type::BITSTRING_SIZE chunks and does not
// have single index entries, but only bitstrings and interval entries.
template< typename IndexT >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& container )
{
	os << std::hex; // use hex throughout

	bool first = true;
	for( const auto& pos: container.storage() )	{ if( first ) { first = false; os << pos; } else { os << "," << pos; } }

	os << std::dec;

	return os;

}

template< typename IndexT >
inline Hybrid_bitset_range_element<IndexT,true> fuse_range(
		const Hybrid_bitset_range_element<IndexT,true>& a,
		const Hybrid_bitset_range_element<IndexT,true>& b
	)
{
	using element_t = Hybrid_bitset_range_element<IndexT,true>;

	// check overlap
	if( a.range_end() < b() || b.range_end() < a() ) { return element_t(); } // elements do not overlap

	return element_t().set_range( std::min(a(),b()), std::max(a.range_end(),b.range_end()) );
}

template< typename IndexT >
inline Hybrid_bitset_range_element<IndexT,true> range_intersection(
		const Hybrid_bitset_range_element<IndexT,true>& a,
		const Hybrid_bitset_range_element<IndexT,true>& b
	)
{
	using element_t = Hybrid_bitset_range_element<IndexT,true>;

	// check overlap

	return element_t().set_range( std::max(a(),b()), std::min(a.range_end(),b.range_end()) );
}

template< typename IndexT >
std::pair<IndexT,IndexT> range_intersection_endpoints(
		const Hybrid_bitset_range_element<IndexT,true>& a,
		const Hybrid_bitset_range_element<IndexT,true>& b
	)
{
	const auto intersection( range_intersection(a,b) );
	return intersection ? std::make_pair(intersection.range_begin(),intersection.range_end()) : std::make_pair(IndexT(0),IndexT(0));
}

// gather for aligned container
template< typename IndexT, typename RealT >
RealT gather( const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a, const RealT *weights )
{
	using std::cbegin; using std::cend;
	using value_type = Hybrid_bitset_range_element<IndexT,true>;

	return std::accumulate( cbegin(a.storage()), cend(a.storage()),	RealT(0),
			[weights]( RealT sum, const auto& element ) { return sum + ( element.is_range() ? value_type::gather( weights, element ) : value_type::gather( weights+element.range_begin(), element.get_bitfield() ) ); }
	);
}

template< typename IndexT >
inline void range_difference(
		const Hybrid_bitset_range_element<IndexT,true>& a,
		const Hybrid_bitset_range_element<IndexT,true>& b,
		const Hybrid_bitset_range_element<IndexT,true>& intersection,
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& anb,
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& bna
	)
{
	using element_t = Hybrid_bitset_range_element<IndexT,true>;

	const auto begin = intersection.range_begin();
	const auto end = intersection.range_end();

	// add any range content that's lower-valued than the intersection (it's in either a, b or in neither by design)
	if( a.range_begin() < begin )
	{
		anb.emplace_back( std::move( element_t( a.range_begin() ).set_range_end( begin ) ) );
	}
	else
	{
		if( b.range_begin() < begin )
		{
			bna.emplace_back( std::move( element_t( b.range_begin() ).set_range_end( begin ) ) );
		}
	}

	// add any range content that's higher-valued than the intersection (it's in either a, b or in neither by design)
	if( end < a.range_end() )
	{
		anb.emplace_back( std::move( element_t( end ).set_range_end( a.range_end() ) ) );
	}
	else
	{
		if( end < b.range_end() )
		{
			bna.emplace_back( std::move( element_t( end ).set_range_end( b.range_end() ) ) );
		}
	}
}

template< typename IndexT >
inline void pre_intersection_range_difference(
		const Hybrid_bitset_range_element<IndexT,true>& a,
		const Hybrid_bitset_range_element<IndexT,true>& b,
		const Hybrid_bitset_range_element<IndexT,true>& intersection,
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& anb,
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& bna
	)
{
	using element_t = Hybrid_bitset_range_element<IndexT,true>;

	const auto begin = intersection.range_begin();
	const auto end = intersection.range_end();

	// add any range content that's lower-valued than the intersection (it's in either a, b or in neither by design)
	if( a.range_begin() < begin )
	{
		anb.emplace_back( std::move( element_t( a.range_begin() ).set_range_end( begin ) ) );
	}
	else
	{
		if( b.range_begin() < begin )
		{
			bna.emplace_back( std::move( element_t( b.range_begin() ).set_range_end( begin ) ) );
		}
	}
}

template< typename IndexT >
inline bool inplace_range_difference_low(
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& aitr,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& anb,
		IndexT begin
	)
{
	bool modified = false;

	if( anb >= aitr ) {
		anb = a.storage().emplace( anb );
		aitr = anb+1;
		modified = true;
	}
	anb->set_range( aitr->range_begin(), begin );
	a.append_block_mask(*anb);
	++anb;

	return modified;
}

template< typename IndexT >
inline bool inplace_range_difference_high(
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& aitr,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& anb,
		IndexT end
	)
{
	bool modified = false;

	if( anb >= aitr ) {
		anb = a.storage().emplace( anb );
		aitr = anb+1;
		modified = true;
	}
	anb->set_range( end, aitr->range_end() );
	a.append_block_mask(*anb);
	++anb;

	return modified;
}

template< typename IndexT >
inline bool inplace_range_difference(
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& aitr,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& bitr,
		const Hybrid_bitset_range_element<IndexT,true>& intersection,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& anb,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& bna
	)
{
	bool modified = false;

	const auto begin = intersection.range_begin();
	const auto end = intersection.range_end();

	// add any range content that's lower-valued than the intersection (it's in either a, b or in neither by design)
	if( aitr->range_begin() < begin )
	{
		modified = inplace_range_difference_low(a,aitr,anb,begin);
	}
	else
	{
		if( bitr->range_begin() < begin )
		{
			modified = inplace_range_difference_low(b,bitr,bna,begin);
		}
	}

	// add any range content that's higher-valued than the intersection (it's in either a, b or in neither by design)
	if( end < aitr->range_end() )
	{
		modified = inplace_range_difference_high(a,aitr,anb,end);
	}
	else
	{
		if( end < bitr->range_end() )
		{
			modified = inplace_range_difference_high(b,bitr,bna,end);
		}
	}

	return modified;
}

template< typename IndexT >
inline bool pre_intersection_inplace_range_difference(
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& aitr,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& bitr,
		const Hybrid_bitset_range_element<IndexT,true>& intersection,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& anb,
		typename IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >::container_t::iterator& bna
	)
{
	const auto begin = intersection.range_begin();
	const auto end = intersection.range_end();

	bool modified = false;

	// add any range content that's lower-valued than the intersection (it's in either a, b or in neither by design)
	if( aitr->range_begin() < begin )
	{
		modified = inplace_range_difference_low(a,aitr,anb,begin);
	}
	else
	{
		if( bitr->range_begin() < begin )
		{
			modified = inplace_range_difference_low(b,bitr,bna,begin);
		}
	}
// /*
	// adjust any remaining higher-valued range content to begin where the intersection ends (it's in either a, b or in neither by design)
	if( end < aitr->range_end() )
	{
		aitr->set_range( end, aitr->range_end() );
		//++bitr;
	}
	else
	{
		//++aitr;
		if( end < bitr->range_end() ) {
			bitr->set_range( end, bitr->range_end() );
		}
		else
		{
			//++bitr;
		}
	}
// */
	return modified;
}

template< typename IteratorT, typename ValueT, typename CompT >
inline bool advance_while_false( IteratorT& first, const IteratorT& last, const ValueT value, CompT comp=std::greater<ValueT>() )
{
	if( !comp(first->range_end(),value) ) { do { ++first; } while( first != last && !comp(first->range_end(),value) ); return first != last; }
	else { return false; }
}

template< typename IteratorT, typename ValueT, typename CompT >
inline bool advance_while_true( IteratorT& first, const IteratorT& last, const ValueT value, CompT comp=std::less<ValueT>() )
{
	if( comp(first->range_end(),value) ) { do { ++first; } while( first != last && comp(first->range_end(),value) ); return true; }
	else { return false; }
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

	container_t result;

	// nothing to intersect?
	if( a.has_overlap(b) ) // this test should suffice; covers cases where either a or b, or both a and b, are empty
	{
		auto aitr = cbegin(a.storage());
		const auto aenditr = cend(a.storage());

		auto bitr = cbegin(b.storage());
		const auto benditr = cend(b.storage());

		while( aitr != aenditr && bitr != benditr )
		{
			//if( aitr->range_end() <= bitr->range_begin() ) { do { ++aitr; } while( aitr != aenditr && aitr->range_end() <= bitr->range_begin() ); }
			if( aitr->range_end() <= bitr->range_begin() ) { while( ++aitr != aenditr && aitr->range_end() <= bitr->range_begin() ); }
			//if( aitr->range_end() <= bitr->range_begin() ) { while( ++aitr != aenditr && aitr->range_end() < bitr->range_begin() ); }
			else
			{
				//if( bitr->range_end() <= aitr->range_begin() ) { do { ++bitr; } while( bitr != benditr && bitr->range_end() <= aitr->range_begin() ); }
				if( bitr->range_end() <= aitr->range_begin() ) { while( ++bitr != benditr && bitr->range_end() <= aitr->range_begin() ); }
				//if( bitr->range_end() <= aitr->range_begin() ) { while( ++bitr != benditr && bitr->range_end() < aitr->range_begin() ); }
				else // *aitr and *bitr overlap in some way
				{
					if( aitr->is_range() )
					{
						if( bitr->is_range() ) // both aitr and bitr are ranges
						{
							// fast-forward to the end of shared range
							auto intersection( range_intersection( *aitr, *bitr ) ); // will always succeed, since overlap is guaranteed in this branch. There will by design not be overlap with anything already in storage.

							// update positions
							if( aitr->range_end() == intersection.range_end() ) { ++aitr; }
							if( bitr->range_end() == intersection.range_end() ) { ++bitr; }

							result.push_back( std::move(intersection) );
						}
						else // aitr is a range, but bitr is not
						{
							result.push_back(*bitr);

							// update positions
							if( aitr->range_end() == bitr->range_end() ) { ++aitr; }
							++bitr;
						}
					}
					else
					{
						if( bitr->is_range() ) // bitr is a range, but aitr is not
						{
							result.push_back(*aitr);

							// update positions
							if( bitr->range_end() == aitr->range_end() ) { ++bitr; }
							++aitr;
						}
						else // neither aitr nor bitr are ranges, so (*aitr)() == (*bitr)() must hold; equivalent to if( !aitr->is_range() && !bitr->is_range() )
						{
							const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
							if( bitfield ) { result.emplace_back( value_type(*aitr, bitfield) ); }

							// update positions
							++aitr;
							++bitr;
						}
					}
				}
			}
		}
		result.shrink_to_fit();
	}

	return result;
}

// fast-ish check whether two sets intersect
template< typename IndexT >
bool has_set_intersection(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	using std::cbegin; using std::cend;

	// nothing intersect?
	if( !a.has_overlap(b) ) { return false; } // this test should suffice; covers cases where either a or b, or both a and b, are empty

	auto aitr = cbegin(a.storage());
	const auto aenditr = cend(a.storage());

	auto bitr = cbegin(b.storage());
	const auto benditr = cend(b.storage());

	while( aitr != aenditr && bitr != benditr )
	{
        if( !(aitr->range_end() > (*bitr)()) ) { do { ++aitr; } while( aitr != aenditr && !(aitr->range_end() > (*bitr)())); }
        else
		//if( !advance_while_false( aitr, aenditr, (*bitr)(), std::greater<index_t>() ) ) // { } // linear search
        {
            if( !(bitr->range_end() > (*aitr)()) ) { do { ++bitr; } while( bitr != benditr && !(bitr->range_end() > (*aitr)()) ); }
            else // *aitr and *bitr /might/ overlap in some way
    		//if( !advance_while_false( bitr, benditr, (*aitr)(), std::greater<index_t>() ) ) // { } // linear search
            {
            	if( aitr->is_range() || bitr->is_range() )
                {
            		// overlap is guaranteed when either aitr or bitr (or both) is a range
            		return true;
                }
				else // neither aitr nor bitr are ranges, so (*aitr)() == (*bitr)() must hold; equivalent to if( !aitr->is_range() && !bitr->is_range() )
				{
					// we need to test whether there is overlap or not
					if( aitr->get_bitfield() & bitr->get_bitfield() ) { return true; }

					// no overlap has been found so far: update positions and do another iteration
					++aitr;
					++bitr;
				}
            }
        }
    }
	return false;
}

// return the size of the intersection
template< typename IndexT >
std::size_t set_intersection_size(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	using std::cbegin; using std::cend;

	std::size_t size(0);

	// nothing to intersect?
	if( !a.has_overlap(b) ) { return size; } // this test should suffice; covers cases where either a or b, or both a and b, are empty

	auto aitr = cbegin(a.storage());
	const auto aenditr = cend(a.storage());

	auto bitr = cbegin(b.storage());
	const auto benditr = cend(b.storage());

	while( aitr != aenditr && bitr != benditr )
	{
        if( !(aitr->range_end() > (*bitr)()) ) { do { ++aitr; } while( aitr != aenditr && !(aitr->range_end() > (*bitr)())); }
        else
        {
        	if( !(bitr->range_end() > (*aitr)()) ) { do { ++bitr; } while( bitr != benditr && !(bitr->range_end() > (*aitr)()) ); }
            else // *aitr and *bitr overlap in some way
            {
            	if( aitr->is_range() )
                {
                	if( bitr->is_range() ) // both aitr and bitr are ranges
                	{
                        // fast-forward to the end of shared range
                		const auto intersection( range_intersection( *aitr, *bitr ) ); // will always succeed, since overlap is guaranteed in this branch. There will by design not be overlap with anything already in storage.
                		size += intersection.size();

                		// update positions
                		if( !(aitr->range_end() > intersection.range_end() ) ) { ++aitr; }
                		if( !(bitr->range_end() > intersection.range_end() ) ) { ++bitr; }
                 	}
                	else // aitr is a range, but bitr is not
                	{
                		size += bitr->size();

                		// update positions
                		if( ++bitr != benditr && !(aitr->range_end() > (*bitr)()) ) { ++aitr; }
                	}
                }
            	else
            	{
					if( bitr->is_range() ) // bitr is a range, but aitr is not
					{
						size += aitr->size();

						// update positions
						if( ++aitr != aenditr && !(bitr->range_end() > (*aitr)()) ) { ++bitr; }
					}
					else // neither aitr nor bitr are ranges, so (*aitr)() == (*bitr)() must hold; equivalent to if( !aitr->is_range() && !bitr->is_range() )
					{
						size += apegrunt::popcnt( aitr->get_bitfield() & bitr->get_bitfield() );

						// update positions
						++aitr;
						++bitr;
					}
            	}
            }
        }
    }
	return size;
}

// convenience/syntactic sugar
template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > operator&(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	return set_intersection(a,b);
}

// intersect and gather for aligned container
template< typename IndexT, typename UnaryFunction >
UnaryFunction for_each_in_set_intersection(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b,
		UnaryFunction f
	)
{
	using value_type = Hybrid_bitset_range_element<IndexT,true>;
	using std::cbegin; using std::cend;

	// nothing to intersect?
	if( a.has_overlap(b) ) // this test should suffice; covers cases where either a or b, or both a and b, are empty
	{
		auto aitr = cbegin(a.storage());
		const auto aenditr = cend(a.storage());

		auto bitr = cbegin(b.storage());
		const auto benditr = cend(b.storage());

		while( aitr != aenditr && bitr != benditr )
		{
			if( aitr->range_end() <= bitr->range_begin() ) { do { ++aitr; } while( aitr != aenditr && aitr->range_end() <= bitr->range_begin() ); }
			else
			{
				if( bitr->range_end() <= aitr->range_begin() ) { do { ++bitr; } while( bitr != benditr && bitr->range_end() <= aitr->range_begin() ); }
				else // *aitr and *bitr overlap in some way
				{
					if( aitr->is_range() )
					{
						if( bitr->is_range() ) // both aitr and bitr are ranges
						{
							// fast-forward to the end of shared range
							const auto intersection( range_intersection( *aitr, *bitr ) ); // will always succeed, since overlap is guaranteed in this branch. There will by design not be overlap with anything already in storage.
							f(intersection);

							// update positions
							if( aitr->range_end() == intersection.range_end() ) { ++aitr; }
							if( bitr->range_end() == intersection.range_end() ) { ++bitr; }
						}
						else // aitr is a range, but bitr is not
						{
							f(*bitr);

							// update positions
	                		if( aitr->range_end() == bitr->range_end() ) { ++aitr; }
							++bitr;
						}
					}
					else
					{
						if( bitr->is_range() ) // bitr is a range, but aitr is not
						{
							f(*aitr);

							// update positions
	                		if( bitr->range_end() == aitr->range_end() ) { ++bitr; }
							++aitr;
						}
						else // neither aitr nor bitr are ranges, so (*aitr)() == (*bitr)() must hold; equivalent to if( !aitr->is_range() && !bitr->is_range() )
						{
							const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
							if( bitfield ) { f( value_type( *aitr, bitfield ) ); }

							// update positions
							++aitr;
							++bitr;
						}
					}
				}
			}
		}
	}
	return f;
}

// set differences (a-not-b and b-not-a) for aligned container
template< typename IndexT >
std::pair< IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >, IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > >
set_differences(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	using container_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	using return_t = std::pair<container_t,container_t>;
	using value_type = Hybrid_bitset_range_element<IndexT,true>;

	using std::cbegin; using std::cend;

	// nothing to intersect?
	if( !a.has_overlap(b) ) { return std::make_pair(a,b); } // this test should suffice; covers cases where either a or b, or both a and b, are empty

	auto aitr = cbegin(a.storage());
	const auto aenditr = cend(a.storage());

	auto bitr = cbegin(b.storage());
	const auto benditr = cend(b.storage());

	return_t result;

	auto& anb = result.first;
	auto& bna = result.second;

	while( aitr != aenditr && bitr != benditr )
	{
        if( !(aitr->range_end() > (*bitr)()) ) { do { anb.push_back(*aitr); ++aitr; } while( aitr != aenditr && !(aitr->range_end() > (*bitr)())); }
        else
        {
            if( !(bitr->range_end() > (*aitr)()) ) { do { bna.push_back(*bitr); ++bitr; } while( bitr != benditr && !(bitr->range_end() > (*aitr)()) ); }
            else // *aitr and *bitr overlap in some way
            {
            	if( aitr->is_range() )
                {
                	if( bitr->is_range() ) // both aitr and bitr are ranges
                	{
                        // find shared range
                		const auto intersection( range_intersection( *aitr, *bitr ) ); // will always succeed, since overlap is guaranteed in this branch. There will by design not be overlap with anything already in storage.

                		pre_intersection_range_difference( *aitr, *bitr, intersection, anb, bna );

						// update positions
                		if( aitr->range_end() == intersection.range_end() ) { ++aitr; }
                		if( bitr->range_end() == intersection.range_end() ) { ++bitr; }
               	}
                	else // aitr is a range, but bitr is not
                	{
                		anb.emplace_back( std::move( value_type( (*bitr)(), ~bitr->get_bitfield() ) ) );

                		if( bitr->range_end() < aitr->range_end() )
                		{
                			anb.emplace_back( std::move( value_type( bitr->range_end() ).set_range_end( aitr->range_end() ) ) );
                		}

						// update positions
                 		if( ++bitr != benditr && !(aitr->range_end() > (*bitr)()) ) { ++aitr; }
                	}
                }
            	else
            	{
					if( bitr->is_range() ) // bitr is a range, but aitr is not
					{
						bna.emplace_back( value_type( (*aitr)(), ~aitr->get_bitfield() ) );

                		if( aitr->range_end() < bitr->range_end() )
                		{
                			bna.emplace_back( std::move( value_type( aitr->range_end() ).set_range_end( bitr->range_end() ) ) );
                		}

						// update positions
						if( ++aitr != aenditr && !(bitr->range_end() > (*aitr)()) ) { ++bitr; }
					}
					else // neither aitr nor bitr are ranges and (*aitr)() == (*bitr)(); equivalent to if( !aitr->is_range() && !bitr->is_range() )
					{
						const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
						if( bitfield ) // some overlap
						{
							const auto bfa = aitr->get_bitfield() & ~bitfield;
							if( bfa ) { anb.emplace_back( value_type( (*aitr)(), bfa ) ); }

							const auto bfb = bitr->get_bitfield() & ~bitfield;
							if( bfb ) { bna.emplace_back( value_type( (*bitr)(), bfb ) ); }
						}
						else // no overlap
						{
							anb.push_back(*aitr);
							bna.push_back(*bitr);
						}

						// update positions
						++aitr;
						++bitr;
					}
            	}
            }
        }
    }

	// add anything that wasn't already processed above
	for( ; aitr != aenditr; ++aitr ) { anb.push_back(*aitr); }
	anb.shrink_to_fit();

	for( ; bitr != benditr; ++bitr ) { bna.push_back(*bitr); }
	bna.shrink_to_fit();

	return result;
}

// intersect and relative complements for aligned container
template< typename IndexT >
bool inplace_set_differences(
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	using value_type = Hybrid_bitset_range_element<IndexT,true>;

	using std::cbegin; using std::cend;

	// nothing to overlap?
	if( !a.has_overlap(b) ) { return false; } // this test should suffice; covers cases where either a or b, or both a and b, are empty

	auto aitr = begin(a.storage());
	const auto aenditr = end(a.storage());

	auto bitr = begin(b.storage());
	const auto benditr = end(b.storage());

	auto anb = aitr; // a not b
	auto bna = bitr; // b not a

	while( aitr != aenditr && bitr != benditr )
	{
        if( !(aitr->range_end() > (*bitr)()) ) { do { std::swap(*anb,*aitr); ++anb; ++aitr; } while( aitr != aenditr && !(aitr->range_end() > (*bitr)())); }
        else
        {
            if( !(bitr->range_end() > (*aitr)()) ) { do { std::swap(*bna,*bitr); ++bna; ++bitr; } while( bitr != benditr && !(bitr->range_end() > (*aitr)()) ); }
            else // *aitr and *bitr overlap in some way
            {
            	if( aitr->is_range() )
                {
                	if( bitr->is_range() ) // both aitr and bitr are ranges
                	{
                        // fast-forward to the end of shared range
                		const auto intersection( range_intersection( *aitr, *bitr ) ); // will always succeed, since overlap is guaranteed in this branch. There will by design not be overlap with anything already in storage.

                		pre_intersection_inplace_range_difference( a, b, aitr, bitr, intersection, anb, bna );

						// update positions
                		if( !(aitr->range_end() > intersection.range_end() ) ) { ++aitr; }
                		if( !(bitr->range_end() > intersection.range_end() ) ) { ++bitr; }
                	}
                	else // aitr is a range, but bitr is not
                	{
                		*anb = value_type( (*bitr)(), ~bitr->get_bitfield() );
                		++anb;

                		if( bitr->range_end() < aitr->range_end() )
                		{
                			if( anb == aitr ) { anb = a.storage().emplace( ++anb, value_type() ); aitr = anb+1; }
                			*anb = std::move( value_type( bitr->range_end() ).set_range_end( aitr->range_end() ) );
                			++anb;
                		}

						// update positions
                 		if( ++bitr != benditr && !(aitr->range_end() > (*bitr)()) ) { ++aitr; }
                	}
                }
            	else
            	{
					if( bitr->is_range() ) // bitr is a range, but aitr is not
					{
						*bna = value_type( (*aitr)(), ~aitr->get_bitfield() );
						++bna;

                		if( aitr->range_end() < bitr->range_end() )
                		{
                			if( bna == bitr ) { bna = b.storage().emplace( ++bna, value_type() ); bitr = bna+1; }
                			*bna = std::move( value_type( aitr->range_end() ).set_range_end( bitr->range_end() ) );
                			++bna;
                		}

						// update positions
						if( ++aitr != aenditr && !(bitr->range_end() > (*aitr)()) ) { ++bitr; }
					}
					else // neither aitr nor bitr are ranges and (*aitr)() == (*bitr)(); equivalent to if( !aitr->is_range() && !bitr->is_range() )
					{
						const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
						if( bitfield ) // some overlap
						{
							const auto bfa = aitr->get_bitfield() & ~bitfield;
							if( bfa ) { *anb = value_type( (*aitr)(), bfa ); ++anb; }

							const auto bfb = bitr->get_bitfield() & ~bitfield;
							if( bfb ) { *bna = value_type( (*bitr)(), bfb ); ++bna; }
						}
						else // no overlap
						{
							std::swap( *anb, *aitr ); ++anb;
							std::swap( *bna, *bitr ); ++bna;
						}

						// update positions
						++aitr;
						++bitr;
					}
            	}
            }
        }
    }

	// add anything that wasn't already processed above
	for( ; aitr != aenditr; ++aitr ) { std::swap(*anb,*aitr); ++anb; }
	a.storage().erase( anb, aenditr );
	//anb.shrink_to_fit();

	for( ; bitr != benditr; ++bitr ) { std::swap(*bna,*bitr); ++bna; }
	b.storage().erase( bna, benditr );
	//bna.shrink_to_fit();

	return true;
}

// intersect and relative complements for aligned container
template< typename IndexT >
std::tuple< IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >,
			IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >,
			IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > >
set_intersection_and_differences(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	using container_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	//using return_t = std::tuple< container_t, container_t, container_t >;
	using value_type = Hybrid_bitset_range_element<IndexT,true>;
	//using index_t = typename value_type::value_type;

	using std::cbegin; using std::cend;

	// nothing to gather?
	if( !a.has_overlap(b) ) { return std::make_tuple(a,b,container_t()); } // this test should suffice; covers cases where either a or b, or both a and b, are empty

	auto result = std::make_tuple(container_t(),container_t(),container_t());

	auto aitr = cbegin(a.storage());
	const auto aenditr = cend(a.storage());

	auto bitr = cbegin(b.storage());
	const auto benditr = cend(b.storage());

	auto& anb = std::get<0>(result); // a not b
	auto& bna = std::get<1>(result);
	auto& ab = std::get<2>(result);

	while( aitr != aenditr && bitr != benditr )
	{
        if( !(aitr->range_end() > (*bitr)()) ) { do { anb.push_back(*aitr); ++aitr; } while( aitr != aenditr && !(aitr->range_end() > (*bitr)())); }
        else
        {
            if( !(bitr->range_end() > (*aitr)()) ) { do { bna.push_back(*bitr); ++bitr; } while( bitr != benditr && !(bitr->range_end() > (*aitr)()) ); }
            else // *aitr and *bitr overlap in some way
            {
            	if( aitr->is_range() )
                {
                	if( bitr->is_range() ) // both aitr and bitr are ranges
                	{
                        // fast-forward to the end of shared range
                		auto intersection( range_intersection( *aitr, *bitr ) ); // will always succeed, since overlap is guaranteed in this branch. There will by design not be overlap with anything already in storage.

                 		range_difference( *aitr, *bitr, intersection, anb, bna );

						// update positions
                		if( !(aitr->range_end() > intersection.range_end() ) ) { ++aitr; }
                		if( !(bitr->range_end() > intersection.range_end() ) ) { ++bitr; }

                		ab.emplace_back( std::move(intersection) );
                	}
                	else // aitr is a range, but bitr is not
                	{
                		ab.push_back(*bitr);

                		anb.emplace_back( std::move( value_type( (*bitr)(), ~bitr->get_bitfield() ) ) );

                		if( bitr->range_end() < aitr->range_end() )
                		{
                			anb.emplace_back( std::move( value_type( bitr->range_end() ).set_range_end( aitr->range_end() ) ) );
                		}

						// update positions
                 		if( ++bitr != benditr && !(aitr->range_end() > (*bitr)()) ) { ++aitr; }
                	}
                }
            	else
            	{
					if( bitr->is_range() ) // bitr is a range, but aitr is not
					{
						ab.push_back(*aitr);

						bna.emplace_back( value_type( (*aitr)(), ~aitr->get_bitfield() ) );

                		if( aitr->range_end() < bitr->range_end() )
                		{
                			bna.emplace_back( std::move( value_type( aitr->range_end() ).set_range_end( bitr->range_end() ) ) );
                		}

						// update positions
						if( ++aitr != aenditr && !(bitr->range_end() > (*aitr)()) ) { ++bitr; }
					}
					else // neither aitr nor bitr are ranges and (*aitr)() == (*bitr)(); equivalent to if( !aitr->is_range() && !bitr->is_range() )
					{
						const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
						//const index_t bitfield = aitr->get_bitfield() & bitr->get_bitfield();
						if( bitfield ) // some overlap
						{
							ab.emplace_back( value_type((*aitr)(),bitfield) );

							const auto bfa = aitr->get_bitfield() & ~bitfield;
							//index_t bfa = aitr->get_bitfield() & ~bitfield;
							if( bfa ) { 	anb.emplace_back( value_type( (*aitr)(), bfa ) ); }

							const auto bfb = bitr->get_bitfield() & ~bitfield;
							//index_t bfb = bitr->get_bitfield() & ~bitfield;
							if( bfb ) { bna.emplace_back( value_type( (*bitr)(), bfb ) ); }
						}
						else // no overlap
						{
							anb.push_back(*aitr);
							bna.push_back(*bitr);
						}

						// update positions
						++aitr;
						++bitr;
					}
            	}
            }
        }
    }

	// add anything that wasn't already processed above
	for( ; aitr != aenditr; ++aitr ) { anb.push_back(*aitr); }
	for( ; bitr != benditr; ++bitr ) { bna.push_back(*bitr); }

	ab.shrink_to_fit();
	anb.shrink_to_fit();
	bna.shrink_to_fit();

	return result;
}

// gather intersection
template< typename IndexT >
struct gather_intersection
{
	using intersection_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	template< typename ElementT >
	inline void operator()( const ElementT& element ) { intersection.push_back(element); }

	inline intersection_t&& get() { return std::move(intersection); }
	inline bool empty() const { return intersection.empty(); }

	IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > intersection;
};


// intersect and relative complements for aligned container
template< typename IndexT, typename UnaryFunction >
//set_differences_holder< IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > > set_differences(
std::tuple< IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >, IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >, UnaryFunction > set_differences_and_for_each_in_intersection(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b,
		UnaryFunction _f
	)
{
	using container_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	//using return_t = set_differences_holder<container_t>;
	//using return_t = std::tuple<container_t,container_t,UnaryFunction>;
	using value_type = Hybrid_bitset_range_element<IndexT,true>;
	//using index_t = typename value_type::value_type;

	using std::cbegin; using std::cend;

	// nothing to gather?
	if( !a.has_overlap(b) ) { return std::make_tuple(a,b,_f); } // this test should suffice; covers cases where either a or b, or both a and b, are empty

	auto aitr = cbegin(a.storage());
	const auto aenditr = cend(a.storage());

	auto bitr = cbegin(b.storage());
	const auto benditr = cend(b.storage());

	auto result = std::make_tuple(container_t(),container_t(),_f);

	auto& anb = std::get<0>(result); // a not b
	auto& bna = std::get<1>(result); // b not a
	auto& f = std::get<2>(result);

	while( aitr != aenditr && bitr != benditr )
	{
        if( !(aitr->range_end() > (*bitr)()) ) { do { anb.push_back(*aitr); ++aitr; } while( aitr != aenditr && !(aitr->range_end() > (*bitr)())); }
        else
        {
            if( !(bitr->range_end() > (*aitr)()) ) { do { bna.push_back(*bitr); ++bitr; } while( bitr != benditr && !(bitr->range_end() > (*aitr)()) ); }
            else // *aitr and *bitr overlap in some way
            {
            	if( aitr->is_range() )
                {
                	if( bitr->is_range() ) // both aitr and bitr are ranges
                	{
                        // fast-forward to the end of shared range
                		const auto intersection( range_intersection( *aitr, *bitr ) ); // will always succeed, since overlap is guaranteed in this branch. There will by design not be overlap with anything already in storage.

                		f(intersection);

                		pre_intersection_range_difference( *aitr, *bitr, intersection, anb, bna );

						// update positions
                		//if( !(aitr->range_end() > intersection.range_end() ) ) { ++aitr; }
                		//if( !(bitr->range_end() > intersection.range_end() ) ) { ++bitr; }
                		if( aitr->range_end() == intersection.range_end() ) { ++aitr; }
                		if( bitr->range_end() == intersection.range_end() ) { ++bitr; }
                	}
                	else // aitr is a range, but bitr is not
                	{
                		f(*bitr);

                		anb.emplace_back( std::move( value_type( (*bitr)(), ~bitr->get_bitfield() ) ) );

                		if( bitr->range_end() < aitr->range_end() )
                		{
                			anb.emplace_back( std::move( value_type( bitr->range_end() ).set_range_end( aitr->range_end() ) ) );
                		}

						// update positions
                 		if( ++bitr != benditr && !(aitr->range_end() > (*bitr)()) ) { ++aitr; }
                	}
                }
            	else
            	{
					if( bitr->is_range() ) // bitr is a range, but aitr is not
					{
						f(*aitr);

						bna.emplace_back( value_type( (*aitr)(), ~aitr->get_bitfield() ) );

                		if( aitr->range_end() < bitr->range_end() )
                		{
                			bna.emplace_back( std::move( value_type( aitr->range_end() ).set_range_end( bitr->range_end() ) ) );
                		}

						// update positions
						if( ++aitr != aenditr && !(bitr->range_end() > (*aitr)()) ) { ++bitr; }
					}
					else // neither aitr nor bitr are ranges and (*aitr)() == (*bitr)(); equivalent to if( !aitr->is_range() && !bitr->is_range() )
					{
						const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
						//const index_t bitfield = aitr->get_bitfield() & bitr->get_bitfield();
						if( bitfield ) // some overlap
						{
							f( value_type((*aitr)(),bitfield) );

							const auto bfa = aitr->get_bitfield() & ~bitfield;
							//index_t bfa = aitr->get_bitfield() & ~bitfield;
							if( bfa ) { anb.emplace_back( value_type( (*aitr)(), bfa ) ); }

							const auto bfb = bitr->get_bitfield() & ~bitfield;
							//index_t bfb = bitr->get_bitfield() & ~bitfield;
							if( bfb ) { bna.emplace_back( value_type( (*bitr)(), bfb ) ); }
						}
						else // no overlap
						{
							anb.push_back(*aitr);
							bna.push_back(*bitr);
						}

						// update positions
						++aitr;
						++bitr;
					}
            	}
            }
        }
    }

	// add anything that wasn't already processed above
	for( ; aitr != aenditr; ++aitr ) { anb.push_back(*aitr); }
	for( ; bitr != benditr; ++bitr ) { bna.push_back(*bitr); }

	anb.shrink_to_fit();
	bna.shrink_to_fit();

	return result;
}

// intersect and differences for aligned container; input is modified in-place
template< typename IndexT, typename UnaryFunction >
UnaryFunction inplace_set_differences_and_for_each_in_intersection(
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b,
		UnaryFunction f
	)
{
	using value_type = Hybrid_bitset_range_element<IndexT,true>;
	using std::cbegin; using std::cend;

	// nothing to gather?
	if( !a.has_overlap(b) ) { return f; } // this test should suffice; covers cases where either a or b, or both a and b, are empty

	auto& astor = a.storage();
	auto& bstor = b.storage();

	auto aitr = begin(astor);
	auto aenditr = end(astor);

	auto bitr = begin(bstor);
	auto benditr = end(bstor);

	auto anb = aitr; // a-not-b; this is the leading edge of the result, and anb <= aitr will hold
	auto bna = bitr; // b-not-a; this is the leading edge of the result, and bna <= bitr will hold

// /*
	// clear occlusion masks (we'll rebuild them while constructing sets a-not-b and b-not-a)
	a.clear_block_mask();
	b.clear_block_mask();
// */

	while( aitr != aenditr && bitr != benditr )
	{
        if( aitr->range_end() <= bitr->range_begin() ) { do { std::swap(*anb,*aitr); a.append_block_mask(*anb); ++anb; ++aitr; } while( aitr != aenditr && aitr->range_end() <= bitr->range_begin() ); }
        else
        {
            if( bitr->range_end() <= aitr->range_begin() ) { do { std::swap(*bna,*bitr); b.append_block_mask(*bna); ++bna; ++bitr; } while( bitr != benditr && bitr->range_end() <= aitr->range_begin() ); }
            else // *aitr and *bitr overlap in some way
            {
            	if( aitr->is_range() )
                {
                	if( bitr->is_range() ) // both aitr and bitr are ranges
                	{
                        // fast-forward to the end of shared range
                		const auto intersection( range_intersection( *aitr, *bitr ) ); // will always succeed, since overlap is guaranteed in this branch. There will by design not be overlap with anything already in storage.

                		f(intersection);

                    	if( pre_intersection_inplace_range_difference( a, b, aitr, bitr, intersection, anb, bna ) )
                		{
                			aenditr = end(astor); benditr = end(bstor); // need to refresh these, since inplace_range_difference has invalidated one of them
                		}

						// update positions
                		if( aitr->range_end() == intersection.range_end() ) { ++aitr; }
                		if( bitr->range_end() == intersection.range_end() ) { ++bitr; }
                	}
                	else // aitr is a range, but bitr is not
                	{
                		f(*bitr);

						// does the range element extend beyond the current location?
                		if( bitr->range_end() < aitr->range_end() )
                		{
                			// do we need to insert a new range element for the left-over range?
                			if( anb >= aitr ) {
                				anb = astor.emplace( anb ); // emplace default-constructed value_type
                				aitr = anb+1; aenditr = end(astor); // replace invalidated iterators
                			}
    						// set the left-over range
                			aitr->set_range( bitr->range_end(), aitr->range_end() );
                		}
                		else { ++aitr; } // the range was completely consumed

                		anb->set( *bitr, ~bitr->get_bitfield() );
						a.append_block_mask(*anb);

						// update positions
						++anb;
						++bitr; // bitr is not a range; always increment
                	}
                }
            	else
            	{
					if( bitr->is_range() ) // bitr is a range, but aitr is not
					{
						f(*aitr);

						// does the range element extend beyond the current location?
                		if( aitr->range_end() < bitr->range_end() )
                		{
                			// do we need to insert a new range element for the left-over range?
                			if( bna >= bitr ) {
                				//std::cout << "bna >= bitr d=" << std::distance(bna,bitr) << std::endl;
                				bna = bstor.emplace( bna ); // emplace default-constructed value_type
                				bitr = bna+1; benditr = end(bstor); // replace invalidated iterators
                			}
    						// set the left-over range
    						bitr->set_range( aitr->range_end(), bitr->range_end() );
                		}
                		else { ++bitr; } // the range was completely consumed

                		bna->set( *aitr, ~aitr->get_bitfield() );
						b.append_block_mask(*bna);

						// update positions
						++bna;
						++aitr; // aitr is not a range; always increment
					}
					else // neither aitr nor bitr are ranges and (*aitr)() == (*bitr)()
					{
						const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
						if( bitfield ) // some overlap
						{
							f( value_type( *aitr, bitfield ) );

							const auto bfa = aitr->get_bitfield() & ~bitfield;
							if( bfa ) { anb->set( *aitr, bfa ); a.append_block_mask(*anb); ++anb; }

							const auto bfb = bitr->get_bitfield() & ~bitfield;
							if( bfb ) { bna->set( *bitr, bfb ); b.append_block_mask(*bna); ++bna; }
						}
						else // no overlap
						{
							std::swap( *anb, *aitr ); a.append_block_mask(*anb); ++anb;
							std::swap( *bna, *bitr ); b.append_block_mask(*bna); ++bna;
						}

						// update positions
						++aitr;
						++bitr;
					}
            	}
            }
        }
    }

	// add anything that wasn't already processed above
	while( aitr != aenditr ) { std::swap(*anb,*aitr); a.append_block_mask(*anb); ++anb; ++aitr; }
	astor.erase( anb, aenditr );

	while( bitr != benditr ) { std::swap(*bna,*bitr); b.append_block_mask(*bna); ++bna; ++bitr; }
	bstor.erase( bna, benditr );

	//anb.shrink_to_fit();
	//bna.shrink_to_fit();

	return f;
}

// union
template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > set_union(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	using container_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	using std::cbegin; using std::cend;

	if( a.empty() ) { if( b.empty() ) { return container_t(); } else { return container_t(b); } }
	else if( b.empty() ) { return container_t(a); }

	auto aitr = cbegin(a.storage());
	const auto aenditr = cend(a.storage());

	auto bitr = cbegin(b.storage());
	const auto benditr = cend(b.storage());

	container_t theunion;
	auto& storage = theunion.m_storage;

	while( aitr != aenditr && bitr != benditr )
	{
        if( (*aitr)() < (*bitr)() )
        {
        	if( aitr->is_range() )
        	{
    			if( !theunion.empty() && storage.back().is_range() )
    			{
    				const auto extended = fuse_range( storage.back(), *aitr );
    				if( extended ) { storage.back() = std::move(extended); }
    				else { storage.emplace_back( *aitr ); }
    			}
				else { storage.emplace_back( *aitr ); }
        	}
            else if( !((*aitr)() < ( theunion.empty() ? 0 : storage.back().range_end() )) )
			{
				storage.emplace_back( *aitr );
			}
        	++aitr;
        }
        else
        {
            if( (*bitr)() < (*aitr)() )
        	{
            	if( bitr->is_range() )
            	{
        			if( !theunion.empty() && storage.back().is_range() )
        			{
        				const auto extended = fuse_range( storage.back(), *bitr );
        				if( extended ) { storage.back() = std::move(extended); }
        				else { storage.emplace_back( *bitr ); }
        			}
    				else { storage.emplace_back( *bitr ); }
            	}
                else if( !((*bitr)() < ( theunion.empty() ? 0 : storage.back().range_end() )) )
                {
    				storage.emplace_back( *bitr );
    			}
            	++bitr;
        	}
        	else // apos == bpos
        	{
            	if( aitr->is_range() )
            	{
            		if( bitr->is_range() )
            		{
            			const auto abunion = fuse_range( *aitr, *bitr ); // will always succeed, since range overlap is guaranteed
        				if( !theunion.empty() && storage.back().is_range() )
        				{
        					const auto extended = fuse_range( storage.back(), abunion );
        					if( extended ) { storage.back() = std::move(extended); }
        					else { storage.emplace_back( std::move(abunion) ); }
        				}
        				else { storage.emplace_back( std::move(abunion) ); }
            		}
            		else
            		{
            			if( !theunion.empty() && storage.back().is_range() )
            			{
             				const auto extended = fuse_range( storage.back(), *aitr );
            				if( extended ) { storage.back() = std::move(extended); }
            				else { storage.emplace_back( *aitr ); }
            			}
        				else { storage.emplace_back( *aitr ); }
            		}
            	}
            	else if( bitr->is_range() ) // bitr is range, but aitr is not
            	{
        			if( !theunion.empty() && storage.back().is_range() )
        			{
        				const auto extended = fuse_range( storage.back(), *bitr );
        				if( extended ) { storage.back() = std::move(extended); }
        				else { storage.emplace_back( *bitr ); }
        			}
    				else { storage.emplace_back( *bitr ); }
            	}
            	else // neither aitr nor bitr are ranges
            	{
            		storage.emplace_back( (*aitr)(), aitr->get_bitfield() | bitr->get_bitfield() );
            		if( storage.back().all_set() )
            		{
            			storage.back().set_range_end( aitr->range_end() );
            			if( storage.size() > 1 && (storage.end()-2)->is_range() ) // extend an existing range?
            			{
            				const auto extended = fuse_range( *(storage.end()-2), storage.back() );
            				if( extended ) { storage.pop_back(); storage.back() = std::move(extended); }
            			}
            		}
            	}
            	++aitr; ++bitr;
        	}
        }
    }

	// insert the rest
	auto itr = ( aitr != aenditr ? aitr : bitr );
	const auto end = ( aitr != aenditr ? aenditr : benditr );
	while( itr != end )
	{
		if( itr->is_range() )
		{
			if( !theunion.empty() && storage.back().is_range() )
			{
				const auto extended = fuse_range( storage.back(), *itr );
				if( extended ) { storage.back() = std::move(extended); }
				else { storage.emplace_back( *itr ); }
			}
			else { storage.emplace_back( *itr ); }
		}
		else
		{
			storage.emplace_back( *itr );
		}
		++itr;
	}

	// mask update for the set union operation; we can simply combine the block masks of each input container here
	theunion.update_block_mask(a.m_block_mask|b.m_block_mask);

	storage.shrink_to_fit();
	return theunion;
}

template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > operator|(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	return set_union(a,b);
}

template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > set_union(
		const std::vector< std::reference_wrapper< const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > > >& sets
	)
{
	using container_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	using std::cbegin; using std::cend;

	if( sets.size() == 1 ) { return container_t( sets.front() ); }

	typename container_t::block_mask_type block_mask(0);

	// concatenate input sequences
	const auto size_reserve = std::accumulate( cbegin(sets), cend(sets), 0, [](auto sum, const auto& a) { return sum+a.get().storagesize(); } );
	container_t theunion; auto& u = theunion.m_storage; u.reserve(size_reserve);

	for( const auto& src: sets ) { u.insert( u.end(), cbegin(src.get().storage()), cend(src.get().storage()) ); block_mask |= src.get().m_block_mask; }
	theunion.update_block_mask(block_mask);

	// sort
	std::sort( u.begin(), u.end() );

	// compact & clean-up
	auto pos = u.begin();
	auto fwd = pos+1;

	auto range_pos = u.end();
	while( fwd != u.end() )
	{
		if( (*pos)() < (*fwd)() )
		{
			if( pos->range_end() < (*fwd)() )
			{ // pos is an isolated range or a bitfield
				std::swap( *(++pos), *fwd ); ++fwd;
			}
			else if( pos->is_range() )
			{
				range_pos = pos;
				if( fwd->is_range() )
				{
					auto fused = apegrunt::fuse_range(*pos, *fwd);
					if(fused) { *pos = fused; ++fwd; }
					else { ++pos; std::swap(*pos,*fwd); ++fwd; }
				}
				else if( pos->range_end() > (*fwd)() )
				{
					++fwd;
				}
				else
				{
					std::swap( *(++pos), *fwd ); ++fwd;
				}
			}
			else
			{
				std::swap( *(++pos), *fwd ); ++fwd;
			}
		}
		else if( pos->is_range() ) // pos == fwd, but is it a range
		{
			range_pos = pos;
			if( fwd->is_range() )
			{
				auto fused = apegrunt::fuse_range(*pos, *fwd);
				if(fused) { *pos = fused; ++fwd; } // ranges overlap
				else { ++pos; std::swap( *pos, *fwd ); ++fwd; } // we've got two non-overlapping ranges; move fwd adjacent to pos, as there might be free slots (if not, then it's semantically a non-op)
			}
			else
			{
				if( pos->range_end() > (*fwd)() ) // range pos completely covers fwd
				{
					++fwd;
				}
				else // we're done with pos
				{
					std::swap( *(++pos), *fwd ); ++fwd;
				}
			}
		}
		else // pos == fwd, and pos is not a range, but fwd might be
		{
			if( fwd->is_range() )
			{
				*pos = *fwd;
				range_pos = pos;
			}
			else
			{
				pos->merge_bitfield( fwd->get_bitfield() );
				if( pos->all_set() )
				{
					if( range_pos != u.end() && range_pos->range_end() == (*pos)() )
					{
						range_pos->set_range_end( pos->range_end() ); pos = range_pos;
					}
					else
					{
						pos->set_range_end( pos->range_end() ); range_pos = pos;
					}
				}
			}
			++fwd;
		}
	}

	++pos;
	u.erase( pos, fwd ); // trim off excess elements
	u.shrink_to_fit(); // make it a tight fit

	return theunion;
}

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_HPP
