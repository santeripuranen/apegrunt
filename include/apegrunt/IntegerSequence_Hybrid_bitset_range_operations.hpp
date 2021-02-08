/** @file IntegerSequence_Hybrid_bitset_range_operations.hpp
 
	Copyright (c) 2018-2020 Santeri Puranen.

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

#include "IntegerSequence_forward.h"
#include "IntegerSequence_iterator.hpp"
#include "IntegerSequence_operations_forward.h"

#include "misc/Array_view.h"

namespace apegrunt {

// This version is aligned at value_type::BITSTRING_SIZE chunks and does not
// have single index entries, but only bitstrings and interval entries.
template< typename IndexT >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& container )
{
	os << std::hex; // use hex throughout

	bool first = true;
	for( const auto& pos: container.m_storage )	{ if( first ) { first = false; os << pos; } else { os << "," << pos; } }

	os << std::dec;

	return os;
}

template< typename IndexT >
Hybrid_bitset_range_element<IndexT,true> fuse_range(
		const Hybrid_bitset_range_element<IndexT,true>& a,
		const Hybrid_bitset_range_element<IndexT,true>& b
	)
{
	using element_t = Hybrid_bitset_range_element<IndexT,true>;

	//element_t fused;
	// check overlap
	if( a.range_end() < b() || b.range_end() < a() ) { return element_t(); } // elements do not overlap

	return element_t().set_range( std::min(a(),b()), std::max(a.range_end(),b.range_end()) );
}

template< typename IndexT >
Hybrid_bitset_range_element<IndexT,true> intersect_range(
		const Hybrid_bitset_range_element<IndexT,true>& a,
		const Hybrid_bitset_range_element<IndexT,true>& b
	)
{
	using element_t = Hybrid_bitset_range_element<IndexT,true>;

	// check overlap
	//if( a.range_end() < b() || b.range_end() < a() ) { return element_t(); } // elements do not overlap

	return element_t().set_range( std::max(a(),b()), std::min(a.range_end(),b.range_end()) );
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

// intersection
template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > set_intersection(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	//std::cout << "apegrunt::set_intersection()"; std::cout.flush();
	using container_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	using value_type = Hybrid_bitset_range_element<IndexT,true>;
	using index_t = typename value_type::value_type;

	container_t intersection;

	// nothing to intersect?
	//if( a.is_empty() || b.is_empty() ) { return intersection; }
	if( !a.has_overlap(b) ) { return intersection; } // this test should suffice; covers cases where either a or b, or both a and b, are empty

	auto& storage( intersection.m_storage );

	auto aitr = cbegin(a.m_storage);
	const auto aenditr = cend(a.m_storage);

	auto bitr = cbegin(b.m_storage);
	const auto benditr = cend(b.m_storage);

	// premature (and potentially excessive) memory allocation is costly and slows us down
	//intersection.m_storage.reserve(std::min(a.m_storage.size(),b.m_storage.size())); // can't need more than this, but could be way more than we need

	while( aitr != aenditr && bitr != benditr )
	{
        if( !(aitr->range_end() > (*bitr)()) ) { ++aitr; }
        else
        {
        	if( !(bitr->range_end() > (*aitr)()) ) { ++bitr; }
            else // *aitr and *bitr overlap in some way
            {
            	if( aitr->is_range() )
                {
            		if( bitr->is_range() ) // both aitr and bitr are ranges
                	{
            			// fast-forward to the end of shared range
						storage.emplace_back( intersect_range( *aitr, *bitr ) ); // will away be succesful, since overlap is guaranteed in this branch.
						/* There will by design not be overlap with anything already in storage, so we can skip the following conditions
						const auto ab = intersect_range( *aitr, *bitr ); // will away be succesful, since overlap is guaranteed in this branch
						if( !storage.empty() && storage.back().is_range() )
						{
							// can we extend and existing range with ab?
							const auto extended = fuse_range( storage.back(), ab );
							if( extended ) { storage.back() = std::move(extended); }
							else { storage.emplace_back( ab ); }
						}
						else { storage.emplace_back( ab ); }
						*/
						// update block mask
						intersection.update_block_mask( container_t::generate_block_mask( storage.back() ) );
						// update positions
						//if( !(aitr->range_end() > ab.range_end() ) ) { ++aitr; }
						//if( !(bitr->range_end() > ab.range_end() ) ) { ++bitr; }
						if( !(aitr->range_end() > storage.back().range_end() ) ) { ++aitr; }
						if( !(bitr->range_end() > storage.back().range_end() ) ) { ++bitr; }
                 	}
                	else // aitr is a range, but bitr is not
                	{
						// extend intersection
						//storage.emplace_back( (*bitr)(), bitr->get_bitfield() );
						storage.emplace_back( *bitr );
						intersection.update_block_mask( container_t::generate_block_mask( storage.back() ) );
						// update positions
						++bitr;
						if( !(aitr->range_end() > (*bitr)()) ) { ++aitr; }
						//if( !(aitr->range_end() > bitr->range_end()) ) { ++aitr; }
                	}
                }
            	else
            	{
					if( bitr->is_range() ) // bitr is a range, but aitr is not
					{
						// extend intersection
						//storage.emplace_back( (*aitr)(), aitr->get_bitfield() );
						storage.emplace_back( *aitr );
						intersection.update_block_mask( container_t::generate_block_mask( storage.back() ) );
						// update positions
						++aitr;
						if( !(bitr->range_end() > (*aitr)()) ) { ++bitr; }
						//if( !(bitr->range_end() > aitr->range_end()) ) { ++bitr; }
					}
					else // neither aitr nor bitr are ranges, so (*aitr)() == (*bitr)() must hold; equivalent to if( !aitr->is_range() && !bitr->is_range() )
					{
						const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
						if( bitfield )
						{
							storage.emplace_back( (*aitr)(), bitfield );
							intersection.update_block_mask( container_t::generate_block_mask( storage.back() ) );
						}

						// update positions
						++aitr;
						++bitr;
					}
            	}
            }
        }
    }
	storage.shrink_to_fit();
	return intersection;
}

template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > operator&(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	return set_intersection(a,b);
}

// intersect and gather for aligned container
template< typename IndexT, typename RealT >
RealT intersect_and_gather(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b,
		const std::vector<RealT>& weights
	)
{
	using value_type = Hybrid_bitset_range_element<IndexT,true>;
	using index_t = typename value_type::value_type;
	using real_t = RealT;

	using std::cbegin; using std::cend;

	real_t weight(0);

	// nothing to gather?
	//if( a.is_empty() || b.is_empty() ) { return weight; }
	if( !a.has_overlap(b) ) { return weight; } // this test should suffice; covers cases where either a or b, or both a and b, are empty

	auto aitr = cbegin(a.m_storage);
	const auto aenditr = cend(a.m_storage);

	auto bitr = cbegin(b.m_storage);
	const auto benditr = cend(b.m_storage);

	while( aitr != aenditr && bitr != benditr )
	{
        if( !(aitr->range_end() > (*bitr)()) ) { ++aitr; }
        else
        {
            if( !(bitr->range_end() > (*aitr)()) ) { ++bitr; }
            else // *aitr and *bitr overlap in some way
            {
            	if( aitr->is_range() )
                {
                	if( bitr->is_range() ) // both aitr and bitr are ranges
                	{
                        // fast-forward to the end of shared range
                		const auto range( intersect_range( *aitr, *bitr ) ); // will away be succesful, since overlap is guaranteed in this branch. There will by design not be overlap with anything already in storage.
                		weight += value_type::gather( weights.data(), range );

                		if( !(aitr->range_end() > range.range_end() ) ) { ++aitr; }
                		if( !(bitr->range_end() > range.range_end() ) ) { ++bitr; }
                 	}
                	else // aitr is a range, but bitr is not
                	{
                		weight += value_type::gather( &weights[bitr->range_begin()], bitr->get_bitfield() );

                		// update positions
                		++bitr;
                		if( !(aitr->range_end() > (*bitr)()) ) { ++aitr; }
                	}
                }
            	else
            	{
					if( bitr->is_range() ) // bitr is a range, but aitr is not
					{
						weight += value_type::gather( &weights[aitr->range_begin()], aitr->get_bitfield() );

						// update positions
						++aitr;
						if( !(bitr->range_end() > (*aitr)()) ) { ++bitr; }
					}
					else // neither aitr nor bitr are ranges, so (*aitr)() == (*bitr)() must hold; equivalent to if( !aitr->is_range() && !bitr->is_range() )
					{
						const auto bitfield = aitr->get_bitfield() & bitr->get_bitfield();
						if( bitfield )
						{
							weight += value_type::gather( &weights[aitr->range_begin()], bitfield );
						}

						// update positions
						++aitr;
						++bitr;
					}
            	}
            }
        }
    }
	return weight;
}

// union
template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > set_union(
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a,
		const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b
	)
{
	//std::cout << "apegrunt::set_union(Hybrid_bitset_range_element): a.size()=" << a.size() << " a.m_storage.size()=" << a.m_storage.size() << " | b.size()=" << b.size() << " b.m_storage.size()=" << b.m_storage.size() << std::endl;
	using container_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	using range_type = typename container_t::range_t;
	//using value_type = typename range_type::value_type;

	using std::cbegin; using std::cend;

	if( a.is_empty() ) { if( b.is_empty() ) { return container_t(); } else { return container_t(b); } }
	else if( b.is_empty() ) { return container_t(a); }

	auto aitr = cbegin(a.m_storage);
	const auto aenditr = cend(a.m_storage);

	auto bitr = cbegin(b.m_storage);
	const auto benditr = cend(b.m_storage);

	container_t theunion;
	auto& storage = theunion.m_storage;

	while( aitr != aenditr && bitr != benditr )
	{
        if( (*aitr)() < (*bitr)() )
        {
        	if( aitr->is_range() )
        	{
    			if( !theunion.is_empty() && storage.back().is_range() )
    			{
    				const auto extended = fuse_range( storage.back(), *aitr );
    				if( extended ) { storage.back() = std::move(extended); }
    				else { storage.emplace_back( *aitr ); }
    			}
				else { storage.emplace_back( *aitr ); }
				//theunion.update_block_mask( container_t::generate_block_mask( theunion.m_storage.back() ) );
        	}
            else if( !((*aitr)() < ( theunion.is_empty() ? 0 : storage.back().range_end() )) )
			{
				storage.emplace_back( *aitr );
				//theunion.update_block_mask( container_t::generate_block_mask( storage.back() ) );
			}
        	++aitr;
        }
        else
        {
        	//if( bpos < apos )
            if( (*bitr)() < (*aitr)() )
        	{
            	if( bitr->is_range() )
            	{
        			if( !theunion.is_empty() && storage.back().is_range() )
        			{
        				const auto extended = fuse_range( storage.back(), *bitr );
        				if( extended ) { storage.back() = std::move(extended); }
        				else { storage.emplace_back( *bitr ); }
        			}
    				else { storage.emplace_back( *bitr ); }
    				//theunion.update_block_mask( container_t::generate_block_mask( storage.back() ) );
            	}
                else if( !((*bitr)() < ( theunion.is_empty() ? 0 : storage.back().range_end() )) )
                {
    				storage.emplace_back( *bitr );
    				//theunion.update_block_mask( container_t::generate_block_mask( storage.back() ) );
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
        				if( !theunion.is_empty() && storage.back().is_range() )
        				{
        					const auto extended = fuse_range( storage.back(), abunion );
        					if( extended ) { storage.back() = std::move(extended); }
        					else { storage.emplace_back( std::move(abunion) ); }
        				}
        				else { storage.emplace_back( std::move(abunion) ); }
            		}
            		else
            		{
            			if( !theunion.is_empty() && storage.back().is_range() )
            			{
             				const auto extended = fuse_range( storage.back(), *aitr );
            				if( extended ) { storage.back() = std::move(extended); }
            				else { storage.emplace_back( *aitr ); }
            			}
        				else { storage.emplace_back( *aitr ); }
            		}
    				//theunion.update_block_mask( container_t::generate_block_mask( storage.back()) );
            	}
            	else if( bitr->is_range() ) // bitr is range, but aitr is not
            	{
        			if( !theunion.is_empty() && storage.back().is_range() )
        			{
        				const auto extended = fuse_range( storage.back(), *bitr );
        				if( extended ) { storage.back() = std::move(extended); }
        				else { storage.emplace_back( *bitr ); }
        			}
    				else { storage.emplace_back( *bitr ); }
    				//theunion.update_block_mask( container_t::generate_block_mask( storage.back()) );
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
    				//theunion.update_block_mask( container_t::generate_block_mask((*aitr)()) );
            	}
            	++aitr; ++bitr;
        	}
        }
		//theunion.update_block_mask( container_t::generate_block_mask( storage.back() ) );
    }
	//std::cout << "apegrunt::set_union(): finish it off" << std::endl;
	// insert the rest
	auto itr = ( aitr != aenditr ? aitr : bitr );
	const auto end = ( aitr != aenditr ? aenditr : benditr );
	while( itr != end )
	{
		if( itr->is_range() )
		{
			if( !theunion.is_empty() && storage.back().is_range() )
			{
				const auto extended = fuse_range( storage.back(), *itr );
				if( extended ) { storage.back() = std::move(extended); }
				else { storage.emplace_back( *itr ); }
			}
			else { storage.emplace_back( *itr ); }
			//theunion.update_block_mask( container_t::generate_block_mask( storage.back()) );
		}
		else
		{
			storage.emplace_back( *itr );
			//theunion.update_block_mask( container_t::generate_block_mask( storage.back() ) );
		}
		//theunion.update_block_mask( container_t::generate_block_mask( storage.back() ) );
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

// it's cheaper to perform a sequence of union operations when the individual ops are deferred and performed all in one go
template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > set_union(
		const std::vector< std::reference_wrapper< const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > > >& sets
	)
{
	using container_t = IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >;
	using range_type = typename container_t::range_t;
	//using value_type = typename range_type::value_type;

	using std::cbegin; using std::cend;

	if( sets.size() == 1 ) { return container_t( sets.front() ); }

	typename container_t::block_mask_type block_mask(0);

	// calculate required storage size and reserve space
	const auto size_reserve = std::accumulate( cbegin(sets), cend(sets), 0, [](auto sum, const auto& a) { return sum+a.get().storagesize(); } );
	container_t theunion; auto& u = theunion.m_storage; u.reserve(size_reserve);

	// concatenate input sequences
	for( const auto& src: sets ) { u.insert( u.end(), cbegin(src.get().m_storage), cend(src.get().m_storage) ); block_mask = block_mask | src.get().m_block_mask; }
	theunion.update_block_mask(block_mask);

	// sort elements
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
