/** @file StateVector_impl_block_compressed_alignment_storage.hpp

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

#ifndef APEGRUNT_STATEVECTOR_IMPL_BLOCK_COMPRESSED_ALIGNMENT_STORAGE_HPP
#define APEGRUNT_STATEVECTOR_IMPL_BLOCK_COMPRESSED_ALIGNMENT_STORAGE_HPP

#include <cstdlib> // for std::div

//#include <iosfwd>
#include <vector>
#include <memory> // for std::make_shared
#include <algorithm> // for std::find and std::min

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/range/irange.hpp>

#include "Alignment_forward.h"

#include "StateVector.h"
#include "StateVector_impl_base.hpp"
#include "StateVector_state_types.hpp"
#include "StateVector_iterator.h"
#include "StateVector_iterator_impl_block_compressed_alignment_storage.hpp"

#include "StateVector_parser_forward.h"
//#include "StateVector_generator_forward.h"
#include "StateVector_mutator_forward.h"

#include "State_block.hpp"
#include "Alignment_impl_block_compressed_storage.hpp"
#include "aligned_allocator.hpp"

#include "Apegrunt_options.h"

#include "IndexVector.h"
#include "misc/Stopwatch.hpp"

namespace apegrunt {

template< typename StateT >
class StateVector_impl_block_compressed_alignment_storage : public StateVector_impl_base< StateVector_impl_block_compressed_alignment_storage<StateT>, StateT >
{
public:
	using state_t = StateT;
	using my_type = StateVector_impl_block_compressed_alignment_storage<state_t>;
	using base_type = StateVector_impl_base< my_type, state_t >;
	using block_type = typename StateVector<state_t>::block_type;

	enum { N=block_type::N };

	using parent_t = Alignment_impl_block_compressed_storage< my_type >;
	using block_index_t = typename base_type::block_index_t;
	using block_index_container_t = typename base_type::block_index_container_t;
	using block_index_container_ptr = typename base_type::block_index_container_ptr;

	using allocator_t = memory::AlignedAllocator<block_type,alignof(block_type)>;
	using block_storage_t = std::vector< std::vector< block_type, allocator_t > >;
	using block_storage_ptr = std::shared_ptr< block_storage_t >;

	//using block_adder_t = typename parent_t::block_adder_t;
	//using block_adder_ptr = typename parent_t::block_adder_ptr;
	using block_adder_t = Block_adder<state_t,char,32>;
	using block_adder_ptr = std::shared_ptr< block_adder_t >;

	using const_iterator = typename base_type::const_iterator;
	using iterator = typename base_type::iterator ;

	using value_type = typename base_type::value_type;

	//StateVector_impl_block_compressed_alignment_storage() = delete;
	StateVector_impl_block_compressed_alignment_storage() : base_type(), m_dirty(false) { }
	~StateVector_impl_block_compressed_alignment_storage() = default;

	StateVector_impl_block_compressed_alignment_storage( std::size_t id=0 )
	: base_type(id, std::to_string(id)),
	  m_block_storage(),
	  m_block_adder(), // always after m_block_storage
	  m_block_indices(),
	  m_cache(),
	  m_access_cache_block(),
	  m_access_cache_col(std::numeric_limits<std::size_t>::max()),
	  m_block_col(0),
	  m_pos(0),
	  m_size(0),
	  m_dirty(false),
	  m_has_block_indices(false),
	  m_cache_block_indices_mutex(),
	{
	}

	StateVector_impl_block_compressed_alignment_storage( const my_type& other )
	: base_type( other ),
	  m_block_storage(other.m_block_storage),
	  m_block_adder(other.m_block_adder), // always after m_block_storage
	  m_block_indices(other.m_block_indices),
	  m_cache(other.m_cache),
	  m_access_cache_block(other.m_access_cache_block),
	  m_access_cache_col(other.m_access_cache_col),
	  m_block_col(other.m_block_col),
	  m_pos(other.m_pos),
	  m_size(other.m_size),
	  m_dirty(other.m_dirty),
	  m_has_block_indices(other.m_has_block_indices),
	  m_cache_block_indices_mutex(),
	{
	}

	StateVector_impl_block_compressed_alignment_storage( my_type&& other ) noexcept
	: base_type( other ),
	  m_block_storage( std::move(other.m_block_storage) ),
	  m_block_adder( std::move(other.m_block_adder) ), // always after m_block_storage
	  m_block_indices( std::move(other.m_block_indices) ),
	  m_cache( std::move(other.m_cache) ),
	  m_access_cache_block( std::move(other.m_access_cache_block) ),
	  m_access_cache_col(other.m_access_cache_col),
	  m_block_col(other.m_block_col),
	  m_pos(other.m_pos),
	  m_size(other.m_size),
	  m_dirty(other.m_dirty),
	  m_has_block_indices(other.m_has_block_indices),
	  m_cache_block_indices_mutex(),
	{
	}

	my_type& operator=( my_type&& other ) noexcept
	{
		// in base class:
		this->set_id( other.id() );
		this->set_multiplicity( other.multiplicity() );
		this->set_weight( other.weight() );
		this->set_id_string( other.id_string() );

		// in my_type:
		m_block_storage = std::move(other.m_block_storage);
		m_block_indices = std::move(other.m_block_indices);
		m_cache = other.m_cache;
		m_access_cache_block = std::move(other.m_access_cache_block);
		m_access_cache_col = other.m_access_cache_col;
		m_block_col = other.m_block_col;
		m_pos = other.m_pos;
		m_size = other.m_size;
		m_dirty = other.m_dirty;
		m_has_block_indices = other.m_has_block_indices;
		return *this;
	}
// */
// /*
	my_type& operator=( const my_type& other )
	{
		// in base class:
		this->set_id( other.id() );
		this->set_multiplicity( other.multiplicity() );
		this->set_weight( other.weight() );
		this->set_id_string( other.id_string() );

		// in my_type:
		m_block_storage = other.m_block_storage;
		m_block_indices = other.m_block_indices;
		m_cache = other.m_cache;
		m_access_cache_block = other.m_access_cache_block;
		m_access_cache_col = other.m_access_cache_col;
		m_block_col = other.m_block_col;
		m_pos = other.m_pos;
		m_size = other.m_size;
		m_dirty = other.m_dirty;
		m_has_block_indices = other.m_has_block_indices;
		return *this;
	}

	StateVector_impl_block_compressed_alignment_storage( block_adder_ptr& block_adder, std::size_t id, const std::string& id_string="", std::size_t size_hint=0 )
	: base_type( id, id_string ),
	  m_block_storage( block_adder->get_block_storage() ),
	  m_block_adder( block_adder ),
	  m_block_indices(),
	  m_cache(),
	  m_access_cache_block(),
	  m_access_cache_col(0),
	  m_block_col(0),
	  m_pos(0),
	  m_size(0),
	  m_dirty(false),
	  m_has_block_indices(false),
	  m_cache_block_indices_mutex(),
	{
		//if( size_hint > 0 ) { m_block_indices->reserve(size_hint); }
	}

	inline const_iterator cbegin() const { return const_iterator( std::make_unique<const_iterator_impl>( 0, 0, this ) ); }
    //inline const_iterator cend() const { return const_iterator( std::make_unique<const_iterator_impl>( m_pos, (m_pos==0 ? m_block_indices->size() : m_block_indices->size()-1 ), this ) ); }
    inline const_iterator cend() const { return std::make_unique<const_iterator_impl>( apegrunt::get_pos_in_block(m_size), apegrunt::get_last_block_index(m_size)+(apegrunt::get_pos_in_block(m_size)==0), this ); }

    inline iterator begin() { return iterator( std::make_unique<iterator_impl>( 0, 0, this ) ); }
    //inline iterator end() { return iterator( std::make_unique<iterator_impl>( m_pos, (m_pos==0 ? m_block_indices->size() : m_block_indices->size()-1 ), this ) ); }
    inline iterator end() { return std::make_unique<iterator_impl>( apegrunt::get_pos_in_block(m_size), apegrunt::get_last_block_index(m_size)+(apegrunt::get_pos_in_block(m_size)==0), this ); }

    inline value_type operator[]( std::size_t index ) const
    {
    	const auto a = std::div( index, N );

    	if( std::size_t(a.quot) != m_access_cache_col ) // cast to std::size_t to rid compiler nag about different signedness comparison; a.quot is always positive
    	{
    		this->get_block(a.quot);
    		//m_access_cache_col = a.quot;
    		//m_access_cache_block = this->get_block(m_access_cache_col); // get_block will populate the cache
    	}
   		return m_access_cache_block[a.rem];

    	//return this->get_block(a.quot)[a.rem];
    }

    bool operator==( const my_type& rhs ) const
    {
    	using boost::get;

    	if( this->size() == rhs.size() )
    	{
    		if( m_block_storage == rhs.m_block_storage ) // this and rhs belong to the same Alignment
    		{
    			for( auto index_pair: zip_range( *m_block_indices, *rhs.m_block_indices ) )
    			{
    				if( get<0>(index_pair) != get<1>(index_pair) ) { return false; } // compare pair indices
    			}
    		}
    		else
    		{
				for( std::size_t i=0; i < m_block_indices->size(); ++i )
				{
					if( this->get_block(i) != rhs.get_block(i) ) { return false; }
				}
    		}
    	}
    	else
    	{
    		return false;
    	}
    	return true;
    }
    inline block_type get_block( std::size_t index ) const
    {
    	// simplified version
		if( index != m_access_cache_col )
		{
			m_access_cache_col = index;
			m_access_cache_block = m_block_adder->get(index, (*this->get_block_indices())[index]); // Block_adder will return block_type() if block index is not found.
		}
		return m_access_cache_block;
    }

    inline const block_index_container_ptr get_block_indices() const
    {
		if( !m_has_block_indices ) // no unnecessary fighting for the mutex
		{
			m_cache_block_indices_mutex.lock();
			if( !m_block_indices ) { this->cache_block_indices(); }
			m_cache_block_indices_mutex.unlock();
		}
    	return m_block_indices;
    }

	inline std::size_t operator&&( const my_type& rhs ) const
	{
		std::size_t n = 0;

		const auto& lhs_block_indices = *(this->get_block_indices());
		const auto& rhs_block_indices = *(rhs.get_block_indices());

		const std::size_t n_indices = std::min( lhs_block_indices.size(), rhs_block_indices.size() ) -1; // "-1": process the last, potentially partially filled block separately
		if( m_block_storage == rhs.m_block_storage ) // this and rhs belong to the same Alignment
		{
			const auto& block_storage = *(m_block_storage.get());
			for( std::size_t i = 0; i < n_indices; ++i )
			{
				// compare pair indices first; count identical states only if needed
				n += lhs_block_indices[i] == rhs_block_indices[i] ? N : count_identical( block_storage[i][lhs_block_indices[i]], block_storage[i][rhs_block_indices[i]] );
			}
		}
		else
		{
			for( std::size_t i = 0; i < n_indices; ++i ) // compare pair indices
			{
				n += count_identical( this->get_block(i), rhs.get_block(i) );
			}
		}
		{
			const auto my_block = this->get_block(n_indices); // last common block
			const auto rhs_block = rhs.get_block(n_indices); // last common block
			for( std::size_t i=0; i< std::min(m_pos,rhs.m_pos); ++i ) { my_block[i] == rhs_block[i] && ++n; }
		}
		return n;
	}

	inline bool is_similar_to( const my_type& rhs, std::size_t min_identical ) const
	{
	inline bool is_similar_to( const my_type& rhs, std::size_t min_identical ) const
	{
		using boost::get;
		std::size_t n(0);

		const auto lhs_indices = this->get_block_indices(); // hold 'em so we don't lose 'em
		const auto rhs_indices = rhs.get_block_indices(); // hold 'em so we don't lose 'em

		const auto& lhs_block_indices = *lhs_indices;
		const auto& rhs_block_indices = *rhs_indices;

		const auto size = std::min(this->size(),rhs.size());
		const auto lbs = apegrunt::get_last_block_size(size);
		const auto lbi = apegrunt::get_last_block_index(size);
		const std::size_t n_indices = std::min( lhs_block_indices.size(), rhs_block_indices.size() ) - (lbs!=apegrunt::StateBlock_size); // "-1": process the last, potentially partially filled block separately

		if( m_block_storage == rhs.m_block_storage ) // this and rhs belong to the same Alignment
		{
			auto nidentical = [](const auto& blocks, auto i, auto j) { return apegrunt::count_identical( blocks[i], blocks[j] ); };
			auto remain(n_indices);
			for( auto zipped: apegrunt::zip_range(lhs_block_indices,rhs_block_indices,*m_block_storage) )
			{
				n += ( get<0>(zipped) == get<1>(zipped) ? N : nidentical(get<2>(zipped), get<0>(zipped), get<1>(zipped)) );
				if( min_identical > n+N*(--remain)+lbs ) { return false; }
				else if( n >= min_identical ) { return true; }
			}
		}
		else
		{
			for( std::size_t i = 0; i < n_indices; ++i )
			{
				n += apegrunt::count_identical( this->get_block(i), rhs.get_block(i) );
				if( min_identical > n+(n_indices-i)*N+lbs ) { return false; } // +1 for the last (it could be partially filled, but we give it the benefit of doubt)
				else if ( n >= min_identical ) { return true; }
			}
		}
		// we will seldom actually reach this point in practice
		{
			const auto my_block = this->get_block(lbi); // last common block
			const auto rhs_block = rhs.get_block(lbi); // last common block
			for( std::size_t i=0; i< std::min(m_pos,rhs.m_pos); ++i ) { my_block[i] == rhs_block[i] && ++n; }
		}
		return n >= min_identical; // ? true : false;
	}

	inline std::size_t size() const { return m_size; }

	inline std::size_t bytesize() const { return apegrunt::bytesize(*m_block_indices); }

	inline const std::type_info& type() const { return typeid(my_type); }

	block_storage_ptr get_block_storage() const
	{
		return m_block_storage;
	}

private:
	//using iterator_impl = apegrunt::iterator::StateVector_iterator_impl_block_compressed_alignment_storage< State_holder<state_t>, StateVector_impl_block_compressed_alignment_storage<state_t> >;
	//using const_iterator_impl = apegrunt::iterator::StateVector_iterator_impl_block_compressed_alignment_storage< State_holder<state_t>, StateVector_impl_block_compressed_alignment_storage<state_t> >;

	using iterator_impl = apegrunt::iterator::StateVector_iterator_impl_block_compressed_alignment_storage< state_t, StateVector_impl_block_compressed_alignment_storage<state_t> >;
	using const_iterator_impl = apegrunt::iterator::StateVector_iterator_impl_block_compressed_alignment_storage< state_t, StateVector_impl_block_compressed_alignment_storage<state_t> >;

	friend iterator_impl;
	friend const_iterator_impl;

	template< typename IntegerT > static constexpr IntegerT popcnt( IntegerT mask ) { IntegerT n(0); while( mask ) { ++n; mask >>= 1; } return n; }

	block_storage_ptr m_block_storage;
	block_adder_ptr m_block_adder;
	mutable block_index_container_ptr m_block_indices;
	block_type m_cache;
	mutable block_type m_access_cache_block;
	mutable std::size_t m_access_cache_col;
	std::array< std::size_t, number_of_states<state_t>::value > m_frequencies;
	std::size_t m_block_col;
	std::size_t m_pos;
	std::size_t m_size;
	bool m_dirty;
	//stopwatch::stopwatch m_build_timer_parse; // ( Apegrunt_options::verbose() ? Apegrunt_options::get_out_stream() : nullptr ); // for timing statistics
	//stopwatch::stopwatch m_build_timer_store;

	mutable bool m_has_block_indices;
	mutable std::mutex m_cache_block_indices_mutex;


	enum : std::size_t { MODULO_MASK=N-1 };
	enum { MODULO_NBITS=my_type::popcnt<std::size_t>(MODULO_MASK) };

	// Parser interface
	template< typename T >
	inline void push_back( const T& state )
	//inline void push_back( const char& state )
	{
		m_cache[m_pos] = state; ++m_pos; m_dirty=true; ++m_size;
		if( N == m_pos ) { this->flush_block_buffer(); }
	}

	template< typename T >
	inline void push_backN( const T *state )
	//inline void push_back( const char& state )
	{
		for( auto i=0; i < N; ++i )
		{
			m_cache[i] = *(state+i);
		}
		m_dirty=true; m_pos=N; // m_size+=N;
		this->flush_block_buffer();
	}

	inline bool append( const std::string& state_string )
	{
		m_block_indices.reserve( m_block_indices.size() + state_string.size()/N + (state_string.size() % N == 0 ? 0 : 1) );
		for( const auto& state: state_string ) { this->push_back( state ); }
		this->flush_block_buffer();
		return true;
	}

	inline bool appendN( const std::string& state_string )
	{
		const auto Nfull_blocks = state_string.size() / N;
		const auto Nreminder_chars = state_string.size() & MODULO_MASK;

		m_block_indices.reserve( m_block_indices.size() + Nfull_blocks + (Nreminder_chars ? 1 : 0) );
		for( auto i=0; i < Nfull_blocks; ++i ) { this->push_backN( state_string.data()+i*N ); }
		m_size += Nfull_blocks*N;
		if( Nreminder_chars )
		{
			for( auto i=0; i < Nreminder_chars; ++i ) { this->push_back( *(state_string.data()+Nfull_blocks*N+i) ); }
			this->flush_block_buffer();
		}
		return true;
	}

	void assign( const std::string& state_string )
	{
		this->clear();
		this->append(state_string); /* this->flush_block_buffer(); */ // don't flush twice, when the implementation uses append()
	}

	//> clear internal state, with the exception of link to parent Alignment
	void clear()
	{
		m_block_indices->clear();
		m_cache.clear();
		m_frequencies.fill(0);
		m_block_col=0;
		m_pos=0;
		m_dirty=false;
		m_size=0;
	}

	inline void flush_block_buffer()
	{
		if(m_dirty)
		{
/* Original/reference code
			if( m_block_storage->size() < m_block_col+1 ) { m_block_storage->emplace_back( 1, m_cache ); m_block_indices.push_back(0); }
			else
			{
				using std::cbegin; using std::cend;
				auto& block_list = (*m_block_storage)[m_block_col];
				auto list_pos = std::find( cbegin(block_list), cend(block_list), m_cache );
				//auto list_pos = block_list.find( m_cache );
				if( list_pos == cend(block_list) )
				{
					block_list.emplace_back( m_cache ); m_block_indices.push_back( block_list.size()-1 );
				}
				else
				{
					m_block_indices.push_back( std::distance( cbegin(block_list), list_pos ) );
				}
			}
*/
			// new TST-based routine
			m_block_indices.push_back( m_block_adder->insert( m_cache, m_block_col, m_index ) );

			if( N  == m_pos ) // test if cache is full or partially filled
			{
				// update access cache
				m_access_cache_block = m_cache;
				m_access_cache_col = m_block_col;
				// reset construction cache and update column position
				m_cache.clear(); m_pos=0; ++m_block_col;
			}
			m_dirty=false;
		}
	}

	inline void cache_block_indices() const
	{
		if( !m_block_indices ) { m_block_indices = std::make_shared<block_index_container_t>(); }
		auto& bis = *m_block_indices;
		bis.reserve( m_block_storage->size() );

		auto block_accounting = m_block_adder->get_block_accounting();

		const block_index_t my_id(this->id());

		for( const auto& column_acc: *block_accounting )
		{
			auto block_index(0);
			for( const auto& acc: column_acc )
			{
				if( contains(acc, my_id) ) { bis.push_back(block_index); break; }
				++block_index;
			}
		}

		m_has_block_indices = true;
	}

	// Allow parser access to private members
	STATEVECTOR_PARSER_GRAMMAR_FRIENDS(my_type) // macro defined in #include "StateVector_parser_forward.h"

	friend class StateVector_mutator<my_type>;

	/// boost.serialization interface.
	friend class boost::serialization::access;
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(boost::serialization::base_object< base_type >(*this));
        ar & BOOST_SERIALIZATION_NVP(m_block_storage);
        ar & BOOST_SERIALIZATION_NVP(m_block_adder);
        ar & BOOST_SERIALIZATION_NVP(m_block_indices);
        ar & BOOST_SERIALIZATION_NVP(m_cache);
        //ar & BOOST_SERIALIZATION_NVP(m_frequencies);
        ar & BOOST_SERIALIZATION_NVP(m_pos);
        ar & BOOST_SERIALIZATION_NVP(m_dirty);
    }

	/// Helper class for memory management thru std::shared_ptr
	class deleter
	{
	public:
		void operator()( my_type* p )
		{
			delete p;
		}
	};
	friend class deleter;

};

template< typename StateT >
class StateVector_mutator< StateVector_impl_block_compressed_alignment_storage< StateT > >
{
public:
	using statevector_t = StateVector_impl_block_compressed_alignment_storage< StateT >;
	StateVector_mutator( statevector_t *statevector ) : m_statevector(statevector) { }
	~StateVector_mutator() { m_statevector->flush_block_buffer(); } // flush when done

	template< typename StateT2 >
	void operator()( StateT2 state ) { m_statevector->push_back( state ); }

	template< typename RealT >
	inline void set_weight( RealT weight ) { m_statevector->set_weight( weight ); }

private:
	statevector_t *m_statevector;
};

/*
template< typename StateT >
class StateVector_subscript_proxy
{
public:
	StateVector_subscript_proxy() { }
	~StateVector_subscript_proxy() { }

	StateVector_subscript_proxy( const std::vector<StateT>* container ) : m_container(container) { }

	inline StateT operator[]( std::size_t index ) const { return (*m_container)[index]; }

private:
	const std::vector<StateT> *m_container;
};
*/
} // namespace apegrunt

#endif // APEGRUNT_STATEVECTOR_IMPL_BLOCK_COMPRESSED_ALIGNMENT_STORAGE_HPP

