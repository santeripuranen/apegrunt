#include "apegrunt/StateVector_impl_block_compressed_alignment_storage.hpp"

namespace apegrunt {

struct SSE2_tag;

template<>
bool is_similar_to_impl_sse2<nucleic_acid_state_t>(
	const StateVector_impl_block_compressed_alignment_storage< nucleic_acid_state_t >& lhs,
	const StateVector_impl_block_compressed_alignment_storage< nucleic_acid_state_t >& rhs,
	std::size_t min_identical )
{
	return is_similar_to_impl<nucleic_acid_state_t,SSE2_tag>(lhs,rhs,min_identical);
}

} // namespace apegrunt
