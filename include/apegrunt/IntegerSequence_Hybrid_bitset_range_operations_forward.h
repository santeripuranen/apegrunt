/** @file IntegerSequence_Hybrid_bitset_range_operations_forward.h

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

#ifndef APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_FORWARD_H
#define APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_FORWARD_H

#include <ostream>
#include <vector>

#include "IntegerSequence_forward.h"
#include "IntegerSequence_operations_forward.h"
#include "IntegerSequence_Hybrid_bitset_range_forward.h"

namespace apegrunt {

template< typename IndexT, bool Aligned >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& container );

template< typename IndexT >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& container );

template< typename IndexT, bool Aligned >
IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > set_intersection( const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& a, const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& b );

template< typename IndexT, typename RealT, bool Aligned=true >
RealT intersect_and_gather( const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& a, const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& b, const std::vector<RealT>& weights );

template< typename IndexT, typename RealT >
RealT intersect_and_gather( const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a, const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b, const std::vector<RealT>& weights );

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_FORWARD_H
