/** @file IntegerSequence_operations_forward.h

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

#ifndef APEGRUNT_INTEGERSEQUENCE_OPERATIONS_FORWARD_H
#define APEGRUNT_INTEGERSEQUENCE_OPERATIONS_FORWARD_H

#include <vector>
#include <ostream>

#include "IntegerSequence_forward.h"

namespace apegrunt {

// forward declarations

template< typename SequenceT >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< SequenceT >& container );

template< typename ContainerT >
ContainerT set_intersection( const ContainerT& a, const ContainerT& b );

template< typename SequenceT, typename RealT=double >
RealT intersect_and_gather( const IntegerSequence< SequenceT >& a, const IntegerSequence< SequenceT >& b, const std::vector<RealT>& weights );

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_OPERATIONS_FORWARD_H
