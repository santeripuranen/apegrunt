/** @file IntegerSequence_forward.h

	Copyright (c) 2016-2020 Santeri Puranen. All rights reserved.

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

#ifndef APEGRUNT_INTEGERSEQUENCE_FORWARD_H
#define APEGRUNT_INTEGERSEQUENCE_FORWARD_H

#include <memory> // for std::shared_ptr and std::make_shared

namespace apegrunt {

// forward declarations
template< typename SequenceT >
class IntegerSequence;

template< typename SequenceT >
using IntegerSequence_ptr = std::shared_ptr< IntegerSequence< SequenceT > >;

template< typename IntegerSequenceT, typename... Args >
IntegerSequence_ptr< typename IntegerSequenceT::value_type > make_IntegerSequence_ptr( Args&... args )
{
	return std::make_shared< IntegerSequenceT >( args... );
}

template< typename T > struct hash;
template< typename SequenceT > struct hash< IntegerSequence< SequenceT > >;

template< typename T > struct equal_to;
template< typename SequenceT > struct equal_to< IntegerSequence< SequenceT > >;

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_FORWARD_H
