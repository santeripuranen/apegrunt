/** @file IndexVector_forward.h

	Copyright (c) 2016 Santeri Puranen. All rights reserved.

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

#ifndef APEGRUNT_INDEXVECTOR_FORWARD_H
#define APEGRUNT_INDEXVECTOR_FORWARD_H

#include <memory> // for std::shared_ptr and std::make_shared

namespace apegrunt {

// forward declarations
template< typename IndexT >
class IndexVector;

template< typename IndexT >
using IndexVector_ptr = std::shared_ptr< IndexVector< IndexT > >;

template< typename IndexVectorT, typename... Args >
IndexVector_ptr< typename IndexVectorT::value_type > make_IndexVector_ptr( Args&... args )
{
	return std::make_shared< IndexVectorT >( args... );
}

} // namespace apegrunt

#endif // APEGRUNT_INDEXVECTOR_FORWARD_H

