#include "misc/Stopwatch.hpp"

namespace stopwatch
{

std::ostream& operator<< ( std::ostream& os, const stopwatch& timer )
{
	return timer(os);
}

} // namespace stopwatch
