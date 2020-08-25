/*
	This file is part of NSEssentials.

	Use of this source code is granted via a BSD-style license, which can be found
	in License.txt in the repository root.

	@author Nico Schertler
*/

#pragma once

#include <iterator>

namespace nse
{
	namespace util
	{
		template <typename Iterator>
		class IteratorRange
		{
		public:

			typedef Iterator iterator;

			IteratorRange()
			{ }

			IteratorRange(Iterator begin, Iterator end)
				: _begin(begin), _end(end)
			{ }

			Iterator begin() const { return _begin; }
			Iterator end() const { return _end; }

			//Size of the range, only available for random access iterators
			template <typename IIterator = Iterator>
			typename std::enable_if<std::is_base_of<std::random_access_iterator_tag, typename std::iterator_traits<IIterator>::iterator_category>::value,
				typename std::iterator_traits<IIterator>::difference_type>::type
				size() const { return _end - _begin; }

			//Random access operator, only available for random access iterators
			template <typename IIterator = Iterator>
			typename std::enable_if<std::is_base_of<std::random_access_iterator_tag, typename std::iterator_traits<IIterator>::iterator_category>::value,
				typename std::iterator_traits<IIterator>::reference>::type
				operator[](int i) const { return *(_begin + i); }

			//Returns a new range in the reverse direction of this range, only available for bidirectional iterators
			template <typename IIterator = Iterator>
			typename std::enable_if<std::is_base_of<std::bidirectional_iterator_tag, typename std::iterator_traits<IIterator>::iterator_category>::value,
				IteratorRange<std::reverse_iterator<IIterator>>>::type
				Reverse() const
			{
				return IteratorRange<std::reverse_iterator<IIterator>>(std::reverse_iterator<IIterator>(_end), std::reverse_iterator<IIterator>(_begin));
			}

		protected:
			Iterator _begin, _end;
		};

		template <typename Iterator>
		IteratorRange<Iterator> MakeIteratorRange(Iterator begin, Iterator end)
		{
			return IteratorRange<Iterator>(begin, end);
		}
	}
}