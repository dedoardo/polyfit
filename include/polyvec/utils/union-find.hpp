#pragma once

#include <vector>

namespace polyvec
{
	struct EmptyData
	{ };

	template <typename EntryData = EmptyData>
	class UnionFind
	{
		struct UFEntry : public EntryData
		{
			int parent;			
		};

	public:
		UnionFind(size_t size)
		{
			uf.resize(size);
			for (int i = 0; i < uf.size(); ++i)
				uf[i].parent = i;
		}

		int getRepresentative(int i) 
		{
			if (uf[i].parent != i)
				uf[i].parent = getRepresentative(uf[i].parent);
			return uf[i].parent;
		}

		void merge(int i, int j)
		{
			uf[getRepresentative(i)].parent = getRepresentative(j);
		}

		size_t size() const { return uf.size(); }

		UFEntry& operator[](int i) { return uf[i]; }
		const UFEntry& operator[](int i) const { return uf[i]; }

	private:
		std::vector<UFEntry> uf;
	};
}