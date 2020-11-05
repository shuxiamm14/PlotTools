#include "iostream"
#include "TString.h"
#include <map>

class BelongRegion
{
public:
	BelongRegion();
	~BelongRegion();
	std::map<TString, std::vector<TString>> m_region_map;
	std::vector<TString> m_all_region;
	std::vector<TString> m_enabled_region;
	void add(TString region);
	bool have(TString keyword);
	std::vector<TString> all();
	void clear();
	bool isEmpty();
	bool isCategory(TString category);
	bool isEnabled(TString region);
	void enable(TString region);
	void enable(std::vector<TString> regions);
};