#include "region.h"

BelongRegion::BelongRegion(){};

BelongRegion::~BelongRegion(){};

void BelongRegion::add(TString region){
	if(isEnabled(region)) m_all_region.push_back(region);
}

bool BelongRegion::have(TString keyword){
	for(auto reg: m_all_region){
		if(reg.Contains(keyword)) return true;
	}
	return false;
}

std::vector<TString> BelongRegion::all(){
	return m_all_region;
}

void BelongRegion::clear(){
	return m_all_region.clear();
}

bool BelongRegion::isCategory(TString category){
	if(m_region_map.find(category) == m_region_map.end()){
		printf("BelongRegion::isCategory() : ERROR : category %s not found in the region map.\n", category.Data());
	}
	for(auto reg: m_all_region)
		for(auto mapreg: m_region_map.at(category))
			if(reg == mapreg)
				return true;
	return false;
}

bool BelongRegion::isEmpty(){
	return m_all_region.size()==0;
}

void BelongRegion::enable(TString region){
	if(!isEnabled(region)) m_enabled_region.push_back(region);
}

void BelongRegion::enable(std::vector<TString> regions){
	for(auto region: regions)
		enable(region);
}

bool BelongRegion::isEnabled(TString region){
	for(auto reg: m_enabled_region){
		if(reg == region) return true;
	}
	return false;
}